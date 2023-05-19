#include <mitsuba/core/profiler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _sampler-independent:

Independent sampler (:monosp:`independent`)
-------------------------------------------

.. pluginparameters::

 * - sample_count
   - |int|
   - Number of samples per pixel (Default: 4)

 * - seed
   - |int|
   - Seed offset (Default: 0)

The independent sampler produces a stream of independent and uniformly
distributed pseudorandom numbers. Internally, it relies on the
`PCG32 random number generator <https://www.pcg-random.org/>`_
by Melissa O’Neill.

This is the most basic sample generator; because no precautions are taken to avoid
sample clumping, images produced using this plugin will usually take longer to converge.
Looking at the figures below where samples are projected onto a 2D unit square, we see that there
are both regions that don't receive many samples (i.e. we don't know much about the behavior of
the function there), and regions where many samples are very close together (which likely have very
similar values), which will result in higher variance in the rendered image.

This sampler is initialized using a deterministic procedure, which means that subsequent runs
of Mitsuba should create the same image. In practice, when rendering with multiple threads
and/or machines, this is not true anymore, since the ordering of samples is influenced by the
operating system scheduler. Although these should be absolutely negligible, with relative errors
on the order of the machine epsilon (:math:`6\cdot 10^{-8}`) in single precision.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/sampler/independent_1024_samples.svg
   :caption: 1024 samples projected onto the first two dimensions.
.. subfigure:: ../../resources/data/docs/images/sampler/independent_64_samples_and_proj.svg
   :caption: 64 samples projected onto the first two dimensions and their
             projection on both 1D axis (top and right plot).
.. subfigend::
   :label: fig-independent-pattern

.. tabs::
    .. code-tab:: xml
        :name: independent-sampler

        <sampler type="independent">
            <integer name="sample_count" value="64"/>
        </sampler>

    .. code-tab:: python

        'type': 'independent',
        'sample_count': '64'

 */

template <typename Float, typename Spectrum>
class CorrelatedSampler final : public PCG32Sampler<Float, Spectrum> {
public:
    MI_IMPORT_BASE(PCG32Sampler, m_sample_count, m_base_seed, m_rng, seed,
                   seeded, m_samples_per_wavefront, m_wavefront_size,
                   schedule_state, m_dimension_index, current_sample_index,
                   compute_per_sequence_seed)
    MI_IMPORT_TYPES()

    CorrelatedSampler(const Properties &props) : Base(props) {
        m_time_correlate_number = props.get<int>("time_correlate_number", 2);
        m_path_correlate_number = props.get<int>("path_correlate_number", m_time_correlate_number);

        
        m_use_stratified_sampling_for_each_interval = props.get<bool>("use_stratified_sampling_for_each_interval", true);
        m_antithetic_shift = props.get<ScalarFloat>("antithetic_shift", 0.0);

        std::cout << "m_time_correlate_number" << m_time_correlate_number << std::endl;
        std::cout << "m_path_correlate_number" << m_path_correlate_number << std::endl;
    }

    ref<Sampler<Float, Spectrum>> fork() override {
        CorrelatedSampler *sampler = new CorrelatedSampler(Properties());
        sampler->m_sample_count = m_sample_count;
        sampler->m_base_seed = m_base_seed;
        sampler->m_time_correlate_number = m_time_correlate_number;
        sampler->m_path_correlate_number = m_path_correlate_number;
        sampler->m_use_stratified_sampling_for_each_interval = m_use_stratified_sampling_for_each_interval;
        sampler->m_antithetic_shift = m_antithetic_shift;
        return sampler;
    }

    ref<Sampler<Float, Spectrum>> clone() override {
        return new CorrelatedSampler(*this);
    }

    void seed(uint32_t seed, uint32_t wavefront_size) override {
        Base::seed(seed, wavefront_size);
        m_permutation_seed = compute_per_sequence_seed(seed);

        uint32_t seed_value = m_base_seed + seed;

        if constexpr (dr::is_array_v<Float>) {
            UInt32 idx = dr::arange<UInt32>(m_wavefront_size),
                tmp = dr::opaque<UInt32>(seed_value);

            UInt32 time_idx = idx / m_time_correlate_number;
            UInt32 path_idx = idx / m_path_correlate_number;
            
            /* Scramble seed and stream index using the Tiny Encryption Algorithm.
            Just providing a linearly increasing sequence of integers as streams
            does not produce a sufficiently statistically independent set of RNGs */
            auto [v0, v1] = sample_tea_32(tmp + 1, time_idx);
            auto [v2, v3] = sample_tea_32(tmp + 2, path_idx);

            m_rng_time.seed(1, v0, v1);
            m_rng_path.seed(1, v2, v3);

        } else {
            m_rng_time.seed(1, seed_value + 1, PCG32_DEFAULT_STREAM);
            m_rng_path.seed(1, seed_value + 2, PCG32_DEFAULT_STREAM);
        }
    }

    void schedule_state() override {
        Base::schedule_state();
        dr::schedule(m_rng_time.inc, m_rng_time.state);
        dr::schedule(m_rng_path.inc, m_rng_path.state);
        dr::schedule(m_permutation_seed);
    }

    void loop_put(dr::Loop<Mask> &loop) override {
        Base::loop_put(loop);
        loop.put(m_rng_time.state);
        loop.put(m_rng_path.state);
    }

    Float next_1d(Mask active = true) override {
        Assert(seeded());
        // return next_1d_time(active, TIME_SAMPLING_STRATIFIED, 0.0);
        
        return m_rng.template next_float<Float>(active);
    }

    Point2f next_2d(Mask active = true) override {
        Float f1 = next_1d(active),
              f2 = next_1d(active);
        return Point2f(f1, f2);
    }

    Float next_1d_time(Mask active = true, int strategy = TIME_SAMPLING_UNIFORM, ScalarFloat antithetic_shift = 0.0) override {
        Assert(seeded());
        // (1) Uniform sampling -> Just use m_rng
        if(strategy == TIME_SAMPLING_UNIFORM){
            return m_rng.template next_float<Float>(active);
        }
        
        UInt32 sample_indices = current_sample_index();
        // UInt32 perm_seed = m_permutation_seed + m_dimension_index++;

        // // Shuffle the samples order
        // UInt32 p = permute_kensler(sample_indices / 2, m_sample_count, perm_seed, active);

        // // Add a random perturbation
        // Float j = m_rng.template next_float<Float>(active);

        // return (p + j) * dr::rcp(ScalarFloat(m_sample_count));

        // Assert(seeded());
        // int n_stratum = m_sample_count / 2;

        // UInt32 sample_indices = current_sample_index();
        // UInt32 perm_seed = m_permutation_seed + m_dimension_index++;
        // UInt32 sample_indices_temp = permute_kensler(sample_indices / 2, n_stratum, perm_seed, active);
        // UInt32 p = sample_indices_temp;
        
        // Float r = m_rng.template next_float<Float>(active);
        // Float j = (p + r) * dr::rcp(ScalarFloat(n_stratum));

        // UInt32 p2 = sample_indices % 2;//permute_kensler(sample_indices, m_sample_count, perm_seed, active);
        // // r = m_rng.template next_float<Float>(active);
        // return (p2 + j) * 0.5;
        
        Float r = 0.0;

        if ((strategy == TIME_SAMPLING_UNIFORM) || (strategy == TIME_SAMPLING_STRATIFIED)){
            r = m_rng.template next_float<Float>(active);
        } else {
            r = m_rng_time.template next_float<Float>(active);
        }
        //m_use_stratified_sampling_for_each_interval = ;

        if(m_use_stratified_sampling_for_each_interval){
            
            int n_stratum = m_sample_count / m_time_correlate_number;

            if(strategy == TIME_SAMPLING_STRATIFIED){
                UInt32 perm_seed = m_permutation_seed + m_dimension_index++;
                UInt32 p1 = permute_kensler(sample_indices / m_time_correlate_number, n_stratum, perm_seed, active);
                
                perm_seed = m_permutation_seed + m_dimension_index++;
                UInt32 p2 = permute_kensler(sample_indices / m_time_correlate_number, n_stratum, perm_seed, active);
                
                UInt32 p = dr::select(dr::neq(sample_indices % m_time_correlate_number, 0), p1, p2);

                // p = p / m_time_correlate_number;
                r = (p + r) / n_stratum;
            } else {
                UInt32 p = sample_indices / m_time_correlate_number;
                r = (p + r) / n_stratum;
            }
        }

        
        // (2) Stratified
        if(strategy == TIME_SAMPLING_STRATIFIED){
            // return r;
            // UInt32 perm_seed = m_permutation_seed + m_dimension_index++;

            // // Shuffle the samples order
            UInt32 p = sample_indices % m_time_correlate_number;//permute_kensler(sample_indices, m_sample_count, perm_seed, active);
            // // r = m_rng.template next_float<Float>(active);
            return (p + r) * dr::rcp(ScalarFloat(m_time_correlate_number));
        }
        else if(strategy == TIME_SAMPLING_ANTITHETIC){
            Assert(m_time_correlate_number == 2);
            // use stratified sampling
            // if(m_use_stratified_sampling_for_each_interval){
            //     int n_stratum = m_sample_count / m_time_correlate_number;
            //     UInt32 p = sample_indices / m_time_correlate_number;
            //     Float j = m_rng_time.template next_float<Float>(active);
            //     Float r = (p + j) / n_stratum;
            //     Float r2 = r + antithetic_shift;
            //     UInt32 remainder = sample_indices % m_time_correlate_number;
            //     return dr::select(dr::neq(remainder, 1), r, r2);
            // }
            // Float r = m_rng_time.template next_float<Float>(active);
            Float r2 = r + antithetic_shift;
            UInt32 remainder = sample_indices % m_time_correlate_number;
            return dr::select(dr::neq(remainder, 1), r, r2);
        }
        else if(strategy == TIME_SAMPLING_ANTITHETIC_MIRROR){
            Assert(m_time_correlate_number == 2);
            // use stratified sampling
            // if(m_use_stratified_sampling_for_each_interval){
            //     int n_stratum = m_sample_count / m_time_correlate_number;
            //     UInt32 p = sample_indices / m_time_correlate_number;
            //     Float j = m_rng_time.template next_float<Float>(active);
            //     Float r = (p + j) / n_stratum;
            //     Float r2 = 1.0 - r + antithetic_shift;
            //     UInt32 remainder = sample_indices % m_time_correlate_number;
            //     return dr::select(dr::neq(remainder, 1), r, r2);
            // }
            // Float r = m_rng_time.template next_float<Float>(active);
            Float r2 = 1.0 - r + antithetic_shift;
            UInt32 remainder = sample_indices % m_time_correlate_number;
            return dr::select(dr::neq(remainder, 1), r, r2);
        }
        else if(strategy == TIME_SAMPLING_PERIODIC){
            // use stratified sampling
            // if(m_use_stratified_sampling_for_each_interval){
            //     int n_stratum = m_sample_count / m_time_correlate_number;
            //     UInt32 p = sample_indices / m_time_correlate_number;
            //     Float j = m_rng_time.template next_float<Float>(active);
            //     Float r = (p + j) / n_stratum;
                
            //     UInt32 remainder = sample_indices % m_time_correlate_number;
            //     UInt32 sample_indices = current_sample_index();
            //     return r + (Float(remainder)) / Float(m_time_correlate_number);
            // }

            // Float r = m_rng_time.template next_float<Float>(active);
            UInt32 remainder = sample_indices % m_time_correlate_number;
            return r + Float(remainder) / Float(m_time_correlate_number);
        }
        
        // Float r = m_rng.template next_float<Float>(active);
        return r;
    }


    Float next_1d_correlate(Mask active = true, Bool correlate = false) override {
        Assert(seeded());
        Float r1 = m_rng_path.template next_float<Float>(active);
        Float r2 = m_rng.template next_float<Float>(active);
        return dr::select(correlate, r1, r2);

        // if(correlate){
        //     Float r = m_rng_path.template next_float<Float>(active);
        //     return r;
        // }
        // Float r = m_rng.template next_float<Float>(active);
        // return r;
    }

    Point2f next_2d_correlate(Mask active = true, Bool correlate = false) override {
        Float f1 = next_1d_correlate(active, correlate),
              f2 = next_1d_correlate(active, correlate);
        return Point2f(f1, f2);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "CorrelatedSampler[" << std::endl
            << "  base_seed = " << m_base_seed << std::endl
            << "  sample_count = " << m_sample_count << std::endl
            << "  samples_per_wavefront = " << m_samples_per_wavefront << std::endl
            << "  wavefront_size = " << m_wavefront_size << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

protected:
    mitsuba::PCG32<UInt32> m_rng_time;
    int m_time_correlate_number;
    mitsuba::PCG32<UInt32> m_rng_path;
    int m_path_correlate_number;
    ScalarFloat m_antithetic_shift;
    bool m_use_stratified_sampling_for_each_interval;

    /// Per-sequence permutation seed
    UInt32 m_permutation_seed;
    
private:
    CorrelatedSampler(const CorrelatedSampler &sampler) : Base(sampler) {}
};

MI_IMPLEMENT_CLASS_VARIANT(CorrelatedSampler, Sampler)
MI_EXPORT_PLUGIN(CorrelatedSampler, "Independent Sampler");
NAMESPACE_END(mitsuba)
