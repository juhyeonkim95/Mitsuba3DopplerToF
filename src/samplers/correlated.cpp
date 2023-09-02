#include <mitsuba/core/profiler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

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

    Float next_1d_time(Mask active = true, ETimeSampling strategy = ETimeSampling::TIME_SAMPLING_UNIFORM, ScalarFloat antithetic_shift = 0.0) override {
        Assert(seeded());
        // (1) Uniform sampling -> Just use m_rng
        if(strategy == TIME_SAMPLING_UNIFORM){
            return m_rng.template next_float<Float>(active);
        }
        
        UInt32 sample_indices = current_sample_index();
        Float r = 0.0;

        // Prepare random numbers
        if (strategy == ETimeSampling::TIME_SAMPLING_STRATIFIED){
            r = m_rng.template next_float<Float>(active);
        } else {
            r = m_rng_time.template next_float<Float>(active);
        }

        if(m_use_stratified_sampling_for_each_interval){
            int n_stratum = m_sample_count / m_time_correlate_number;
            if(strategy == ETimeSampling::TIME_SAMPLING_STRATIFIED){
                // use permutation
                UInt32 perm_seed = m_permutation_seed + m_dimension_index++;
                UInt32 p1 = permute_kensler(sample_indices / m_time_correlate_number, n_stratum, perm_seed, active);
                
                perm_seed = m_permutation_seed + m_dimension_index++;
                UInt32 p2 = permute_kensler(sample_indices / m_time_correlate_number, n_stratum, perm_seed, active);
                
                UInt32 p = dr::select(dr::neq(sample_indices % m_time_correlate_number, 0), p1, p2);
                r = (p + r) / n_stratum;
            } else {
                UInt32 p = sample_indices / m_time_correlate_number;
                r = (p + r) / n_stratum;
            }
        }
        
        if(strategy == ETimeSampling::TIME_SAMPLING_STRATIFIED){
            UInt32 p = sample_indices % m_time_correlate_number;
            return (p + r) * dr::rcp(ScalarFloat(m_time_correlate_number));
        }
        else if(strategy == ETimeSampling::TIME_SAMPLING_ANTITHETIC){
            if(m_time_correlate_number == 2){
                Float r2 = r + antithetic_shift;
                UInt32 remainder = sample_indices % m_time_correlate_number;
                return dr::select(dr::neq(remainder, 1), r, r2);
            } else {
                UInt32 remainder = sample_indices % m_time_correlate_number;
                return r + Float(remainder) / Float(m_time_correlate_number);
            }
        }
        else if(strategy == ETimeSampling::TIME_SAMPLING_ANTITHETIC_MIRROR){
            Assert(m_time_correlate_number == 2);
            Float r2 = 1.0 - r + antithetic_shift;
            UInt32 remainder = sample_indices % m_time_correlate_number;
            return dr::select(dr::neq(remainder, 1), r, r2);
        }
        else if(strategy == ETimeSampling::TIME_SAMPLING_PERIODIC){
            UInt32 remainder = sample_indices % m_time_correlate_number;
            return r + Float(remainder) / Float(m_time_correlate_number);
        }
    
        return r;
    }


    Float next_1d_correlate(Mask active = true, Bool correlate = false) override {
        Assert(seeded());
        Float r1 = m_rng_path.template next_float<Float>(active);
        Float r2 = m_rng.template next_float<Float>(active);
        return dr::select(correlate, r1, r2);
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
    mitsuba::PCG32<UInt32> m_rng_time;  // correlated time sampler
    int m_time_correlate_number;
    mitsuba::PCG32<UInt32> m_rng_path;  // correlated path sampler
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
