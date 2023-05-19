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
                   schedule_state, m_dimension_index, current_sample_index)
    MI_IMPORT_TYPES()

    CorrelatedSampler(const Properties &props) : Base(props) { }

    ref<Sampler<Float, Spectrum>> fork() override {
        CorrelatedSampler *sampler = new CorrelatedSampler(Properties());
        sampler->m_sample_count = m_sample_count;
        sampler->m_base_seed = m_base_seed;
        return sampler;
    }

    ref<Sampler<Float, Spectrum>> clone() override {
        return new CorrelatedSampler(*this);
    }

    Float next_1d(Mask active = true) override {
        Assert(seeded());
        return m_rng.template next_float<Float>(active);
    }

    Point2f next_2d(Mask active = true) override {
        Float f1 = next_1d(active),
              f2 = next_1d(active);
        return Point2f(f1, f2);
    }

    Float next_1d_time(Mask active = true, int strategy = TIME_SAMPLING_UNIFORM, int time_intervals = 2) override {
        Assert(seeded());

        if(strategy == TIME_SAMPLING_UNIFORM){
            Float r = m_rng.template next_float<Float>(active);
            return r;
        }
        else if(strategy == TIME_SAMPLING_STRATIFIED){
            UInt32 p = current_sample_index();
            p = p % time_intervals;
            Float j = m_rng.template next_float<Float>(active);
            Float r = (p + j) / time_intervals;
            return r;
        }
        else if(strategy == TIME_SAMPLING_REGULAR){
            UInt32 p = current_sample_index();
            p = p % time_intervals;
            Float r = (p + .5f) / time_intervals;
            return r;
        }
        else if(strategy == TIME_SAMPLING_ANTITHETIC || strategy == TIME_SAMPLING_PERIODIC){
            
            Float r = m_rng.template next_float<Float>(active);
            
            if constexpr(dr::is_array_v<Float>){    
                UInt32 sample_indices = current_sample_index();
                return dr::gather<Float>(r, sample_indices / time_intervals);
            }

            return r;

            // UInt32 p = sample_indices / time_intervals;
            // UInt32 remainder = sample_indices % time_intervals;

            // UInt32 rand_uint = (sample_tea_32(p, m_dimension_index++).first & 0x00FFFFFF);
            // Float rand_float = Float(rand_uint) / Float(0x01000000);

            // Float r = (remainder + rand_float) / time_intervals;
            // return r;
        }
        else if(strategy == TIME_SAMPLING_ANTITHETIC_MIRROR){
            Float r = m_rng.template next_float<Float>(active);
            UInt32 sample_indices = current_sample_index();

            UInt32 p = sample_indices / time_intervals;
            UInt32 remainder = sample_indices %  time_intervals;

            UInt32 rand_uint = (sample_tea_32(p, m_dimension_index++).first & 0x00FFFFFF);
            Float rand_float = Float(rand_uint) / Float(0x01000000);
            Float rand_float_minus = 1.0 - rand_float;
            
            return dr::select(dr::neq(remainder, 1), rand_float, rand_float_minus);
        }
        Float r = m_rng.template next_float<Float>(active);
        return r;
    }

    Float next_1d_correlate(Mask active = true, bool correlate = false, int time_intervals = 2) override {
        Assert(seeded());

        Float r = m_rng.template next_float<Float>(active);
        if constexpr(dr::is_array_v<Float>){
            UInt32 sample_indices = current_sample_index();
            r = dr::gather<Float>(r, sample_indices / time_intervals);
        }
        return r;
    }

    Point2f next_2d_correlate(Mask active = true, bool correlate = false, int time_intervals = 2) override {
        Assert(seeded());
        
        Float f1 = m_rng.template next_float<Float>(active);
        Float f2 = m_rng.template next_float<Float>(active);
        if constexpr(dr::is_array_v<Float>){
            UInt32 sample_indices = current_sample_index();
            f1 = dr::gather<Float>(f1, sample_indices / time_intervals);
            f2 = dr::gather<Float>(f2, sample_indices / time_intervals);
        }
        if(correlate){
            
            //UInt32 sample_indices = current_sample_index();
            //UInt32 p = sample_indices / time_intervals;
            //auto r12 = sample_tea_32(p, m_dimension_index++);
            //f1 = Float(r12.first & 0x00FFFFFF) / Float(0x01000000);
            //f2 = Float(r12.second & 0x00FFFFFF) / Float(0x01000000);
        }
        
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

private:
    CorrelatedSampler(const CorrelatedSampler &sampler) : Base(sampler) {}
};

MI_IMPLEMENT_CLASS_VARIANT(CorrelatedSampler, Sampler)
MI_EXPORT_PLUGIN(CorrelatedSampler, "Independent Sampler");
NAMESPACE_END(mitsuba)
