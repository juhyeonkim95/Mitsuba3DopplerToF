#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/random.h>
#include <drjit/loop.h>

// #ifndef TIME_SAMPLING_UNIFORM
// #define TIME_SAMPLING_UNIFORM 0
// #define TIME_SAMPLING_STRATIFIED 1
// #define TIME_SAMPLING_ANTITHETIC 2
// #define TIME_SAMPLING_ANTITHETIC_MIRROR 3
// #define TIME_SAMPLING_PERIODIC 4
// #define TIME_SAMPLING_REGULAR 5
// #define SPATIAL_CORRELATION_NONE 0
// #define SPATIAL_CORRELATION_PIXEL 1
// #define SPATIAL_CORRELATION_SAMPLER 2
// #endif



NAMESPACE_BEGIN(mitsuba)

enum ETimeSampling{
    TIME_SAMPLING_UNIFORM = 0,
    TIME_SAMPLING_STRATIFIED = 1,
    TIME_SAMPLING_ANTITHETIC = 2,
    TIME_SAMPLING_ANTITHETIC_MIRROR = 3,
    TIME_SAMPLING_PERIODIC = 4,
    TIME_SAMPLING_REGULAR = 5
};

/**
 * \brief Base class of all sample generators.
 *
 * A \a sampler provides a convenient abstraction around methods that generate
 * uniform pseudo- or quasi-random points within a conceptual
 * infinite-dimensional unit hypercube \f$[0,1]^\infty$\f. This involves two
 * main operations: by querying successive component values of such an
 * infinite-dimensional point (\ref next_1d(), \ref next_2d()), or by
 * discarding the current point and generating another one (\ref advance()).
 *
 * Scalar and vectorized rendering algorithms interact with the sampler
 * interface in a slightly different way:
 *
 * Scalar rendering algorithm:
 *
 *   1. The rendering algorithm first invokes \ref seed() to initialize the
 *      sampler state.
 *
 *   2. The first pixel sample can now be computed, after which \ref advance()
 *      needs to be invoked. This repeats until all pixel samples have been
 *      generated. Note that some implementations need to be configured for a
 *      certain number of pixel samples, and exceeding these will lead to an
 *      exception being thrown.
 *
 *   3. While computing a pixel sample, the rendering algorithm usually
 *      requests 1D or 2D component blocks using the \ref next_1d() and
 *      \ref next_2d() functions before moving on to the next sample.
 *
 * A vectorized rendering algorithm effectively queries multiple sample
 * generators that advance in parallel. This involves the following steps:
 *
 *   1. The rendering algorithm invokes \ref set_samples_per_wavefront()
 *      if each rendering step is split into multiple passes (in which
 *      case fewer samples should be returned per \ref sample_1d()
 *      or \ref sample_2d() call).
 *
 *   2. The rendering algorithm then invokes \ref seed() to initialize the
 *      sampler state, and to inform the sampler of the wavefront size,
 *      i.e., how many sampler evaluations should be performed in parallel,
 *      accounting for all passes. The initialization ensures that the set of
 *      parallel samplers is mutually statistically independent (in a
 *      pseudo/quasi-random sense).
 *
 *   3. \ref advance() can be used to advance to the next point.
 *
 *   4. As in the scalar approach, the rendering algorithm can request batches
 *      of (pseudo-) random numbers using the \ref next_1d() and \ref next_2d()
 *      functions.
 *
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB Sampler : public Object {
public:
    MI_IMPORT_TYPES()

    /**
     * \brief Create a fork of this sampler.
     *
     * A subsequent call to \c seed is necessary to properly initialize
     * the internal state of the sampler.
     *
     * May throw an exception if not supported.
     */
    virtual ref<Sampler> fork() = 0;

    /**
     * \brief Create a clone of this sampler.
     *
     * Subsequent calls to the cloned sampler will produce the same
     * random numbers as the original sampler.
     *
     * \remark This method relies on the overload of the copy constructor.
     *
     * May throw an exception if not supported.
     */
    virtual ref<Sampler> clone() = 0;

    /**
     * \brief Deterministically seed the underlying RNG, if applicable.
     *
     * In the context of wavefront ray tracing & dynamic arrays, this function
     * must be called with \c wavefront_size matching the size of the wavefront.
     */
    virtual void seed(uint32_t seed,
                      uint32_t wavefront_size = (uint32_t) -1);

    /**
     * \brief Advance to the next sample.
     *
     * A subsequent call to \c next_1d or \c next_2d will access the first
     * 1D or 2D components of this sample.
     */
    virtual void advance();

    /// Retrieve the next component value from the current sample
    virtual Float next_1d_time(Mask active = true, ETimeSampling strategy = ETimeSampling::TIME_SAMPLING_UNIFORM, 
        ScalarFloat antithetic_shift = 0.0, bool use_stratified_sampling_for_each_interval = false) { return next_1d(); }

    /// Retrieve the next component value from the current sample
    virtual Float next_1d(Mask active = true);

    /// Retrieve the next two component values from the current sample
    virtual Point2f next_2d(Mask active = true);

    /// Retrieve the next component value from the current sample
    virtual Float next_1d_correlate(Mask active = true, Bool correlate = false) { return next_1d(active); };

    /// Retrieve the next two component values from the current sample
    virtual Point2f next_2d_correlate(Mask active = true, Bool correlate = false) { return next_2d(active); };

    /// seed for time and path random generator
    // virtual void seed_time_path(uint32_t seed, uint32_t time_correlate_number, uint32_t path_correlate_number);
    
    /// Return the number of samples per pixel
    uint32_t sample_count() const { return m_sample_count; }

    /// Set the number of samples per pixel
    virtual void set_sample_count(uint32_t spp) { m_sample_count = spp; }

    /// Return the size of the wavefront (or 0, if not seeded)
    uint32_t wavefront_size() const { return m_wavefront_size; };

    /// Return whether the sampler was seeded
    bool seeded() const { return m_wavefront_size > 0; }

    /// Set the number of samples per pixel per pass in wavefront modes (default is 1)
    void set_samples_per_wavefront(uint32_t samples_per_wavefront);

    /// dr::schedule() variables that represent the internal sampler state
    virtual void schedule_state();

    /// Register internal state of this sampler with a symbolic loop
    virtual void loop_put(dr::Loop<Mask> &loop);

    MI_DECLARE_CLASS()
protected:
    Sampler(const Properties &props);
    /// Copy state to a new sampler object
    Sampler(const Sampler&);
    virtual ~Sampler();

    /// Generates a array of seeds where the seed values are unique per sequence
    UInt32 compute_per_sequence_seed(uint32_t seed) const;
    /// Return the index of the current sample
    UInt32 current_sample_index() const;

protected:
    /// Base seed value
    uint32_t m_base_seed;
    /// Number of samples per pixel
    uint32_t m_sample_count;
    /// Number of samples per pass in wavefront modes (default is 1)
    uint32_t m_samples_per_wavefront;
    /// Size of the wavefront (or 0, if not seeded)
    uint32_t m_wavefront_size;
    /// Index of the current dimension in the sample
    UInt32 m_dimension_index;
    /// Index of the current sample in the sequence
    UInt32 m_sample_index;
};

/// Interface for sampler plugins based on the PCG32 random number generator
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB PCG32Sampler : public Sampler<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Sampler, m_base_seed, m_wavefront_size)
    MI_IMPORT_TYPES()
    using PCG32 = mitsuba::PCG32<UInt32>;

    virtual void seed(uint32_t seed,
                      uint32_t wavefront_size = (uint32_t) -1) override;
    virtual void schedule_state() override;
    virtual void loop_put(dr::Loop<Mask> &loop) override;

    MI_DECLARE_CLASS()
protected:
    PCG32Sampler(const Properties &props);

    /// Copy state to a new PCG32Sampler object
    PCG32Sampler(const PCG32Sampler &sampler);
protected:
    PCG32 m_rng;
};

MI_EXTERN_CLASS(Sampler)
MI_EXTERN_CLASS(PCG32Sampler)
NAMESPACE_END(mitsuba)
