#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/filesystem.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/integrator.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Abstract integrator base class, which does not make any assumptions
 * with regards to how radiance is computed.
 *
 * In Mitsuba, the different rendering techniques are collectively referred to
 * as \a integrators, since they perform integration over a high-dimensional
 * space. Each integrator represents a specific approach for solving the light
 * transport equation---usually favored in certain scenarios, but at the same
 * time affected by its own set of intrinsic limitations. Therefore, it is
 * important to carefully select an integrator based on user-specified accuracy
 * requirements and properties of the scene to be rendered.
 *
 * This is the base class of all integrators; it does not make any assumptions
 * on how radiance is computed, which allows for many different kinds of
 * implementations.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB TransientIntegrator : public Integrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Integrator, should_stop, aov_names,
                    m_stop, m_timeout, m_render_timer, m_hide_emitters)
    MI_IMPORT_TYPES(Scene, Sensor)


    MI_DECLARE_CLASS()
protected:
    /// Create an integrator
    TransientIntegrator(const Properties & props);

    /// Virtual destructor
    virtual ~TransientIntegrator() { }
};

/** \brief Abstract integrator that performs Monte Carlo sampling starting from
 * the sensor
 *
 * Subclasses of this interface must implement the \ref sample() method, which
 * performs Monte Carlo integration to return an unbiased statistical estimate
 * of the radiance value along a given ray.
 *
 * The \ref render() method then repeatedly invokes this estimator to compute
 * all pixels of the image.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB TransientSamplingIntegrator : public TransientIntegrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(TransientIntegrator, should_stop, aov_names,
                    m_stop, m_timeout, m_render_timer, m_hide_emitters)
    MI_IMPORT_TYPES(Scene, Sensor, Film, ImageBlock, Medium, Sampler)

    /**
     * \brief Sample the incident radiance along a ray.
     *
     * \param scene
     *    The underlying scene in which the radiance function should be sampled
     *
     * \param sampler
     *    A source of (pseudo-/quasi-) random numbers
     *
     * \param ray
     *    A ray, optionally with differentials
     *
     * \param medium
     *    If the ray is inside a medium, this parameter holds a pointer to that
     *    medium
     *
     * \param aov
     *    Integrators may return one or more arbitrary output variables (AOVs)
     *    via this parameter. If \c nullptr is provided to this argument, no
     *    AOVs should be returned. Otherwise, the caller guarantees that space
     *    for at least <tt>aov_names().size()</tt> entries has been allocated.
     *
     * \param active
     *    A mask that indicates which SIMD lanes are active
     *
     * \return
     *    A pair containing a spectrum and a mask specifying whether a surface
     *    or medium interaction was sampled. False mask entries indicate that
     *    the ray "escaped" the scene, in which case the the returned spectrum
     *    contains the contribution of environment maps, if present. The mask
     *    can be used to estimate a suitable alpha channel of a rendered image.
     *
     * \remark
     *    In the Python bindings, this function returns the \c aov output
     *    argument as an additional return value. In other words:
     *    <pre>
     *        (spec, mask, aov) = integrator.sample(scene, sampler, ray, medium, active)
     *    </pre>
     */
    virtual std::pair<Spectrum, Mask> sample(const Scene *scene,
                                             Sampler *sampler,
                                             const RayDifferential3f &ray,
                                             const Medium *medium = nullptr,
                                             Float *aovs = nullptr,
                                             Mask active = true) const;

    // =========================================================================
    //! @{ \name Integrator interface implementation
    // =========================================================================

    TensorXf render(Scene *scene,
                    Sensor *sensor,
                    uint32_t seed = 0,
                    uint32_t spp = 0,
                    bool develop = true,
                    bool evaluate = true) override;

    //! @}
    // =========================================================================

    MI_DECLARE_CLASS()
protected:
    TransientSamplingIntegrator(const Properties &props);
    virtual ~TransientSamplingIntegrator();

    virtual void render_block(const Scene *scene,
                              const Sensor *sensor,
                              Sampler *sampler,
                              ImageBlock *block,
                              Float *aovs,
                              uint32_t sample_count,
                              uint32_t seed,
                              uint32_t block_id,
                              uint32_t block_size) const;

    void render_sample(const Scene *scene,
                       const Sensor *sensor,
                       Sampler *sampler,
                       ImageBlock *block,
                       Float *aovs,
                       const Vector2f &pos,
                       ScalarFloat diff_scale_factor,
                       Mask active = true) const;

protected:

    /// Size of (square) image blocks to render in parallel (in scalar mode)
    uint32_t m_block_size;

    /**
     * \brief Number of samples to compute for each pass over the image blocks.
     *
     * Must be a multiple of the total sample count per pixel.
     * If set to (uint32_t) -1, all the work is done in a single pass (default).
     */
    uint32_t m_samples_per_pass;
};

/** \brief Abstract integrator that performs *recursive* Monte Carlo sampling
 * starting from the sensor
 *
 * This class is almost identical to \ref TransientSamplingIntegrator. It stores two
 * additional fields that are helpful for recursive Monte Carlo techniques:
 * the maximum path depth, and the depth at which the Russian Roulette path
 * termination technique should start to become active.
 */
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB TransientMonteCarloIntegrator
    : public TransientSamplingIntegrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(TransientSamplingIntegrator)

protected:
    /// Create an integrator
    TransientMonteCarloIntegrator(const Properties &props);

    /// Virtual destructor
    virtual ~TransientMonteCarloIntegrator();

    MI_DECLARE_CLASS()
protected:
    uint32_t m_max_depth;
    uint32_t m_rr_depth;
};

MI_EXTERN_CLASS(TransientIntegrator)
MI_EXTERN_CLASS(TransientSamplingIntegrator)
MI_EXTERN_CLASS(TransientMonteCarloIntegrator)
NAMESPACE_END(mitsuba)