#include <tuple>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-path:

Path tracer (:monosp:`path`)
----------------------------

.. pluginparameters::

 * - max_depth
   - |int|
   - Specifies the longest path depth in the generated output image (where -1
     corresponds to :math:`\infty`). A value of 1 will only render directly
     visible light sources. 2 will lead to single-bounce (direct-only)
     illumination, and so on. (Default: -1)

 * - rr_depth
   - |int|
   - Specifies the path depth, at which the implementation will begin to use
     the *russian roulette* path termination criterion. For example, if set to
     1, then path generation many randomly cease after encountering directly
     visible surfaces. (Default: 5)

 * - hide_emitters
   - |bool|
   - Hide directly visible emitters. (Default: no, i.e. |false|)

This integrator implements a basic path tracer and is a **good default choice**
when there is no strong reason to prefer another method.

To use the path tracer appropriately, it is instructive to know roughly how
it works: its main operation is to trace many light paths using *random walks*
starting from the sensor. A single random walk is shown below, which entails
casting a ray associated with a pixel in the output image and searching for
the first visible intersection. A new direction is then chosen at the intersection,
and the ray-casting step repeats over and over again (until one of several
stopping criteria applies).

.. image:: ../../resources/data/docs/images/integrator/integrator_path_figure.png
    :width: 95%
    :align: center

At every intersection, the path tracer tries to create a connection to
the light source in an attempt to find a *complete* path along which
light can flow from the emitter to the sensor. This of course only works
when there is no occluding object between the intersection and the emitter.

This directly translates into a category of scenes where a path tracer can be
expected to produce reasonable results: this is the case when the emitters are
easily "accessible" by the contents of the scene. For instance, an interior
scene that is lit by an area light will be considerably harder to render when
this area light is inside a glass enclosure (which effectively counts as an
occluder).

Like the :ref:`direct <integrator-direct>` plugin, the path tracer internally
relies on multiple importance sampling to combine BSDF and emitter samples. The
main difference in comparison to the former plugin is that it considers light
paths of arbitrary length to compute both direct and indirect illumination.

.. note:: This integrator does not handle participating media

.. tabs::
    .. code-tab::  xml
        :name: path-integrator

        <integrator type="path">
            <integer name="max_depth" value="8"/>
        </integrator>

    .. code-tab:: python

        'type': 'path',
        'max_depth': 8

 */

template <typename Float, typename Spectrum>
class DopplerPathIntegratorPotisional : public MonteCarloIntegrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(MonteCarloIntegrator, m_max_depth, m_rr_depth, m_hide_emitters)
    MI_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr)

    DopplerPathIntegratorPotisional(const Properties &props) : Base(props) {
        m_time = props.get<ScalarFloat>("time", 0.0015f);

        m_illumination_modulation_frequency_mhz = props.get<ScalarFloat>("w_g", 30.0f);
        m_illumination_modulation_scale = props.get<ScalarFloat>("g_1", 0.5f);
        m_illumination_modulation_offset = props.get<ScalarFloat>("g_0", 0.5f);

        m_sensor_modulation_frequency_mhz = props.get<ScalarFloat>("w_f", 30.0f);
        m_hetero_frequency = props.get<ScalarFloat>("hetero_frequency", 1.0f);
        m_sensor_modulation_frequency_mhz = m_illumination_modulation_frequency_mhz + 1 / m_time * 1e-6;

        m_sensor_modulation_scale = props.get<ScalarFloat>("f_1", 0.5f);
        m_sensor_modulation_offset = props.get<ScalarFloat>("f_0", 0.5f);
        m_sensor_modulation_phase_offset = props.get<ScalarFloat>("f_phase_offset", 0.0f);

        m_low_frequency_component_only = props.get<bool>("low_frequency_component_only", true);
    }


    Float evalModulationWeight(Float ray_time, Float path_length) const
    {
        //Float w_g = 2 * M_PI * m_illumination_modulation_frequency_mhz * 1e6;
        //Float w_f = 2 * M_PI * m_sensor_modulation_frequency_mhz * 1e6;

        Float w_d = 2 * M_PI / m_time * m_hetero_frequency;

        Float phi = (2 * M_PI * m_illumination_modulation_frequency_mhz) / 300 * path_length;
        Float fg_t = 0.25 * dr::cos(w_d * ray_time + phi);
        return fg_t;

        // if(m_low_frequency_component_only){
        //     Float fg_t = 0.25 * dr::cos(w_d * ray_time + phi);
        //     return fg_t;
        // } 
        
        // Float g_t = 0.5 * dr::cos(w_g * ray_time - phi) + 0.5;
        // Float f_t = dr::cos(w_f * ray_time);
        // return f_t * g_t;
    }

    struct PathTracePart{
        SurfaceInteraction3f si;
        Float path_length;
    };

    std::pair<Spectrum, Bool> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray_,
                                     const Medium * /* medium */,
                                     Float * /* aovs */,
                                     Bool active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        if (unlikely(m_max_depth == 0))
            return { 0.f, false };

        // --------------------- Configure loop state ----------------------

        Ray3f ray                     = Ray3f(ray_);
        Spectrum throughput           = 1.f;
        Spectrum result               = 0.f;
        Float path_length             = 0.f;
        Float eta                     = 1.f;
        UInt32 depth                  = 0;

        // If m_hide_emitters == false, the environment emitter will be visible
        Mask valid_ray                = !m_hide_emitters && dr::neq(scene->environment(), nullptr);

        // Variables caching information from the previous bounce
        Interaction3f prev_si         = dr::zeros<Interaction3f>();
        Float         prev_bsdf_pdf   = 1.f;
        Bool          prev_bsdf_delta = true;
        BSDFContext   bsdf_ctx;
        
        PathTracePart ppart;

        // antithetic sample
        Float ray2_time = ray.time + 0.5 * m_time;
        ray2_time = dr::select(ray2_time < m_time, ray2_time, ray2_time - m_time);
        Ray3f ray2 = Ray3f(ray.o, ray.d, ray2_time);


        Float path_length2 = 0;
        Spectrum throughput2           = 1.f;
        Interaction3f prev_si2         = dr::zeros<Interaction3f>();
        Float         prev_bsdf_pdf2   = 1.f;
        Bool          prev_bsdf_delta2 = true;
        BSDFContext   bsdf_ctx2;

        SurfaceInteraction3f si =
                scene->ray_intersect(ray,
                                     /* ray_flags = */ +RayFlags::All,
                                     /* coherent = */ dr::eq(depth, 0u));
        SurfaceInteraction3f si2 =
                scene->ray_intersect(ray2,
                                     /* ray_flags = */ +RayFlags::All,
                                     /* coherent = */ dr::eq(depth, 0u));


        /* Set up a Dr.Jit loop. This optimizes away to a normal loop in scalar
           mode, and it generates either a a megakernel (default) or
           wavefront-style renderer in JIT variants. This can be controlled by
           passing the '-W' command line flag to the mitsuba binary or
           enabling/disabling the JitFlag.LoopRecord bit in Dr.Jit.

           The first argument identifies the loop by name, which is helpful for
           debugging. The subsequent list registers all variables that encode
           the loop state variables. This is crucial: omitting a variable may
           lead to undefined behavior. */
        dr::Loop<Bool> loop("Time Gated Path Tracer", ppart, si, si2, sampler, ray, throughput, result,
                            eta, depth, valid_ray, prev_si, prev_bsdf_pdf,
                            prev_bsdf_delta, active, path_length, ray2, throughput2, path_length2, prev_si2, prev_bsdf_pdf2, prev_bsdf_delta2);

        /* Inform the loop about the maximum number of loop iterations.
           This accelerates wavefront-style rendering by avoiding costly
           synchronization points that check the 'active' flag. */
        loop.set_max_iterations(m_max_depth);
        

        while (loop(active)) {
            /* dr::Loop implicitly masks all code in the loop using the 'active'
               flag, so there is no need to pass it to every function */

            //SurfaceInteraction3f temp = si.adjust_time(ray2.time);
            //result = dr::select(si.is_valid(), dr::abs(dr::dot(temp.p - si.p, Vector3f(0.0,0.0,1.0))), 0.0);
            //break;

            // SurfaceInteraction3f si =
            //     scene->ray_intersect(ray,
            //                          /* ray_flags = */ +RayFlags::All,
            //                          /* coherent = */ dr::eq(depth, 0u));

            // SurfaceInteraction3f si2;

            // if(dr::any_or<true>(depth > 0)){
            //     BSDFPtr bsdf2 = si2.bsdf(ray2);

            //     SurfaceInteraction3f next_si2 = si.adjust_time(ray2.time);
                
            //     ray2 = si2.spawn_ray(dr::normalize(next_si2.p - si2.p));
            //     Float t = dr::dot(next_si2.p - si2.p, next_si2.p - si2.p);

            //     Float G = dr::select(si.is_valid(), dr::abs(dr::dot(ray.d, si.n)) / (si.t * si.t), 1.0);
            //     Float G2 = dr::select(si.is_valid(), dr::abs(dr::dot(ray2.d, si2.n)) / (si2.t * si2.t), 1.0);

            //     auto [bsdf_val_2, bsdf_pdf_2] = bsdf2->eval_pdf(bsdf_ctx2, si2, ray2.d, active);
                
            //     Float p = prev_bsdf_pdf * G / G2;
            //     throughput2 *= bsdf_val_2 / p;
            //     si2 = next_si2;
            // } else {
            //     si2 = scene->ray_intersect(ray2,
            //                             /* ray_flags = */ +RayFlags::All,
            //                             /* coherent = */ dr::eq(depth, 0u));
            // }


            // SurfaceInteraction3f 

            //valid_ray |= active && si2.is_valid();
            //result = si2.t - si.t;
            //break;

            //SurfaceInteraction3f si2 = si.adjust_time(ray2.time);
            //si2.t = (prev_si2.p - si2.p).length();

            path_length += dr::select(si.is_valid(), si.t, 0);
            path_length2 += dr::select(si2.is_valid(), si2.t, 0);

            // valid_ray |= active && si2.is_valid();
            // result = path_length - path_length2;
            // break;

            // ---------------------- Direct emission ----------------------

            /* dr::any_or() checks for active entries in the provided boolean
               array. JIT/Megakernel modes can't do this test efficiently as
               each Monte Carlo sample runs independently. In this case,
               dr::any_or<..>() returns the template argument (true) which means
               that the 'if' statement is always conservatively taken. */
            // if (dr::any_or<true>(dr::neq(si.emitter(scene), nullptr))) {
            //     DirectionSample3f ds(scene, si, prev_si);
            //     Float em_pdf = 0.f;

            //     if (dr::any_or<true>(!prev_bsdf_delta))
            //         em_pdf = scene->pdf_emitter_direction(prev_si, ds,
            //                                               !prev_bsdf_delta);

            //     // Compute MIS weight for emitter sample from previous bounce
            //     Float mis_bsdf = mis_weight(prev_bsdf_pdf, em_pdf);
                
            //     Float length_weight = evalModulationWeight(ray.time, path_length);

            //     // Accumulate, being careful with polarization (see spec_fma)
            //     result = spec_fma(
            //         throughput,
            //         ds.emitter->eval(si, prev_bsdf_pdf > 0.f) * mis_bsdf * length_weight,
            //         result);
                
            //     DirectionSample3f ds2(scene, si2, prev_si2);
            //     Float em_pdf2 = 0.f;

            //     if (dr::any_or<true>(!prev_bsdf_delta2))
            //         em_pdf2 = scene->pdf_emitter_direction(prev_si2, ds2,
            //                                               !prev_bsdf_delta2);

            //     // Compute MIS weight for emitter sample from previous bounce
            //     Float mis_bsdf2 = mis_weight(prev_bsdf_pdf2, em_pdf2);
                
            //     Float length_weight2 = evalModulationWeight(ray2.time, path_length2);

            //     // Accumulate, being careful with polarization (see spec_fma)
            //     result = spec_fma(
            //         throughput2,
            //         ds2.emitter->eval(si2, prev_bsdf_pdf2 > 0.f) * mis_bsdf2 * length_weight2,
            //         result);
            // }
            

            // Continue tracing the path at this point?
            Bool active_next = (depth + 1 < m_max_depth) && si.is_valid();

            if (dr::none_or<false>(active_next))
                break; // early exit for scalar mode

            // ---------------------- Emitter sampling ----------------------

            // Perform emitter sampling?
            BSDFPtr bsdf = si.bsdf(ray);
            Mask active_em = active_next && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            if (dr::any_or<true>(active_em)) {
                // Sample the emitter
                auto [ds, em_weight] = scene->sample_emitter_direction(
                    si, sampler->next_2d(), true, active_em);
                active_em &= dr::neq(ds.pdf, 0.f);

                /* Given the detached emitter sample, recompute its contribution
                   with AD to enable light source optimization. */
                if (dr::grad_enabled(si.p)) {
                    ds.d = dr::normalize(ds.p - si.p);
                    Spectrum em_val = scene->eval_emitter_direction(si, ds, active_em);
                    em_weight = dr::select(dr::neq(ds.pdf, 0), em_val / ds.pdf, 0);
                }

                // Evaluate BSDF * cos(theta)
                Vector3f wo = si.to_local(ds.d);
                auto [bsdf_val, bsdf_pdf] =
                    bsdf->eval_pdf(bsdf_ctx, si, wo, active_em);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // Compute the MIS weight
                Float mis_em =
                    dr::select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));

                // compute path length
                Float em_path_length = path_length + ds.dist;
                
                Float length_weight = evalModulationWeight(ray.time, em_path_length);

                // Accumulate, being careful with polarization (see spec_fma)
                result[active_em] = spec_fma(
                    throughput, 0.5 * bsdf_val * em_weight * mis_em * length_weight, result);
            }

            // Continue tracing the path at this point?
            Bool active_next2 = (depth + 1 < m_max_depth) && si2.is_valid();

            // Perform emitter sampling?
            BSDFPtr bsdf2 = si2.bsdf(ray2);
            Mask active_em2 = active_next2 && has_flag(bsdf2->flags(), BSDFFlags::Smooth);

            if (dr::any_or<true>(active_em2)) {
                // Sample the emitter
                auto [ds2, em_weight2] = scene->sample_emitter_direction(
                    si2, sampler->next_2d(), true, active_em2);
                active_em2 &= dr::neq(ds2.pdf, 0.f);

                /* Given the detached emitter sample, recompute its contribution
                   with AD to enable light source optimization. */
                if (dr::grad_enabled(si2.p)) {
                    ds2.d = dr::normalize(ds2.p - si2.p);
                    Spectrum em_val2 = scene->eval_emitter_direction(si2, ds2, active_em2);
                    em_weight2 = dr::select(dr::neq(ds2.pdf, 0), em_val2 / ds2.pdf, 0);
                }

                // Evaluate BSDF * cos(theta)
                Vector3f wo2 = si2.to_local(ds2.d);
                auto [bsdf_val2, bsdf_pdf2] =
                    bsdf2->eval_pdf(bsdf_ctx2, si2, wo2, active_em2);
                bsdf_val2 = si2.to_world_mueller(bsdf_val2, -wo2, si2.wi);

                // Compute the MIS weight
                Float mis_em2 =
                    dr::select(ds2.delta, 1.f, mis_weight(ds2.pdf, bsdf_pdf2));

                // compute path length
                Float em_path_length2 = path_length2 + ds2.dist;
                
                Float length_weight2 = evalModulationWeight(ray2.time, em_path_length2);

                // Accumulate, being careful with polarization (see spec_fma)
                result[active_em2] = spec_fma(
                    throughput2, 0.5 * bsdf_val2 * em_weight2 * mis_em2 * length_weight2, result);
            }

            // ---------------------- BSDF sampling ----------------------

            Float sample_1 = sampler->next_1d();
            Point2f sample_2 = sampler->next_2d();

            auto [bsdf_sample, bsdf_weight] =
                bsdf->sample(bsdf_ctx, si, sample_1, sample_2, active_next);
            bsdf_weight = si.to_world_mueller(bsdf_weight, -bsdf_sample.wo, si.wi);

            ray = si.spawn_ray(si.to_world(bsdf_sample.wo));

            /* When the path tracer is differentiated, we must be careful that
               the generated Monte Carlo samples are detached (i.e. don't track
               derivatives) to avoid bias resulting from the combination of moving
               samples and discontinuous visibility. We need to re-evaluate the
               BSDF differentiably with the detached sample in that case. */
            if (dr::grad_enabled(ray)) {
                ray = dr::detach<true>(ray);

                // Recompute 'wo' to propagate derivatives to cosine term
                Vector3f wo = si.to_local(ray.d);
                auto [bsdf_val, bsdf_pdf] = bsdf->eval_pdf(bsdf_ctx, si, wo, active);
                bsdf_weight[bsdf_pdf > 0.f] = bsdf_val / dr::detach(bsdf_pdf);
            }

            // ------ Update loop variables based on current interaction ------

            
            throughput *= bsdf_weight;
            eta *= bsdf_sample.eta;
            valid_ray |= active && si.is_valid() &&
                         !has_flag(bsdf_sample.sampled_type, BSDFFlags::Null);


            // Information about the current vertex needed by the next iteration
            prev_si = si;
            prev_bsdf_pdf = bsdf_sample.pdf;
            prev_bsdf_delta = has_flag(bsdf_sample.sampled_type, BSDFFlags::Delta);

            // -------------------- Stopping criterion ---------------------

            dr::masked(depth, si.is_valid()) += 1;

            Float throughput_max = dr::max(unpolarized_spectrum(throughput));

            Float rr_prob = dr::minimum(throughput_max * dr::sqr(eta), .95f);
            Mask rr_active = depth >= m_rr_depth,
                 rr_continue = sampler->next_1d() < rr_prob;

            /* Differentiable variants of the renderer require the the russian
               roulette sampling weight to be detached to avoid bias. This is a
               no-op in non-differentiable variants. */
            throughput[rr_active] *= dr::rcp(dr::detach(rr_prob));

            active = active_next && (!rr_active || rr_continue) &&
                     dr::neq(throughput_max, 0.f);

            si = scene->ray_intersect(ray,
                                     /* ray_flags = */ +RayFlags::All,
                                     /* coherent = */ dr::eq(depth, 0u));

            SurfaceInteraction3f next_si2 = si.adjust_time(ray2.time);//dr::select(si.is_valid(), si.adjust_time(ray2.time), si2);
            ray2 = dr::select(si.is_valid(), si2.spawn_ray_to(next_si2.p), ray2);
            //si2 = next_si2;
            //throughput2 *= bsdf_weight;

            
            
            auto [bsdf_val_2, bsdf_pdf_2] = bsdf2->eval_pdf(bsdf_ctx2, si2, si2.to_local(ray2.d), active_next && si.is_valid());
            // //result = dr::select(dr::neq(bsdf_pdf_2, 0), bsdf_val_2 / bsdf_pdf_2, 0);
            // //break;

            //bsdf_pdf_2 = dr::select(si.is_valid(), bsdf_pdf_2, 0.0);
            //bsdf_val_2 = dr::select(si.is_valid(), bsdf_val_2, 0.0);

            // // bsdf_val_2 = si2.to_world_mueller(bsdf_val_2, -si2.to_local(ray2.d), si2.wi);

            Mask occluded = scene->ray_test(ray2, active_next && si.is_valid());
            throughput2 = dr::select(occluded, 0.0, throughput2);
            
            next_si2.t = dr::norm(next_si2.p - si2.p);

            Float G = dr::select(si.is_valid(), dr::abs(dr::dot(ray.d, si.n)) / (si.t * si.t), 1.0);
            Float G2 = dr::select(next_si2.is_valid(), dr::abs(dr::dot(ray2.d, next_si2.n)) / (next_si2.t * next_si2.t), 1.0);
            Float p = prev_bsdf_pdf * G / G2;
            throughput2 *= dr::select(dr::neq(p, 0), bsdf_val_2 / p, 0);
            si2 = next_si2;
            
            //result = bsdf_pdf_2;
            //break;
            //result = dr::norm(next_si2.p - si.p);
            //break;
        }

        return {
            /* spec  = */ result, //dr::select(valid_ray, result, 0.f),
            /* valid = */ valid_ray
        };
    }

    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("PathIntegrator[\n"
            "  max_depth = %u,\n"
            "  rr_depth = %u\n"
            "]", m_max_depth, m_rr_depth);
    }

    /// Compute a multiple importance sampling weight using the power heuristic
    Float mis_weight(Float pdf_a, Float pdf_b) const {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        Float w = pdf_a / (pdf_a + pdf_b);
        return dr::detach<true>(dr::select(dr::isfinite(w), w, 0.f));
    }

    /**
     * \brief Perform a Mueller matrix multiplication in polarized modes, and a
     * fused multiply-add otherwise.
     */
    Spectrum spec_fma(const Spectrum &a, const Spectrum &b,
                      const Spectrum &c) const {
        if constexpr (is_polarized_v<Spectrum>)
            return a * b + c;
        else
            return dr::fmadd(a, b, c);
    }

    MI_DECLARE_CLASS()
private:
    //ScalarFloat m_tau_ps;
    //PathLenghImportanceFunctionType m_path_length_importance_function_type;
    //ScalarFloat m_delta_tau_ps;
    ScalarFloat m_illumination_modulation_frequency_mhz;
    ScalarFloat m_illumination_modulation_scale;
    ScalarFloat m_illumination_modulation_offset;
    ScalarFloat m_sensor_modulation_frequency_mhz;
    ScalarFloat m_sensor_modulation_scale;
    ScalarFloat m_sensor_modulation_offset;
    ScalarFloat m_sensor_modulation_phase_offset;
    ScalarFloat m_time;
    ScalarFloat m_hetero_frequency;
    bool m_low_frequency_component_only;
};

MI_IMPLEMENT_CLASS_VARIANT(DopplerPathIntegratorPotisional, MonteCarloIntegrator)
MI_EXPORT_PLUGIN(DopplerPathIntegratorPotisional, "Doppler Path Tracer integrator");
NAMESPACE_END(mitsuba)
