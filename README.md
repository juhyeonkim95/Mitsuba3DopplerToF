Mitsuba3 Doppler Time-of-Flight Renderer
===================================
<font size="10"> [Paper](https://arxiv.org/abs/2309.16163) / [Project](https://juhyeonkim95.github.io/project-pages/dopplertof) </font>
## About
![visualization](assets/teaser.gif)
<!-- This repository is the official Mitsuba3 implementation of "Doppler Time-of-Flight Rendering" (SIGGRAPH Asia 2023, journal paper).
Please also check Mitsuba0.5 implementation at [here](https://anonymous.4open.science/r/MitsubaDopplerToF-EAC7/README.md). -->
This repository is the official Mitsuba3 implementation of "Doppler Time-of-Flight Rendering" by Juhyeon Kim, Wojciech Jarosz, Ioannis Gkioulekas, Adithya Pediredla (SIGGRAPH Asia 2023, journal paper).
Please also check Mitsuba0.5 implementation at [here](https://github.com/juhyeonkim95/MitsubaDopplerToF).

## Install
To compile, follow the original Mitsuba3's compliation guide at [here](https://github.com/mitsuba-renderer/mitsuba3).
Our implementation only works on CUDA, so make sure to include CUDA as a backend in `mistuba.conf`.
We recommend to include `cuda_rgb` configuration.

## Parameter Explanation
New integrator named `dopplertofpath` is added for Doppler ToF rendering.
Followings are explanation for each parameter.

### ToF Related
* `time` : Exposure time in sec. (default : 0.0015) 
    * Because of precision problem the output image is not multiplied with `time`. You should multiply `time` to get a correct result.
* `w_g` : Illumination modulation frequency in MHz. (default : 30)
* `g_1` : Illumination modulation scale. (default : 0.5)
* `g_0` : Illumination modulation offset. (default : 0.5)
* `w_s` : Sensor frequency in MHz. (default : 30)
* `sensor_phase_offset` : Sensor phase offset in radian. (default : 0)
* We also provide some syntactic sugar parameters with normalization.
    * `hetero_frequency` : Relative heterodyne frequency. 0 for perfect homodyne and 1 for perfect heterodyne. This is a syntactic sugar for `w_s`. If this value is set, `w_s` is calculated from this value. (default : not used)
    * `hetero_offset` : Relative heterodyne offset. 0 for 0 radian and 1 for 2pi radian. This is a syntactic sugar for `sensor_phase_offset`. If this value is set, `sensor_phase_offset` is calculated from this value. (default : not used)

* `wave_function_type` : Modulation waveform. Refer following table for exact configuration. (default : sinusoidal)

| `wave_function_type` | Sensor Modulation | Light Modulation | Low Pass Filtered |
|-------------|-------------------|------------------|-------------------|
| `sinusoidal`  | sinusoidal        | sinusoidal       | sinusoidal        |
| `rectangular` | rectangular       | rectangular      | triangular        |
| `triangular`  | triangular        | triangular       | Corr(tri, tri)    |
| `trapezoidal` | trapezoidal       | delta            | trapezoidal       |

* `low_frequency_component_only` : Whether to use low pass filtering for modulation functions. (default : true)


### Sampling Related
* `time_sampling_method` : Times sampling method.
    * `uniform` : uniform sampling
    * `stratified` : stratified sampling
    * `antithetic` : shifted antithetic sampling
    * `antithetic_mirror` : mirrored antithetic sampling

* `antithetic_shift` : User defined antithetic shift. (default : 0.5 for `antithetic`, 0.0 for `antithetic_mirror`)

Unlike Mitsuba0.5, we only support sampler-level correlation (temporal random replay).
This exploits the advantage of parallelization at the most.
In other words, we do not explicitly correlate ray-by-ray fashion. (e.g `ray_position` or `ray_sampler` in Mitsuba0.5 version)
We are planning to add this implementation in future.


### Correlated Sampler
We implemented path correlation with custom repeated sampler, so use `correlated` sampler in all cases.
Followings are parameters used in `correlated` sampler.

* `time_correlate_number` : The number of correlated time random generator. (default : 2)
* `path_correlate_number` : The number of correlated path random generator. (default : `time_correlate_number`)
* `use_stratified_sampling_for_each_interval` : Whether to use full stratification over time. If set to `true`, it works differently by `time_sampling_mode`. (default : true)
    * `stratified` : correlated randomly over different stratum (Fig.8-(b) in the main paper)
    * `antithetic` : use stratification for primal sample (Fig.8-(e) in the main paper)
    * `antithetic_mirror` : use stratification for primal sample (Fig.8-(d) in the main paper)

Note that correlated sampler generate repeated random numbers, and `time_sampling_method` and `antithetic_shift` actually decide how to transform this into specific time samples.

## Implemention of Motion Blur
To simulate Doppler ToF rendering, we implemented motion blur which is not implemented in original Mitsuba3.
We exploit OptiX motion blur functionalities for implementation.
The default Mitsuba3 interpolation on transformation matrix seems to be inaccurate, we used a simple linear interpolation (check `transform.h` line 466).
For scene formatting, we tried to be similar with Mitsuba0.5 version.

## Usage
```
cd configs_example
mitsuba scene.xml -m cuda_rgb
```
We also included exhaustive examples in `doppler_tutorials` folder.
You can simulate various experiments in the main paper.

## Citation
If you find this useful for your research, please consider to cite:
```
(TBA)
```
