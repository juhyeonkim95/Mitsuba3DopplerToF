Doppler Time-of-Flight Renderer
===================================
## About
![visualization](assets/teaser.gif)

This repository is the official Mitsuba3 implementation of "Doppler Time-of-Flight Rendering" by Juhyeon Kim, Wojciech Jarosz, Ioannis Gkioulekas, Adithya Pediredla (SIGGRAPH Asia 2023, journal paper).
Please also check Mitsuba0.5 implementation at [here](https://github.com/juhyeonkim95/MitsubaDopplerToF).

## Install
To compile, follow the original Mitsuba3's compliation guide at [here](https://github.com/mitsuba-renderer/mitsuba3).
Our implementation only works on CUDA, so make sure to include CUDA as a backend in `mistuba.conf`.

## Parameter Explanation
New integrator named `dopplertofpath` is added for Doppler ToF rendering.
Followings are explanation for each parameter.

### ToF Related
* `time` : Exposure time in sec. (default : 0.0015)
* `w_g` : Illumination modulation frequency in MHz. (default : 30)
* `g_1` : Illumination modulation scale. (default : 0.5)
* `g_0` : Illumination modulation offset. (default : 0.5)
* `w_s` : Sensor frequency in MHz. (default : 30)
* `sensor_phase_offset` : Sensor phase offset in radian. (default : 0)

We also provide some syntactic sugar parameters.
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
* `time_sampling_mode` : Times sampling method.
    * `uniform` : uniform sampling
    * `stratified` : stratified sampling
    * `antithetic` : shifted antithetic sampling
    * `antithetic_mirror` : mirrored antithetic sampling
    * `analytic` : analytic integration (biased)

* `antithetic_shifts` : User defined antithetic shifts. Multiple input is available separated by underbar. (e.g 0.5 for single antithetic sample or 0.12_0.35 two antithetic samples) (default : 0.5 for `antithetic`, 0.0 for `antithetic_mirror`)
* `antithetic_shifts_number` : Number of antithetic shifts with equal interval. If this value is set, this is used instead of `antithetic_shifts`. This is also used for number of stratum for `stratified`. (default : 0)
* `m_use_full_time_stratification` : Whether to use full stratification over time. If set to `true`, it works differently by `time_sampling_mode`. (default : false)
    * `stratified` : correlated randomly over different stratum (Fig.8-(b) in the main paper)
    * `antithetic` : use stratification for primal sample (Fig.8-(e) in the main paper)
    * `antithetic_mirror` : use stratification for primal sample (Fig.8-(d) in the main paper)

* `spatial_correlation_mode` : Spatial correlation methods. Note that methods start with `ray` explicitly correlate two paths in ray-by-ray style.
    * `none` : no correlation
    * `pixel` : repeat camera ray (use same pixel coordinate) between multiple rays
    * `sampler` : repeat sampler between multiple rays
    * `ray_position` : correlate intersection point between two rays
    * `ray_sampler` : repeat sampler between two rays
    * `ray_selective` : select one of `ray_position`, `ray_sampler` based on material property

* `force_constant_attenuation` : Whether to force constant attenuation as [Heide, 2015] did. `true` is using zero-order Taylor approximation while `false` is using first-order Taylor approximation. Only used for `analytic`. (default : false)
* `primal_antithetic_mis_power` : MIS power for primal and antithetic sample. Only used for other than `analytic`. Refer Sec 4.1 in supplementary material for details. (default : 1.0)

### Others
* `image_offset` :  An output image can be negative, so add offset to make it positive.

## Usage
```
mitsuba -L error config_example/example.xml
```
We also included exhaustive example configurations with result image.

## Citation
If you find this useful for your research, please consider to cite:
```
(TBA)
```
