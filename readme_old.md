Doppler Time-of-Flight Renderer
===================================
## About
This is Mitsuba3 implementation of "Doppler Time-of-Rendering" submitted to SIGGRAPH Asia 2023.
Please also check Mitsuba0.5 implementation at []here.

New integrator named `dopplerpath` is added.
Because we implemented path correlation with custom repeated sampler, use `correlated` sampler in all cases.


To compile, follow the original Mitsuba3's compliation guide.

Followings are some input parameters.

* `time` : exposure time in sec(default : 0.0015)
* `w_g` : illumination frequency in MHz (default : 30)
* `w_f` : sensor frequency in MHz (default : 30)
* `hetero_frequency` : relative heterodyne frequency (1.0 : perfect heterodyne, 0.0 : perfec homodyne). only one of `w_f` or `hetero_frequency` should set.
* `sensor_phase_offset` : sensor phase offset in radian (default : 0)
* `sensor_modulation_function_type` : sensor modulation waveform (default : sinusoidal)
* `illumination_modulation_function_type` : illumination modulation waveform (default : sinusoidal)
* `antithetic_shift` : antithetic shift (default : 0.5)
* `time_sampling_method` : time sampling method. one of uniform, stratified, antithetic, antithetic_mirror
* `path_correlation_depth` : number of correlated path depth (0 : None, 1 : pixel) Currently only sampler correlation is supported.


We also included example configurations with result image in `config_example` folder.

The code is still not refactored yet. Also some of notations are different from Mitsuba0.5 version.
We will work on this later.
