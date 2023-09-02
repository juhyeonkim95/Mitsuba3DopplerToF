
import mitsuba as mi
mi.set_variant('cuda_rgb')
import os
import numpy as np
from utils.image_utils import *
from tqdm import tqdm, trange


TIME_SAMPLING_UNIFORM = 0
TIME_SAMPLING_STRATIFIED = 1
TIME_SAMPLING_ANTITHETIC = 2
TIME_SAMPLING_ANTITHETIC_MIRROR = 3
TIME_SAMPLING_PERIODIC = 4
TIME_SAMPLING_REGULAR = 5


SPATIAL_CORRELATION_NONE = 0
SPATIAL_CORRELATION_PIXEL = 1
SPATIAL_CORRELATION_SAMPLER = 2

def str_to_int(value):
    if isinstance(value, int):
        return value

    str_dict = {
        "uniform": TIME_SAMPLING_UNIFORM,
        "stratified": TIME_SAMPLING_STRATIFIED,
        "antithetic": TIME_SAMPLING_ANTITHETIC,
        "antithetic_mirror": TIME_SAMPLING_ANTITHETIC_MIRROR,
        "periodic": TIME_SAMPLING_PERIODIC,
        "regular": TIME_SAMPLING_REGULAR,

        "none": SPATIAL_CORRELATION_NONE,
        "pixel": SPATIAL_CORRELATION_PIXEL,
        "sampler": SPATIAL_CORRELATION_SAMPLER
    }
    return str_dict[value]


def render_image_multi_pass(scene, integrator, single_pass_spp, total_pass, show_progress=False):
    img_sum = None
    if show_progress:
        for i in trange(total_pass):
            img_i = integrator.render(scene, seed=i, spp=single_pass_spp)
            if i == 0:
                img_sum = img_i
            else:
                img_sum += img_i
    else:
        img_i_old = 0
        for i in range(total_pass):
            img_i = integrator.render(scene, seed=i, spp=single_pass_spp)
            if i == 0:
                img_sum = img_i
            else:
                img_sum += img_i
            img_i_old = img_i

    img = img_sum / total_pass
    return img

def run_scene_velocity(scene, scene_name, **kwargs):
    velocity_rendering_config_dict = {
        'type': 'velocity'
    }
    total_spp = kwargs.get("total_spp", 1024)

    output_path = os.path.join(kwargs.get("base_dir"), scene_name, "velocity_gt")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file_name = kwargs.get("output_file_name")
    numpy_output_file_name = os.path.join(output_path, "%s.npy" % output_file_name)

    if os.path.exists(numpy_output_file_name) and kwargs.get("exit_if_file_exists", False):
       print("File already exists!")
       return

    integrator_velocity = mi.load_dict(velocity_rendering_config_dict)
    single_pass_spp = min(1024, total_spp)
    show_progress = kwargs.get("show_progress", False)
    img_velocity = render_image_multi_pass(scene, integrator_velocity, single_pass_spp, total_spp // single_pass_spp, show_progress=show_progress)
    np.save(numpy_output_file_name, img_velocity)
    save_speed_image(img_velocity[:,:,0], output_path, "%s.png" % output_file_name,  **kwargs)
    

def run_scene_radiance(scene, scene_name, **kwargs):
    config_dict = {
        'type': 'path',
        'max_depth': kwargs.get("max_depth", 4), 
    }
    total_spp = kwargs.get("total_spp", 1024)

    output_path = os.path.join(kwargs.get("base_dir"), scene_name, "radiance")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file_name = kwargs.get("output_file_name")
    numpy_output_file_name = os.path.join(output_path, "%s.npy" % output_file_name)

    if os.path.exists(numpy_output_file_name) and kwargs.get("exit_if_file_exists", False):
       print("File already exists!")
       return

    integrator_velocity = mi.load_dict(config_dict)
    single_pass_spp = min(1024, total_spp)
    
    show_progress = kwargs.get("show_progress", False)
    img = render_image_multi_pass(scene, integrator_velocity, single_pass_spp, total_spp // single_pass_spp, show_progress=show_progress)
    np.save(numpy_output_file_name, img)
    save_hdr_image(img, output_path, "%s.png" % output_file_name)

def run_scene(scene_name, light_setting, expname, total_spp, 
        hetero_frequency=1.0, hetero_offset=0.0, 
        antithetic_shift=None, 
        n_time_samples=2, 
        time_sampling_method="antithetic",
        path_correlation_depth=16,
        use_path_correlation=True,
        xml_file=None, scene=None, integrator_config=None, 
        exit_if_file_exists=True,
        **kwargs):
    
    if antithetic_shift is None:
        if time_sampling_method == "antithetic":
            antithetic_shift = 0.5
        else:
            antithetic_shift = 0.0
            
    output_file_name = kwargs.get("output_file_name", expname)
    
    output_path = os.path.join(kwargs.get("base_dir"), scene_name, "freq_%.3f_offset_%.3f" % (hetero_frequency, hetero_offset) )
    output_file = os.path.join(output_path, "%s.npy" % output_file_name)
    # print(output_path)
    
    if os.path.exists(output_file) and exit_if_file_exists:
        print("File already exists!")
        return np.load(output_file)


    if scene is None:
        scene = mi.load_file(os.path.join("../scenes", scene_name, "%s.xml"%light_setting))
    # scene.integrator()

    # pprint(inspect.getmembers(scene))
    #mi.set_log_level(mi.LogLevel.Error)
    # params = mi.traverse(scene)
    #print(params)
    #return

    # img = mi.render(scene)
    config_dict = {
        'type': 'dopplertofpath',
        'max_depth': kwargs.get("max_depth", 4), 
        'w_g': kwargs.get("w_g", 30.0), 
        'hetero_frequency': hetero_frequency,
        'hetero_offset': hetero_offset,
        'is_doppler_integrator': True,
        'time': kwargs.get("exposure_time", 0.0015),
        'antithetic_shift': antithetic_shift,
        'time_sampling_method': time_sampling_method,
        'path_correlation_depth': path_correlation_depth,
        'low_frequency_component_only': kwargs.get("low_frequency_component_only", True),
        'wave_function_type': kwargs.get("wave_function_type", "sinusoidal")
    }
    integrator_doppler = mi.load_dict(config_dict)


    # sampler = mi.load_dict({
    #     'type': 'correlated',
    #     'antithetic_shift': antithetic_shift,
    #     'time_correlate_number': n_time_samples,
    #     'path_correlate_number': n_time_samples,
    #     'use_stratified_sampling_for_antithetic': False
    # })

    # scene.sensors()[0].sampler = sampler

    # single_pass_spp = 1024

    output_path = os.path.join(kwargs.get("base_dir"), scene_name, "freq_%.3f_offset_%.3f" % (hetero_frequency, hetero_offset) )
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    single_pass_spp = min(1024, total_spp)
    
    show_progress = kwargs.get("show_progress", False)
    img_dop = render_image_multi_pass(scene, integrator_doppler, single_pass_spp, total_spp // single_pass_spp, show_progress=show_progress)
    np.save(os.path.join(output_path, "%s.npy" % output_file_name), img_dop)

    if kwargs.get("export_png", False):
        # img_dop = np.load(os.path.join(output_path, "%s.npy" % output_file_name))
        # visualize_tof_image(img_dop, vmin=kwargs.get("vmin", None), vmax=kwargs.get("vmax", None))
        # plt.colorbar()
        # plt.axis('off')
        # plt.savefig(os.path.join(output_path, "%s.png" % output_file_name), bbox_inches='tight', pad_inches=0)
        # plt.close('all')
        # print("AAAAAAAAAAAAAa")
        save_tof_image(to_tof_image(img_dop, kwargs.get("exposure_time", 0.0015)), output_path, "%s.png" % output_file_name)
    
    return np.asarray(img_dop)


    # plt.colorbar()
    # plt.savefig(os.path.join(output_path, "%s.png" % output_file_name))
    # plt.close()
