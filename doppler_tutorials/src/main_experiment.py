from program_runner import *
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
from PIL import Image
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from tqdm.contrib import itertools
import inspect
from pprint import pprint

from utils.image_utils import *
from utils.common_configs import *
import configargparse
import gc


def main():
    N = 10 + 1
    heterodyne_frequencies = np.linspace(0.0, 1.0, N)
    heterodyne_offsets = np.linspace(0.0, 1.0, N)
    antithetic_shifts = np.linspace(0.0, 1.0, N)
    n_time_samples_list = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    heterodyne_frequencies = np.linspace(0.0, 1.0, N)
    
    parser = configargparse.ArgumentParser()
    parser.add_argument('--config', is_config_file=True, help='config file path')
    parser.add_argument("--scene_name", help="scene name")
    parser.add_argument("--expnumber", type=int, default=0, help="expnumber")
    parser.add_argument("--wave_function_type", type=str, default="sinusoidal", help="waveform")
    parser.add_argument("--low_frequency_component_only", type=bool, default=True, help="low frequency only")
    parser.add_argument("--part", type=int, default=0, help="divide part for faster simulation only support by 2")
    parser.add_argument("--basedir", type=str, default="../", help="base directory")

    args = parser.parse_args()

    # split configuration into two part for speed up single scene
    if args.part == 1:
        heterodyne_offsets = np.linspace(0.0, 0.5, 6)
    elif args.part == 2:
        heterodyne_offsets = np.linspace(0.6, 1.0, 5)
    
    scene_name = args.scene_name
    wave_function_type = args.wave_function_type
    basedir = args.basedir

    # Load scene
    scene = mi.load_file(os.path.join(basedir, "scenes", scene_name, "doppler_point_correlated_sampler.xml"))

    # scene configs
    scene_configs = get_scene_configs()

    scene_config = scene_configs.get(scene_name)

    frequencies, offsets = np.meshgrid(heterodyne_frequencies, heterodyne_offsets)
    frequencies = frequencies.flatten()
    offsets = offsets.flatten()

    # render all configurations
    for f, o in itertools.product(heterodyne_frequencies, heterodyne_offsets):
        common_configs = {
            "scene_name": scene_name,
            "wave_function_type": args.wave_function_type,
            "low_frequency_component_only": args.low_frequency_component_only,
            "hetero_frequency": f, "hetero_offset": o,
            "scene": scene,
            "max_depth": scene_config.get("max_depth")
        }

        # Experiment 0. --> create reference image
        if args.expnumber == 0:
            run_scene_doppler_tof(
                expname="reference", 
                total_spp=scene_config.get("reference_spp"),
                time_sampling_method="antithetic",
                path_correlation_depth=16,
                base_dir=os.path.join(basedir, "results/gt_images"),
                exit_if_file_exists=False,
                export_png=True,
                **common_configs
            )
        
        # Experiment 1. --> different methods with different correlation depths
        elif args.expnumber == 1:
            time_sampling_methods = ["uniform", "stratified", "antithetic", "antithetic_mirror"]
            path_correlation_depths = [0, 1, 2, 16]   # 16 : full
            for time_sampling_method in time_sampling_methods:
                for path_correlation_depth in path_correlation_depths:
                    expname = "%s_path_corr_depth_%d" % (time_sampling_method, path_correlation_depth)
                    run_scene_doppler_tof(
                        expname=expname, 
                        total_spp=scene_config.get("spp"),
                        time_sampling_method=time_sampling_method,
                        path_correlation_depth=path_correlation_depth,
                        base_dir=os.path.join(basedir, "results/time_spatial_sampling_comparison_v2"),
                        exit_if_file_exists=True,
                        export_png=True,
                        **common_configs
                    )
        
        # Experiment 2. --> different methods with different correlation depths WITHOUT further stratification
        elif args.expnumber == 2:
            time_sampling_methods = ["stratified", "antithetic", "antithetic_mirror"]
            path_correlation_depths = [0, 1, 2, 16]   # 16 : full
            for time_sampling_method in time_sampling_methods:
                for path_correlation_depth in path_correlation_depths:
                    expname = "%s_path_corr_depth_%d_no_further_stratification" % (time_sampling_method, path_correlation_depth)
                    run_scene_doppler_tof(
                        expname=expname, 
                        total_spp=scene_config.get("spp"),
                        time_sampling_method=time_sampling_method,
                        path_correlation_depth=path_correlation_depth,
                        base_dir=os.path.join(basedir, "results/time_spatial_sampling_comparison"),
                        exit_if_file_exists=True,
                        use_stratified_sampling_for_each_interval=False,
                        export_png=True,
                        **common_configs
                    )

        # Experiment 3. --> different antithetic shifts
        elif args.expnumber == 3:
            time_sampling_methods = ["antithetic", "antithetic_mirror"]
            for time_sampling_method in time_sampling_methods:
                for antithetic_shift in antithetic_shifts:
                    expname = "%s_shift_%.1f" % (time_sampling_method, antithetic_shift)
                    run_scene_doppler_tof(
                        expname=expname, 
                        total_spp=scene_config.get("spp"),
                        time_sampling_method=time_sampling_method,
                        path_correlation_depth=16,
                        base_dir=os.path.join(basedir, "results/antithetic_shift_comparison"),
                        exit_if_file_exists=True,
                        antithetic_shift=antithetic_shift,
                        export_png=True,
                        **common_configs
                    )

        # Experiment 4. --> different correlation depth
        elif args.expnumber == 4:
            pass
                    


if __name__ == "__main__":
    main()
