from program_runner import *

# mi.set_log_level(mi.LogLevel.Info)
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
from PIL import Image
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm, trange
import inspect
from pprint import pprint

from utils.image_utils import *
import argparse
import gc
from utils.common_configs import *


def main():
    heterodyne_frequencies = [1.0]
    heterodyne_offsets = [0.0]
    
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-n", "--scene_name", help="your name")
    argParser.add_argument("-p", "--part", type=int, default=0, help="part")
    # argParser.add_argument("-fp", "--fpart", type=int, default=0, help="fpart")
    argParser.add_argument("-e", "--expnumber", type=int, default=0, help="expnumber")

    args = argParser.parse_args()
    scene_name = args.scene_name

    
    wave_target_configs = {
        "sinusoidal_low_freq_only" : {
            "low_frequency_component_only": True
        }
    }

    scene_configs = get_animation_scene_configs()

    if not scene_name in scene_configs:
        return
    
    scene_config = scene_configs[scene_name]

    animation_length = (scene_config["animation_length"] - 1) * scene_config["intervals"]
    start = 0
    end = animation_length
    if args.part == 1:
        end = animation_length // 2
    elif args.part == 2:
        start = animation_length // 2

    base_dir = "../scenes_animation"

    for frame_number in trange(start, end, 1):
        scene = mi.load_file(os.path.join(base_dir, scene_name, "animation_%d.xml" % frame_number))
        scene_no_animation = mi.load_file(os.path.join(base_dir, scene_name, "no_animation_%d.xml" % frame_number))
    
        for wave_name, wave_config in wave_target_configs.items():
            output_scene_name = "%s/%s" % (scene_name, wave_name)
            exit_if_file_exists = False

            #if args.part == 2:
            run_scene_velocity(
                scene, 
                scene_name,
                base_dir="../results_animation/linear_interpolate",
                output_file_name="frame_%d" % frame_number,
                exit_if_file_exists=exit_if_file_exists,
                **kwargs
            )
            # else:
            #     pass
            run_scene_radiance(
                scene_no_animation, 
                scene_name,
                base_dir="../results_animation/linear_interpolate",
                output_file_name="frame_%d" % frame_number,
                exit_if_file_exists=exit_if_file_exists
            )
            # continue

            for f in heterodyne_frequencies:
                for o in heterodyne_offsets:
                    total_spp = scene_config.get("reference_spp")
                    # print(f, o, "freq, offset")
                    common_configs = {
                        "hetero_frequency": f,
                        "hetero_offset": o,
                        "scene": scene,
                        "w_g": scene_config.get("w_g"),
                        "wave_config": wave_config,
                        "max_depth": scene_config.get("max_depth"),
                        "exit_if_file_exists": exit_if_file_exists
                    }
                    vrange = 1e-3
                    if f > 0.5:
                        vrange *= 1e-3
                    # GT Images
                    run_scene(
                        output_scene_name, "doppler_point_correlated_sampler", "reference", 
                        total_spp,
                        n_time_samples=2,
                        time_sampling_method="antithetic",
                        path_correlation_depth=16,
                        base_dir="../results_animation/linear_interpolate",
                        output_file_name = "frame_%d" % frame_number,
                        export_png=True,
                        vmin=-vrange,
                        vmax=vrange,
                        exposure_time=scene_config.get("exposure_time", 0.0015)
                        **common_configs
                    )
            # GT Images
            # run_scene(
            #     output_scene_name, "doppler_point_correlated_sampler", "reference", 
            #     scene_config.get("reference_spp"),
            #     n_time_samples=2,
            #     time_sampling_method="antithetic",
            #     path_correlation_depth=16,
            #     base_dir="../results_animation/linear_interpolate",
            #     exit_if_file_exists=True,
            #     **common_configs
            # )
                        
if __name__ == "__main__":
    main()
