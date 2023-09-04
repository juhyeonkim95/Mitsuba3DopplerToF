from program_runner import *

import numpy as np
from PIL import Image
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm, trange
import inspect
from pprint import pprint

from utils.image_utils import *
import configargparse
import gc
from utils.common_configs import *


def main():
    parser = configargparse.ArgumentParser()
    parser.add_argument('--config', is_config_file=True, help='config file path')
    parser.add_argument("--scene_name", help="your name")
    parser.add_argument("--part", type=int, default=0, help="part")
    parser.add_argument("--wave_function_type", type=str, default="sinusoidal", help="waveform")
    parser.add_argument("--basedir", type=str, default="./", help="base directory")

    args = parser.parse_args()
    scene_name = args.scene_name
    wave_function_type = args.wave_function_type

    # (time_sampling, correlation_depth)
    methods = [
        ("uniform", 0),
        ("stratified", 16),
        ("antithetic", 16)
    ]

    # get scene_config
    scene_configs = get_animation_scene_configs()
    if not scene_name in scene_configs:
        return
    scene_config = scene_configs[scene_name]


    end = scene_config.get("animation_end_frame", scene_config["animation_length"]) - 1
    start = scene_config.get("animation_start_frame", 0)
    animation_length = (end - start) * scene_config["intervals"]
    
    start = 0
    end = animation_length

    # split configuration into two part for speed up single scene
    if args.part == 1:
        end = animation_length // 2 + 1
    elif args.part == 2:
        start = animation_length // 2 + 1

    scene_base_dir = os.path.join(args.basedir, "scenes_animation")
    hetero_offsets = (0.0, 0.25)  # phase offset of 0 and 90 degree
    
    # render all frames
    for frame_number in trange(start, end, 1):
        scene = mi.load_file(os.path.join(scene_base_dir, scene_name, "animation_%d.xml" % frame_number))
        scene_no_animation = mi.load_file(os.path.join(scene_base_dir, scene_name, "no_animation_%d.xml" % frame_number))
        
        output_file_name = "frame_%d" % (frame_number)
        exit_if_file_exists = True
        
        output_base_dir =  os.path.join(args.basedir, "results_animation")
        
        # Render GT radial velocity
        run_scene_velocity(
            scene, 
            scene_name,
            base_dir=output_base_dir,
            output_file_name="frame_%d" % frame_number,
            exit_if_file_exists=exit_if_file_exists,
            **scene_config
        )

        # Render scene radiance (standard rendering)
        run_scene_radiance(
            scene_no_animation, 
            scene_name,
            base_dir=output_base_dir,
            output_file_name="frame_%d" % frame_number,
            exit_if_file_exists=exit_if_file_exists
        )
        
        exposure_time = scene_config.get("exposure_time", 0.0015)
        common_configs = {
            "scene_name": scene_name,
            "wave_function_type": args.wave_function_type,
            "low_frequency_component_only": True,
            "scene": scene,
            "w_g": scene_config.get("w_g"),
            "exposure_time": exposure_time,
            "max_depth": scene_config.get("max_depth"),
            "exit_if_file_exists": exit_if_file_exists
        }

        # (1) Render homodyne (no variation is used for homodyne)
        homodyne_images = []
        for hetero_offset in hetero_offsets:
            output_file_name = "frame_%d" % (frame_number)
            
            vrange = 1e-3
            homodyne_image = run_scene_doppler_tof(
                time_sampling_method="antithetic",
                path_correlation_depth=16,
                base_dir=output_base_dir,
                expname=output_file_name,
                export_png=True,
                hetero_offset=hetero_offset,
                hetero_frequency=0.0,
                vmin=-vrange,
                vmax=vrange,
                **common_configs
            )
            homodyne_image = to_tof_image(homodyne_image, exposure_time=exposure_time)
            homodyne_images.append(homodyne_image)

        # (2) Render heterodyne with variations
        for m in methods:
            heterodyne_images = []
            for i, hetero_offset in enumerate(hetero_offsets):
                time_sampling_method, path_correlation_depth = m
                output_file_name = "frame_%d" % (frame_number)
                output_path = "%s/%s/%s_path_corr_depth_%d" % (scene_name, args.wave_function_type, time_sampling_method, path_correlation_depth) 
                total_spp = scene_config.get("total_spp")
                
                vrange = 1e-3
                vrange *= 1e-3
                heterodyne_image = run_scene_doppler_tof(
                    time_sampling_method=time_sampling_method,
                    path_correlation_depth=path_correlation_depth,
                    base_dir=output_base_dir,
                    expname=output_file_name,
                    export_png=True,
                    hetero_offset=hetero_offset,
                    hetero_frequency=1.0,
                    vmin=-vrange,
                    vmax=vrange,
                    total_spp=total_spp,
                    output_path=output_path,
                    **common_configs
                )
                
                heterodyne_image = to_tof_image(heterodyne_image, exposure_time=exposure_time)
                heterodyne_images.append(heterodyne_image)

                # calculate velocity using single phase
                velocity_map = calc_velocity_from_homo_hetero(homodyne_images[i], heterodyne_image, **scene_config)
                save_speed_image(velocity_map, os.path.join(output_base_dir, output_path, "velocity_%.3f" % hetero_offset), "frame_%d.png" % frame_number, **scene_config)
            
            # calculate velocity using multiple phases
            velocity_map = calc_velocity_from_homo_heteros(homodyne_images, heterodyne_images, **scene_config)
            save_speed_image(velocity_map, os.path.join(output_base_dir, output_path, "velocity"), "frame_%d.png" % frame_number, **scene_config)
            
                
                        
if __name__ == "__main__":
    main()
