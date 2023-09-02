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
import configargparse
import gc


def main():
    N = 11
    heterodyne_frequencies = np.linspace(0.0, 1.0, N)
    heterodyne_offsets = np.linspace(0.0, 1.0, N)
    antithetic_shifts = np.linspace(0.0, 1.0, N)
    n_time_samples_list = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    
    parser = configargparse.ArgumentParser()
    parser.add_argument('--config', is_config_file=True, help='config file path')
    parser.add_argument("--scene_name", help="scene name")
    parser.add_argument("--expnumber", type=int, default=0, help="expnumber")
    parser.add_argument("--wave_function_type", type=str, default="sinusoidal", help="waveform")
    parser.add_argument("--low_frequency_component_only", type=bool, default=True, help="low frequency only")
    parser.add_argument("--part", type=int, default=0, help="divide part for faster simulation only support by 2")
    parser.add_argument("--basedir", type=str, default="./", help="base directory")

    args = parser.parse_args()

    # split configuration into two part
    if args.part == 1:
        heterodyne_offsets = np.linspace(0.0, 0.5, 6)
    elif args.part == 2:
        heterodyne_offsets = np.linspace(0.6, 1.0, 5)
    
    scene_name = args.scene_name
    wave_function_type = args.wave_function_type
    basedir = args.basedir

    basedir = '/media/juhyeon/Data1/Mitsuba3Python_final/'
    output_scene_name = "%s/%s" % (scene_name, wave_function_type)

    # Load scene
    scene = mi.load_file(os.path.join(basedir, "scenes", scene_name, "doppler_point_light_correlated_sampler.xml"))

    # scene configs
    scene_configs = {
        "cornell-box" : {
            "max_depth": 4,
            "reference_spp": 4096 * 4,
            "spp": 1024
        },
        "living-room-2" : {
            "max_depth": 4,
            "reference_spp": 4096 * 32,
            "spp": 1024
        },
        "veach-ajar" : {
            "max_depth": 8,
            "reference_spp": 4096 * 32,
            "spp": 1024
        },
        "soccer-ball" : {
            "max_depth": 8,
            "reference_spp": 4096 * 32,
            "spp": 1024
        },
        "bedroom" : {
            "max_depth": 8,
            "reference_spp": 4096 * 32,
            "spp": 1024
        },
        "kitchen" : {
            "max_depth": 8,
            "reference_spp": 4096 * 32,
            "spp": 1024
        }
    }
    scene_config = scene_configs.get(scene_name)

    frequencies, offsets = np.meshgrid(heterodyne_frequencies, heterodyne_offsets)
    frequencies = frequencies.flatten()
    offsets = offsets.flatten()

    for f, o in itertools.product(heterodyne_frequencies, heterodyne_offsets):
        common_configs = {
            "hetero_frequency": f,
            "hetero_offset": o,
            "scene": scene,
            "low_frequency_component_only": args.low_frequency_component_only,
            "wave_function_type": args.wave_function_type,
            "max_depth": scene_config.get("max_depth")
        }

        # GT Images
        if args.expnumber == 0:
            run_scene(
                output_scene_name, "doppler_point_correlated_sampler", "reference", 
                scene_config.get("reference_spp"),
                n_time_samples=2,
                time_sampling_method="antithetic",
                path_correlation_depth=16,
                base_dir=os.path.join(basedir, "results/gt_images"),
                exit_if_file_exists=False,
                export_png=True,
                **common_configs
            )
                        
if __name__ == "__main__":
    main()
