import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import skimage
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
from skimage.transform import rescale, resize, downscale_local_mean
from scipy.interpolate import make_interp_spline, BSpline
import configargparse
from utils.image_utils import *

def show_image(
    expnames,
    **kwargs
):
    base_dir = os.path.join(kwargs.get('base_dir'), kwargs.get('scene_name'))
    exposure_time = kwargs.get('exposure_time', 0.0015)
    # heterodyne_frequencies = np.linspace(0.0, 1.0, 10 + 1)
    heterodyne_frequencies = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    n_columns = len(heterodyne_frequencies)
    n_rows = len(expnames)
    output_base_dir = os.path.join(kwargs.get('output_base_dir'), kwargs.get('scene_name'))

    fig, axis = plt.subplots(n_rows, n_columns, figsize=(2 * n_columns, 2 * n_rows))

    # load reference image
    reference_base_dir = os.path.join(kwargs.get('reference_base_dir', kwargs.get('base_dir')), kwargs.get('scene_name'))
    
    for i, freq in enumerate(heterodyne_frequencies):
        offset = 0.0
        ref_output_path = os.path.join(reference_base_dir, "freq_%.3f_offset_%.3f" % (freq, offset))
        reference_image = np.load(os.path.join(ref_output_path, "reference.npy")) * exposure_time
        reference_image = rgb2luminance(reference_image)
        
        images = []
        for j, expname in enumerate(expnames):
            output_path = os.path.join(base_dir, "freq_%.3f_offset_%.3f" % (freq, offset))
            image = np.load(os.path.join(output_path, "%s.npy" % expname)) * exposure_time     

            image = rgb2luminance(image)
            images.append(image)

        images= np.asarray(images)

        vmin = np.percentile(images, 10)
        vmax = np.percentile(images, 90)

        for j, expname in enumerate(expnames):
            image = images[j]
            rmse_error = np.sqrt(np.mean((image - reference_image)**2))
            relative_rmse = rmse_error / np.sqrt(np.mean((reference_image)**2))

            ax = axis[j,i]
            ax.imshow(image, vmin=vmin, vmax=vmax)
            # ax.axis('off')
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title("%.4f" % relative_rmse)
    
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.05, hspace=0.25)

    plt.savefig(os.path.join(output_base_dir, "image_over_w_r.svg"), dpi=200)
    plt.savefig(os.path.join(output_base_dir, "image_over_w_r.png"), dpi=200)
    plt.show()

if __name__=="__main__":
    parser = configargparse.ArgumentParser()
    parser.add_argument('--config', is_config_file=True, help='config file path')
    parser.add_argument("--expnumber", type=int, default=1, help="expnumber")
    parser.add_argument("--basedir", type=str, default="../", help="base directory")

    args = parser.parse_args()

    sampling_methods = [("uniform", 0), ("stratified", 16), ("antithetic", 16), ("antithetic_mirror", 16)]
    
    expnames = []
    base_dir = os.path.join(args.basedir, "results/time_spatial_sampling_comparison")
    output_base_dir = os.path.join(args.basedir, "results/images_over_hetero_frequency")
    reference_base_dir = os.path.join(args.basedir, "results/gt_images")
    for i, ts in enumerate(sampling_methods):
        t, s = ts
        expname = "%s_path_corr_depth_%d" % (t, s)
        expnames.append(expname)
    
    scene_names = ["cornell-box", "bedroom", "kitchen", "living-room-2", "soccer-ball", "veach-ajar"]
    scene_name = "cornell-box/sinusoidal"
    # scene_name = "cornell-box/sinusoidal"

    show_image(expnames, scene_name=scene_name, base_dir=base_dir,output_base_dir=output_base_dir, reference_base_dir=reference_base_dir)