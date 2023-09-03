import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import skimage
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
from skimage.transform import rescale, resize, downscale_local_mean
from scipy.interpolate import make_interp_spline, BSpline

target_name_dict = {
    "freq" : "$\omega$",
    "offset" : "$\phi$",
}
freq_str = "$\omega$"
offset_str = "$\phi$"


def export_error(**kwargs):
    scene_name = kwargs.get('scene_name')
    base_dir = os.path.join(kwargs.get('base_dir'), kwargs.get('scene_name'))
    reference_base_dir = os.path.join(kwargs.get('reference_base_dir', kwargs.get('base_dir')), kwargs.get('scene_name'))
    expnames = kwargs.get('expnames')
    N_heterodyne_frequencies = kwargs.get('N_heterodyne_frequencies', 10)
    N_heterodyne_offsets = kwargs.get('N_heterodyne_offsets', 10)
    exposure_time = kwargs.get('exposure_time', 0.0015)

    heterodyne_frequencies = np.linspace(0.0, 1.0, N_heterodyne_frequencies + 1)
    heterodyne_offsets = np.linspace(0.0, 1.0, N_heterodyne_offsets + 1)
    
    df = pd.DataFrame()

    for freq in heterodyne_frequencies:
        for offset in heterodyne_offsets:

            # load reference image
            ref_output_path = os.path.join(reference_base_dir, "freq_%.3f_offset_%.3f" % (freq, offset))
            reference_image = np.load(os.path.join(ref_output_path, "reference.npy")) * exposure_time
            reference_image_mean = np.mean(np.abs(reference_image))
            reference_image_max = np.max(np.abs(reference_image))
            reference_image_rmax = np.percentile(np.abs(reference_image), 95)
            
            # load image for each experiment
            output_path = os.path.join(base_dir, "freq_%.3f_offset_%.3f" % (freq, offset))

            for expname in expnames:
                image = np.load(os.path.join(output_path, "%s.npy" % expname)) * exposure_time

                # Evaluate errors
                mae_error = np.mean(np.abs(image - reference_image))                
                rmse_error = np.sqrt(np.mean((image - reference_image)**2))
                psnr = skimage.metrics.peak_signal_noise_ratio(reference_image, image, data_range=reference_image.max() - reference_image.min())
                relative_mae_error = mae_error / reference_image_mean
                relative_rmse_error = rmse_error / reference_image_mean
                snr = -10 * np.log10(relative_rmse_error)

                df2 = {
                    "freq": freq, 
                    "offset": offset, 
                    "expname":expname, 
                    "MAE": mae_error,
                    "RMSE": rmse_error,
                    "PSNR": psnr,
                    "RelativeMAE": relative_mae_error,
                    "RelativeRMSE": relative_rmse_error,
                    "SNR": snr
                    }
                df = df.append(df2, ignore_index = True)
    
    output_base_dir = os.path.join(kwargs.get('output_base_dir'), kwargs.get('scene_name'))
    df.to_csv(os.path.join(output_base_dir, "result.csv"), index=False)

def smooth(x, y, n=100):
    xnew = np.linspace(x.min(), x.max(), n)
    y_spl = make_interp_spline(x, y, k=3)
    y_smooth = y_spl(xnew)
    return xnew, y_smooth


def plot_2d_freq_vs_error_by_expname_subplot(target="freq", other_value="mean", error_type="MAE", **kwargs):
    output_base_dir = os.path.join(kwargs.get('output_base_dir'), kwargs.get('scene_name'))
    df = pd.read_csv(os.path.join(output_base_dir, "result.csv"))
    expnames = kwargs.get('expnames')
    ax = kwargs.get("ax")

    if target == "freq":
        other_target = "offset"
    else:
        other_target = "freq"

    if other_value == "mean":
        df_std = df.groupby(['expname', target])[error_type].std().reset_index(name=error_type)
        df = df.groupby(['expname', target])[error_type].mean().reset_index(name=error_type)
    else:
        df = df[df[other_target] == other_value]

    for i, expname in enumerate(expnames):
        if kwargs.get('display_names') is not None:
            label = kwargs.get('display_names')[i]
        else:
            label = expname

        df_exp = df[df['expname'] == expname]
        df_std_exp = df_std[df_std['expname'] == expname]
        if kwargs.get('alphas', None) is not None:
            alpha = kwargs.get('alphas', None).get(expname, 1.0)
        else:
            alpha = 1.0

        if kwargs.get('line_styles', None) is not None:
            if isinstance(kwargs.get('line_styles', None), dict):
                line_style = kwargs.get('line_styles', None)[expname]
            else:
                line_style = kwargs.get('line_styles', None)
            x = np.asarray(df_exp[target])
            y = np.asarray(df_exp[error_type])
            error = np.asarray(df_std_exp[error_type])
            if kwargs.get("use_smooth", True):
                xnew = np.linspace(x.min(), x.max(), 100)
                y_spl = make_interp_spline(x, y, k=3)
                y_smooth = y_spl(xnew)
                error_spl = make_interp_spline(x, error, k=3)
                error_smooth = error_spl(xnew)
                ax.plot(xnew, y_smooth, line_style, label=label, alpha=alpha)
                if(kwargs.get("plot_std", False)):
                    ax.fill_between(xnew, y_smooth-error_smooth, y_smooth+error_smooth, facecolor=line_style[0], alpha=0.2)
            else:
                ax.plot(x, y, line_style, label=label, alpha=alpha)
                if(kwargs.get("plot_std", False)):
                    ax.fill_between(x, y-error, y+error, facecolor=line_style[0], alpha=0.2)
        else:
            ax.plot(np.asarray(df_exp[target]), np.asarray(df_exp[error_type]), label=label, alpha=alpha)

    ax.locator_params(axis='x', nbins=2)
    ax.set_xlim(0.0, 1.0)
    ax.locator_params(axis='y', nbins=6)

    if "Relative" in error_type:
        ax.set_yscale('log')
    else:
        ax.ticklabel_format(style='sci', scilimits=(-3,4), axis='y')
    
    ax.set_xlabel(None)
    ax.set_ylabel(None)

def plot_3d_freq_vs_a_vs_error_by_expname(target="freq", other_value='mean', error_type='MAE', **kwargs):
    output_base_dir = os.path.join(kwargs.get('output_base_dir'), kwargs.get('scene_name'))
    df = pd.read_csv(os.path.join(output_base_dir, "result.csv"))

    if target == "freq":
        other_target = "offset"
    else:
        other_target = "freq"

    if other_value == "mean":
        df = df.groupby(['expname', target])[error_type].mean().reset_index(name=error_type)
    else:
        df = df[(df[other_target] - other_value).abs() < 1e-5]

    expnames = kwargs.get("expnames")
    expvals = kwargs.get("expvals")
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    error = np.zeros((len(expvals), 11))
    heterodyne_frequencies = np.linspace(0.0, 1.0, 10 + 1)
    a, b = np.meshgrid(heterodyne_frequencies, np.array(expvals))

    for i, val in enumerate(expvals):
        expname = expnames[i]
        df_exp = df[df['expname'] == expname]
        error[i, :] = df_exp[error_type]

    ax.plot_surface(a, b, error, cmap='coolwarm')
    
    ax.ticklabel_format(style='sci', scilimits=(-3,4), axis='z')

    if kwargs.get("hide_labels", False):
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
    else:
        ax.set_xlabel(freq_str)
        ax.set_ylabel(kwargs.get("valname"))
        ax.set_zlabel(error_type, rotation=90)

    plt.tight_layout()
    if other_value == "mean":
        output_dir = os.path.join(output_base_dir, "%s_%s" % (other_target, other_value))
    else:
        output_dir = os.path.join(output_base_dir, "%s_%.1f" % (other_target, other_value))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(os.path.join(output_dir, "plot_3d_freq_vs_a_%s_%s.png" % (error_type, kwargs.get("output_file_name", ""))), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(error)
    plt.colorbar()
    plt.savefig(os.path.join(output_dir, "plot_2d_freq_vs_a_%s_%s.png" % (error_type, kwargs.get("output_file_name", ""))), bbox_inches='tight')
    

def plot_experiment1(
    scene_names=["cornell-box"],
    wave_function_type="sinusoidal",
    base_dir=None,
    reference_base_dir=None,
    output_base_dir=None,
    **kwargs
):
    time_sampling_methods = ["uniform", "stratified", "antithetic" ,"antithetic_mirror"]
    path_correlation_depths = [16]

    color_dict = {
        "uniform": "k",
        "stratified": "r",
        "antithetic": "g",
        "antithetic_mirror": "b"
    }
    mark_dict = {
        0: "-",
        1: "--",
        2: "-.",
        3: ":",
        16:"-"
    }
    
    expnames = []
    display_names = []
    
    line_styles = {}
    alphas = {}
    
    for t in time_sampling_methods:
        for s in path_correlation_depths:
            expname = "%s_path_corr_depth_%d" % (t, s)
            expnames.append(expname)
            line_styles[expname]= "%s-" % (color_dict[t])
            alphas[expname] = 1.0
            display_name = t.replace("_", " ")
            display_names.append(display_name)

    total_scene_names = []

    for s in scene_names:
        total_scene_names.append("%s/%s" % (s, wave_function_type))

    # ["MAE", "RMSE", "PSNR", "RelativeMAE", "RelativeRMSE", "SNR"]
    target_errors = ["RMSE", "PSNR"]

    n_columns = len(scene_names)
    n_rows = len(target_errors)

    fig, axis = plt.subplots(n_rows, n_columns, figsize=(6 * n_columns, 6 * n_rows))
    
    if output_base_dir is None:
        output_base_dir = base_dir + "_plot"

    for i, scene_name in enumerate(total_scene_names):
        if not os.path.exists(os.path.join(output_base_dir, scene_name)):
            os.makedirs(os.path.join(output_base_dir, scene_name))
        
        export_error(
            base_dir=base_dir,
            reference_base_dir=reference_base_dir,
            output_base_dir=output_base_dir,
            scene_name=scene_name,
            expnames=expnames,
            **kwargs
        )

        for j, error_type in enumerate(target_errors):
            if len(axis.shape) == 2:
                ax = axis[j][i]
            elif n_columns == 1:
                ax = axis[j]
            else:
                ax = axis[i]

            plot_2d_freq_vs_error_by_expname_subplot(
                target="freq",
                base_dir=base_dir,
                output_base_dir=output_base_dir,
                scene_name=scene_name,
                expnames=expnames,
                error_type=error_type,
                line_styles=line_styles,
                alphas = alphas,
                display_names=display_names,
                plot_std=True,
                ax=ax
            )

    plt.tight_layout()
    plt.savefig(os.path.join(output_base_dir, "plot_total.svg"), dpi=600)
    plt.savefig(os.path.join(output_base_dir, "plot_total.png"), dpi=600)
    plt.show()

def plot_experiment2(
    scene_name="cornell-box",
    time_sampling_method = "antithetic",
    wave_function_type = "sinusoidal",
    base_dir=None,
    reference_base_dir=None,
    output_base_dir=None,
    **kwargs
):
    N_antithetic_shifts = kwargs.get("N_antithetic_shifts", 10)
    antithetic_shifts = np.linspace(0.0, 1.0, N_antithetic_shifts + 1)

    expnames = []
    for antithetic_shift in antithetic_shifts:
        expname = "%s_shift_%.1f" % (time_sampling_method, antithetic_shift)
        expnames.append(expname)

    scene_name = "%s/%s" % (scene_name, wave_function_type)
    if output_base_dir is None:
        output_base_dir = base_dir + "_plot"
    if not os.path.exists(os.path.join(output_base_dir, scene_name)):
        os.makedirs(os.path.join(output_base_dir, scene_name))

    export_error(
        base_dir=base_dir,
        reference_base_dir=reference_base_dir,
        output_base_dir=output_base_dir,
        scene_name=scene_name,
        expnames=expnames,
        **kwargs
    )

    for error_type in ["MAE", "RMSE", "PSNR"]:
        plot_3d_freq_vs_a_vs_error_by_expname(
            target="freq",
            other_value="mean",
            base_dir=base_dir,
            output_base_dir=output_base_dir,
            scene_name=scene_name,
            error_type=error_type,
            expnames=expnames,
            expvals=antithetic_shifts,
            valname="$a$",
            output_file_name=time_sampling_method,
            hide_labels=True
        )
        plt.close('all')
        for a in antithetic_shifts:
            plot_3d_freq_vs_a_vs_error_by_expname(
                target="freq",
                other_value=a,
                base_dir=base_dir,
                output_base_dir=output_base_dir,
                scene_name=scene_name,
                error_type=error_type,
                expnames=expnames,
                expvals=antithetic_shifts,
                valname="$a$",
                output_file_name=time_sampling_method,
                hide_labels=True
            )
            plt.close('all')


if __name__ == "__main__":
    project_dir = "/media/juhyeon/Data1/Mitsuba3Python_final/"
    
    # Experiment 1
    # reference_base_dir = os.path.join(project_dir, "results/gt_images")
    # base_dir = os.path.join(project_dir, "results/time_spatial_sampling_comparison")
    # scene_names = ["cornell-box", "living-room-2", "bedroom", "kitchen", "soccer-ball", "veach-ajar"]
    # plot_experiment1(
    #     scene_names = scene_names,
    #     reference_base_dir=reference_base_dir,
    #     base_dir = base_dir
    # )

    # Experiment 2
    # reference_base_dir = os.path.join(project_dir, "results/gt_images")
    # base_dir = os.path.join(project_dir, "results/antithetic_shift_comparison")

    # plot_experiment2(
    #     scene_name = "cornell-box",
    #     time_sampling_method = "antithetic",
    #     reference_base_dir=reference_base_dir,
    #     base_dir = base_dir
    # )

    # plot_experiment2(
    #     scene_name = "cornell-box",
    #     time_sampling_method = "antithetic_mirror",
    #     reference_base_dir=reference_base_dir,
    #     base_dir = base_dir
    # )