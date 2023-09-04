import numpy as np
import matplotlib.pyplot as plt
import os
import cv2

def ToneMap(c, limit):
    luminance = 0.3*c[:,:,0] + 0.6*c[:,:,1] + 0.1*c[:,:,2]
    luminance = np.dstack([luminance]*3)
    col = c * 1.0 / (1.0 + luminance / limit)
    return col


def LinearToSrgb(c):
    kInvGamma = 1.0 / 2.2
    return np.power(c, kInvGamma)

def to_ldr_image(img):
    return LinearToSrgb(ToneMap(img, 1.5))

def rgb2luminance(img):
    return (0.2126 * img[:,:,0]) + (0.7152 * img[:,:,1]) + (0.0722 * img[:,:,2])

def load_tof_image(path):
    image = np.load(path)
    return to_tof_image(image)

def to_tof_image(img, exposure_time = 0.0015):
    img = np.asarray(img)
    img = rgb2luminance(img)
    img = img * exposure_time
    return img

def to_tof_image_0_5(img):
    img = np.asarray(img) - 1
    img = rgb2luminance(img)
    return img

def visualize_tof_image(img, vmin=None, vmax=None, vmin_percentile=10, vmax_percentaile=90):
    plt.figure()
    img = to_tof_image(img)
    if vmin is None:
        vmin = np.percentile(img, vmin_percentile)
    if vmax is None:
        vmax = np.percentile(img, vmax_percentaile)
    
    plt.imshow(img, vmax = vmax, vmin = vmin)

def visualize_normal_image(img):
    plt.figure()
    img = np.asarray(img)
    img = np.power(img, 1/2.2)
    plt.imshow(img)


def get_predefined_vmin_vmax(scene_name, scene_setting):
    configs = {}
    configs["cornell-box doppler_point"] = (-1e-6, 1e-6)
    configs["cornell-box doppler_area"] = (-1e-7, 1e-7)
    
    return configs["%s %s" % (scene_name, scene_setting)]

def save_radiance_image(image, output_path, filename, resize=1):
    image = np.asarray(image)
    image = to_ldr_image(image)
    fig = plt.figure()
    plt.axis('off')
    plt.imshow(image)
    fig.savefig(os.path.join(output_path, filename), bbox_inches='tight', pad_inches=0)
    plt.close('all')


def save_hdr_image(image, output_path, filename, colorbar_also=False, resize=1, 
    **kwargs):
    image = np.asarray(image)
    image = to_ldr_image(image)
    image = image[:, :, 0:3]
    image = (image * 255.0)
    image = np.clip(image, 0, 255)
    image = image.astype(np.uint8)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    if resize != 1:
        image = cv2.resize(image, (0,0), fx=resize, fy=resize) 
    cv2.imwrite(os.path.join(output_path, filename), image)

def save_speed_image(image, output_path, filename, colorbar_also=False,
     velocity_range=5, resize=1, **kwargs):
    cm = plt.get_cmap('RdBu')
    norm = plt.Normalize(-velocity_range, velocity_range)
    if len(image.shape) == 3:
        image = image[:, :, 0]

    image = cm(norm(image))
    image = image[:, :, 0:3]
    image = (image * 255.0).astype(np.uint8)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    if resize != 1:
        image = cv2.resize(image, (0,0), fx=resize, fy=resize) 
    cv2.imwrite(os.path.join(output_path, filename), image)

def save_tof_image(image, output_path, filename, 
    vmin_percentile=5, vmax_percentaile=95, 
    vmin=None, vmax=None, colorbar_also=False, resize=1):
    if colorbar_also:
        visualize_tof_image(image, vmin, vmax, vmin_percentile, vmax_percentaile)
        plt.colorbar()
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        plt.savefig(os.path.join(output_path, filename))
        plt.close('all')
        return

    cm = plt.get_cmap('viridis')
    if vmin is None:
        vmin = np.percentile(image, vmin_percentile)
    if vmax is None:
        vmax = np.percentile(image, vmax_percentaile)
    norm = plt.Normalize(vmin, vmax)

    if len(image.shape) == 3:
        image = rgb2luminance(image)
        
    image = cm(norm(image))
    image = image[:, :, 0:3]
    image = (image * 255.0).astype(np.uint8)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    if resize != 1:
        image = cv2.resize(image, (0,0), fx=resize, fy=resize) 
    
    cv2.imwrite(os.path.join(output_path, filename), image)

def calc_velocity_from_homo_hetero(homodyne, heterodyne, **kwargs):
    # homodyne_avg = np.mean(np.abs(homodyne))
    # epsilon = 0.2 * homodyne_avg
    # heterodyne_avg = np.mean(np.abs(heterodyne))
    # r = np.mean(np.abs(homodyne) < epsilon)
    
    epsilon = 5e-3 * 0.0015
    # epsilon2 = 5e-6

    # homodyne = np.where((np.abs(homodyne) < epsilon) & (homodyne > 0), epsilon, homodyne)
    # homodyne = np.where((np.abs(homodyne) < epsilon) & (homodyne < 0), -epsilon, homodyne)

    ratio = np.divide(heterodyne, homodyne, out=np.zeros_like(homodyne), where=np.abs(homodyne) > 0)
    
    # ratio = np.where((np.abs(homodyne) < epsilon) & (np.abs(heterodyne) < epsilon2), 0, ratio)

    T = kwargs.get("exposure_time", 0.0015)
    ratio = np.clip(ratio, -1, 0.999)
    delta_w = ratio * (1 / T) / (ratio - 1)
    w_g = kwargs.get("w_g", 30) * 1e6
    speed_of_light = 3e8

    velocity_map = 0.5 * delta_w * speed_of_light / w_g 
    # velocity_map = np.abs(heterodyne * homodyne) * 1e12
    # velocity_map = np.clip(velocity_map, -5, 5)

    # velocity_map = np.where(np.abs(homodyne) < epsilon, 0, 1)

    return -velocity_map

def calc_velocity_from_homo_heteros(homodynes, heterodynes, **kwargs):
    ratio_sum = 0
    ratio_confidence_sum = 0
    
    for i in range(len(homodynes)):
        homodyne = homodynes[i]
        heterodyne = heterodynes[i]
        
        ratio = np.divide(heterodyne, homodyne, out=np.zeros_like(homodyne), where=np.abs(homodyne) > 0)
        ratio_confidence = np.abs(homodyne) + 1e-5 * 0.0015

        ratio_sum += ratio * ratio_confidence
        ratio_confidence_sum += ratio_confidence
    
    ratio = ratio_sum / ratio_confidence_sum


    T = kwargs.get("exposure_time", 0.0015)
    ratio = np.clip(ratio, -1, 0.999)
    delta_w = ratio * (1 / T) / (ratio - 1)
    w_g = kwargs.get("w_g", 30) * 1e6
    speed_of_light = 3e8

    velocity_map = 0.5 * delta_w * speed_of_light / w_g 
    # velocity_map = np.abs(heterodyne * homodyne) * 1e12
    # velocity_map = np.clip(velocity_map, -5, 5)

    # velocity_map = np.where(np.abs(homodyne) < epsilon, 0, 1)

    return -velocity_map

def export_video_from_images(images, outputdir, outuput_file_name):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    image_0 = images[0]
    height = image_0.shape[0]
    width = image_0.shape[1]
    size = (width, height)
    out  = cv2.VideoWriter(os.path.join(outputdir, "%s.mp4" % outuput_file_name),  cv2.VideoWriter_fourcc(*'mp4v'), 24, size)


    for image in images:
        # image = cv2.imread(os.path.join(file_dir, file_name))
        # print(image.shape)

        out.write(image)
    
    out.release()