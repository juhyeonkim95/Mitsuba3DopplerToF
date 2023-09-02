import os


def option_to_folder(**kwargs):
    force_collocated_point_light = kwargs.get("force_collocated_point_light", True)
    max_depth = kwargs.get("max_depth", 2)
    # homogeneous = kwargs.get("homogeneous", True)

    heterodyne_frequency = kwargs.get("heterodyne_frequency", 1.0)
    heterodyne_offset = kwargs.get("heterodyne_offset", 0.0)

    option_point_light = "point_light" if force_collocated_point_light else "original_light"

    if kwargs.get("velocity_dictionary", None) is not None:
        velocity_dictionary_scene = kwargs.get("velocity_dictionary").get(kwargs.get("scene_name"), kwargs.get("velocity_dictionary")["default"])
        
        if "sensor" in velocity_dictionary_scene:
            option_point_light += "_moving_camera"


    option_max_bounce = "max_depth_%d" % max_depth
    # option_homogeneous = "homo" if homogeneous else "hetero"

    heterodyne_precision = kwargs.get("heterodyne_precision", 2)
    if heterodyne_precision == 2:
        folder_name = "freq_%.2f_offset_%.2f" % (heterodyne_frequency, heterodyne_offset)
    elif heterodyne_precision == 3:
        folder_name = "freq_%.3f_offset_%.3f" % (heterodyne_frequency, heterodyne_offset)

    if "low_frequency_component_only" in kwargs:
        folder_name += "_low_only" if kwargs.get("low_frequency_component_only") else "_full"

    output_folder = os.path.join(option_point_light, option_max_bounce, folder_name)
    return output_folder

def float_to_no_dot(v):
    a = int((v * 100)) // 100
    b = int((v * 100)) % 100

    return "%d_%02d" % (a, b)

def float_array_to_no_dot(vs):
    s = ""
    str_no_dots = []
    for v in vs:
        str_no_dots.append(float_to_no_dot(v))
    return "_".join(str_no_dots)

    #     s= "%s_%s" % (s, float_to_no_dot(v))
    # return s
