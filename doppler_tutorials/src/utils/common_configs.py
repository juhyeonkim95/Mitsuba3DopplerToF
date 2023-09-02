def get_animation_scene_configs():
    scene_configs = {
        "falling_box" : {
            "max_depth": 4,
            "reference_spp": 1024,
            "spp": 1024,
            "animation_length": 50,
            "intervals": 5,
            "w_g": 150
        }, "bowling" : {
            "max_depth": 4,
            "reference_spp": 1024,
            "spp": 1024,
            "animation_length": 100,
            "animation_start_frame": 30,
            "animation_end_frame": 80,
            "intervals": 5,
            "w_g": 150
        }, "domino" : {
            "max_depth": 4,
            "reference_spp": 1024,
            "spp": 1024,
            "animation_length": 150,
            "intervals": 2,
            "w_g": 150
        }, "staircase2" : {
            "max_depth": 4,
            "reference_spp": 1024,
            "spp": 1024,
            "animation_length": 100,
            "intervals": 2,
            "w_g": 150
        }, "merrygoround" : {
            "max_depth": 4,
            "reference_spp": 1024,
            "spp": 1024,
            "animation_length": 80,
            "intervals": 4,
            "w_g": 150
        }
    }

    return scene_configs

def get_animation_scene_configs_shorter():
    scene_configs = {
        "falling_box" : {
            "max_depth": 4,
            "total_spp": 1024 * 4,
            "animation_length": 50,
            "intervals": 1,
            "w_g": 150
        }, "bowling" : {
            "max_depth": 4,
            "total_spp": 1024 * 64,
            "animation_length": 100,
            "animation_start_frame": 30,
            "animation_end_frame": 80,
            "intervals": 1,
            "w_g": 150
        }, "domino" : {
            "max_depth": 4,
            "total_spp": 1024 * 4,
            "animation_length": 150,
            # "animation_start_frame": 40,
            # "animation_end_frame": 90,
            "intervals": 1,
            "w_g": 150,
            #"exposure_time": 0.0015,
            #"phase_offset": 0.0,
            #"scale": 1,
            #"velocity_range": 0.5
        }, "staircase2" : {
            "max_depth": 4,
            "total_spp": 1024 * 16,
            "animation_length": 100,
            "intervals": 1,
            "w_g": 150
        }, "merrygoround" : {
            "max_depth": 4,
            "total_spp": 1024 * 16,
            "animation_length": 80,
            "intervals": 1,
            "w_g": 150
        }
    }

    return scene_configs