def get_animation_scene_configs():
    scene_configs = {
        "falling_box" : {
            "max_depth": 4,
            "total_spp": 1024 * 4,
            "animation_length": 50,
            "intervals": 1,
            "w_g": 150
        }, "domino" : {
            "max_depth": 4,
            "total_spp": 1024 * 4,
            "animation_length": 150,
            "intervals": 1,
            "w_g": 150,
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

def get_scene_configs():
    scene_configs = {
        "cornell-box" : {
            "max_depth": 4,
            "reference_spp": 4096 * 32,
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
    return scene_configs