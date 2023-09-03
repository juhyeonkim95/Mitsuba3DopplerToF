basedir=/media/juhyeon/Data1/Mitsuba3Python_final/
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type sinusoidal --basedir $basedir --expnumber 1
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type rectangular --basedir $basedir --expnumber 1
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type triangular --basedir $basedir --expnumber 1
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type trapezoidal --basedir $basedir --expnumber 1
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type sinusoidal --basedir $basedir --expnumber 2
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type rectangular --basedir $basedir --expnumber 2
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type triangular --basedir $basedir --expnumber 2
python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type trapezoidal --basedir $basedir --expnumber 2
