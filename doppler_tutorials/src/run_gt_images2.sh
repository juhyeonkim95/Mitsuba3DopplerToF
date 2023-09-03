wave_function_type=sinusoidal
basedir=/media/juhyeon/Data1/Mitsuba3Python_final/
#python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type sinusoidal --basedir $basedir
#python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type rectangular --basedir $basedir
#python main_compare_time_spatial_sampling_methods.py --scene_name cornell-box --wave_function_type triangular --basedir $basedir
max=5
for i in `seq 2 $max`
do
    #python main_compare_time_spatial_sampling_methods.py --scene_name kitchen --wave_function_type sinusoidal --basedir $basedir
    #python main_compare_time_spatial_sampling_methods.py --scene_name soccer-ball --wave_function_type sinusoidal --basedir $basedir
    python main_compare_time_spatial_sampling_methods.py --scene_name living-room-2 --wave_function_type sinusoidal --basedir $basedir
    python main_compare_time_spatial_sampling_methods.py --scene_name veach-ajar --wave_function_type sinusoidal --basedir $basedir
    python main_compare_time_spatial_sampling_methods.py --scene_name bedroom --wave_function_type sinusoidal --basedir $basedir
done