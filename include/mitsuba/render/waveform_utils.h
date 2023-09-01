NAMESPACE_BEGIN(mitsuba)

enum EWaveformType {
    WAVE_TYPE_SINUSOIDAL = 0,
    WAVE_TYPE_RECTANGULAR = 1,
    WAVE_TYPE_TRIANGULAR = 2,
    WAVE_TYPE_TRAPEZOIDAL = 3
};

Float evalModulationFunctionValue(Float _t, uint32_t function_type) const{
    Float t = dr::fmod(_t, 2 * M_PI);
    switch(function_type){
        case WAVE_TYPE_SINUSOIDAL: return dr::cos(t);
        case WAVE_TYPE_RECTANGULAR: return dr::select(dr::abs(t-M_PI) > 0.5 * M_PI, 1, -1); //return dr::sign(dr::cos(t));
        case WAVE_TYPE_TRIANGULAR: return dr::select(t < M_PI, 1 - 2 * t / M_PI, -3 + 2 * t / M_PI);
    }
    return dr::cos(t);
}

Float evalModulationFunctionValueLowPass(Float _t, uint32_t function_type) const{
    Float t = dr::fmod(_t, 2 * M_PI);
    switch(function_type){
        case WAVE_TYPE_SINUSOIDAL: return dr::cos(t);
        case WAVE_TYPE_RECTANGULAR: {
            Float a = t / M_PI;
            Float b = 2 - a;
            Float c = dr::select(a < b, a, b);
            return 2 - 4 * c; //a < 1 ? 1 - 2 * a : 1 - 2 * b;
        }
        case WAVE_TYPE_TRIANGULAR: {    
            Float a = t / M_PI;
            Float b = 2 - a;
            Float c = dr::select(a < b, a, b);
            return (4 * c * c * c - 6 * c * c + 1) * 2.0 / 3.0;
        }
        case WAVE_TYPE_TRAPEZOIDAL: {    
            Float a = t / M_PI;
            Float b = 2 - a;
            Float c = dr::select(a < b, a, b);
            Float r = 2 - 4 * c;
            return dr::clamp(2.0 * r, -2.0, 2.0);
        }
    }
    return dr::cos(t);
}