#include "matrix.h"

#include "sdf.h"

matrix_t::matrix_t(float * ms){
    if (ms != nullptr){
        for (int i = 0; i < 9; i++){
            m[i] = ms[i];
        }
    }
}

point_t
matrix_t::operator*(point_t v){
    return point_t(
        m[0] * v.get_x() + m[1] * v.get_y() + m[2] * v.get_z(),
        m[3] * v.get_x() + m[4] * v.get_y() + m[5] * v.get_z(),
        m[6] * v.get_x() + m[7] * v.get_y() + m[8] * v.get_z()
    );
}

matrix_t
matrix_t::transpose(){
    float ms[9] = {
        m[0], m[3], m[6],
        m[1], m[4], m[7],
        m[2], m[5], m[8]
    };
    return matrix_t(ms);
}


std::vector<float>
matrix_t::stack(){
    std::vector<float> result;

    for (int x = 0; x < 3; x++){
        for (int y = 0; y < 3; y++){
            result.push_back(get(x, y));
        }
    }

    return result;
}

matrix_t
matrix_t::hessian(function_t f, point_t p){
    return matrix_t(nullptr);
}

float
matrix_t::get(int x, int y){
    return y * 3 + x;
}
