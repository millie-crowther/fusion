#include "matrix.h"

#include "sdf.h"
#include <iostream>

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
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z
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
matrix_t::hessian(function_t<float> f, point_t p){
    //TODO: check order of differentiation is correct
    float ms[9];
    for (int i = 0; i < 3; i++){
	for (int j = 0; j < 3; j++){
	    auto dxj = f.differentiate(j);
	    auto dxixj = dxj.differentiate(i);
            ms[i + j * 3] = dxixj(p);
	}
    }
    return matrix_t(ms);
}

matrix_t
matrix_t::jacobian(function_t<point_t> f, point_t p){
    function_t<point_t> fs[3] = {
        f.differentiate(0),
        f.differentiate(1),
        f.differentiate(2)
    };

    point_t j[3] = {
        fs[0](p),
        fs[1](p),
        fs[2](p)
    };

    float ms[9] = {
        j[0].x, j[1].x, j[2].x,
        j[0].y, j[1].y, j[2].y,
        j[0].z, j[1].z, j[2].z
    };

    return matrix_t(ms);    
}

float
matrix_t::get(int x, int y){
    return m[y * 3 + x];
}
