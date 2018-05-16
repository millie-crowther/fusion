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

void
matrix_t::print(){
    for (int i = 0; i < 9; i++){
        if (i % 3 == 0) std::cout << '[';
        std::cout << m[i];
        if (i % 3 != 2) std::cout << ' ';
        if (i % 3 == 2) std::cout << ']';
    }
    std::cout << std::endl;
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
    float ms[9];
    for (int i = 0; i < 3; i++){
	for (int j = 0; j < 3; j++){
	    auto dxj = f.differentiate(j);
	    auto dxixj = dxj.differentiate(i);
            ms[i + j * 3] = dxixj(p);
	}
    }

    auto r = matrix_t(ms);    
    std::cout << "Hessian: " << std::endl;
    r.print();
    return r;
}

matrix_t
matrix_t::jacobian(function_t<point_t> f, point_t p){
    point_t j[3] = {
        f.differentiate(0)(p),
        f.differentiate(1)(p),
        f.differentiate(2)(p)
    };

    float ms[9];
    for (int i = 0; i < 9; i++){
        ms[i] = j[i % 3][i / 3];
    }

    auto r = matrix_t(ms);    
    std::cout << "Jacobian: " << std::endl;
    r.print();
    return r;
}

float
matrix_t::get(int x, int y){
    return m[y * 3 + x];
}
