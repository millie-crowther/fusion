#include "matrix.h"

#include "sdf.h"
#include <iostream>
#include <glm/gtx/string_cast.hpp>

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
        if (i % 3 == 2) std::cout << ']' << std::endl;
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
matrix_t::hessian(std::function<float(point_t)> f, point_t p, float l){
    float ms[9];

    point_t axes[3] = {
        point_t(l, 0, 0),
        point_t(0, l, 0),
        point_t(0, 0, l),
    };

    for (int x = 0; x < 3; x++){
        for (int y = 0; y < 3; y++){
            point_t a = p + axes[x] + axes[y];
            point_t b = p + axes[x] - axes[y];
            point_t c = p - axes[x] + axes[y];
            point_t d = p - axes[x] - axes[y];

            float p = f(a);
            float q = f(b);
            float r = f(c);
            float s = f(d);

            float u = (p - q) / (2.0f * l);
            float v = (r - s) / (2.0f * l);
         
            float res =  (u - v) / (2.0f * l);
            ms[x + y * 3] = res;
        }
    }

    return matrix_t(ms);    
}

matrix_t
matrix_t::jacobian(std::function<point_t(point_t)> f, point_t p, float l){
    point_t j[3] = {
        (f(p + point_t(l, 0, 0)) - f(p - point_t(l, 0, 0))) / (2.0f * l),
        (f(p + point_t(0, l, 0)) - f(p - point_t(0, l, 0))) / (2.0f * l),
        (f(p + point_t(0, 0, l)) - f(p - point_t(0, 0, l))) / (2.0f * l)
    };

    float ms[9];
    for (int i = 0; i < 9; i++){
        ms[i] = j[i % 3][i / 3];
    }

    return matrix_t(ms);    
}

float
matrix_t::get(int x, int y){
    return m[y * 3 + x];
}
