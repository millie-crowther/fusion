#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "function.h"

class matrix_t {
private:
    // private constructor to prevent construction
    matrix_t(float * ms);    

    // data
    float m[9];

public:
    // factory methods
    static matrix_t hessian(function_t<float> f, point_t v);
    static matrix_t jacobian(function_t<point_t> f, point_t v);

    // overridden operator
    point_t operator*(point_t v);

    // public methods
    matrix_t transpose();    
    std::vector<float> stack();
    float get(int x, int y);
    void print();
};

#endif
