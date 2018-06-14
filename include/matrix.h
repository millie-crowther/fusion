#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <functional>
#include "function.h"

class matrix_t {
private:
    // private constructor to prevent construction
    matrix_t(float * ms);    

    // data
    float m[9];

public:
    // factory methods
    static matrix_t hessian(std::function<float(point_t)> f, point_t v, float l);
    static matrix_t jacobian(std::function<point_t(point_t)> f, point_t v, float l);

    // overridden operator
    point_t operator*(point_t v);

    // public methods
    matrix_t transpose();    
    std::vector<float> stack();
    float get(int x, int y);
    void print();
};

#endif
