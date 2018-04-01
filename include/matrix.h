#ifndef MATRIX_H
#define MATRIX_H

#include "point.h"
#include <vector>

class sdf_t;

class matrix_t {
private:
    // private constructor to prevent construction
    matrix_t(float * ms);    

    // data
    float m[9];

public:
    // factory methods
    static matrix_t hessian(sdf_t * sdf, point_t v);
    static matrix_t jacobian(sdf_t * sdf, point_t v);

    // overridden operator
    point_t operator*(point_t v);

    // public methods
    matrix_t transpose();    
    std::vector<float> stack();
    float get(int x, int y);
};

#endif
