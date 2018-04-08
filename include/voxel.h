#ifndef VOXEL_H
#define VOXEL_H

#include "min_params.h"
#include "point.h"

class sdf_t;

class voxel_t {
public:
    voxel_t(sdf_t * sdf, point_t p);

private:
    sdf_t * sdf;
    point_t p;
    point_t u;

    bool update(sdf_t * canon, min_params_t * ps);

    // energy functions used for gradient descent
    point_t energy_gradient(float omega_k, float omega_s, float gamma, float epsilon);
    point_t data_gradient(sdf_t * canon);
    point_t killing_gradient(float gamma);
    point_t level_set_gradient(float epsilon);
};

#endif
