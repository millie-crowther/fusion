#ifndef VOXEL_H
#define VOXEL_H

struct min_params_t {
    // step size of gradient descent
    float alpha;

    // relative weighting of killing condition
    float omega_k;

    // relative weighting of level set condition
    float omega_s;

    // killing condition purity
    float gamma;

    // minimum size of quotient in level set energy
    float epsilon;

    // threshold for terminating registration
    float threshold;
};

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
    point_t energy_gradient(point_t u, float omega_k, float omega_s, float gamma, float epsilon);
    point_t data_gradient(sdf_t * canon, point_t u);
    point_t killing_gradient(point_t u, float gamma);
    point_t level_set_gradient(point_t u, float epsilon);
};

#endif
