#ifndef MIN_PARAMS_H
#define MIN_PARAMS_H

struct min_params_t {
    // learning rates for rigid and non-rigid alignment
    float eta_rigid;
    float eta_nonrigid;

    // relative weighting of killing condition
    float omega_k;

    // relative weighting of level set condition
    float omega_s;

    // killing condition purity
    float gamma;

    // prevention of division by zero in level set gradient
    float epsilon;

    // threshold for terminating registration in mm
    float threshold_rigid;
    float threshold_nonrigid;
};

#endif
