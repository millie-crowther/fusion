#ifndef MIN_PARAMS_H
#define MIN_PARAMS_H

#include "point.h"

namespace fusion_mode {
    const int NIL = 0;
    const int CPU = 1;
    const int CPU_MULTITHREAD = 2;
    const int GPU = 3;
}

struct min_params_t {
    // learning rates for rigid and non-rigid alignment
    float eta;

    // relative weighting of killing condition
    float omega_k;

    // relative weighting of level set condition
    float omega_s;

    // killing condition purity
    float gamma;

    // prevention of division by zero in level set gradient
    float epsilon;

    // threshold for terminating registration in mm
    float threshold;

    int mode;
    int frames;

    point_t size;
  
    float voxel_length;

    float delta;
   
    float camera_fx;
    float camera_fy;
    float camera_cx;
    float camera_cy;

    float sdf_eta; //TODO: find a good value for this
};

#endif
