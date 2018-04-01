#ifndef SDF_H
#define SDF_H

#include <vector>

#include "voxel.h"
#include "point.h"

class sdf_t {
public:
    // constructor
    sdf_t(float * depths, float delta, point_t size, float l);

    // distance functions
    float distance(point_t p);
    point_t distance_gradient(point_t p);

    // surface fusion
    void fuse(sdf_t * canon, min_params_t * ps);

private:
    // clamp range of distance function
    float delta;

    // size of voxel grid
    point_t size;

    // side length of each voxel
    float l;

    // depth map of frame
    float * depths;

    // voxel grid
    std::vector<voxel_t> voxels;
};

#endif
