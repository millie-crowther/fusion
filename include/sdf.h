#ifndef SDF_H
#define SDF_H

#include <vector>

#include "voxel.h"
#include "point.h"

class sdf_t {
public:
    // constructor
    sdf_t(int id, std::vector<std::vector<float>> * depths, float delta, point_t size, float l);

    // distance functions
    float distance(point_t p);
    point_t distance_gradient(point_t p);

private:
    // access voxel for point
    voxel_t * voxel_at(point_t p);

    // private fields:

    // unique id, ascending
    int id;

    // clamp range of distance function
    float delta;

    // size of voxel grid
    point_t size;

    // side length of each voxel
    float voxel_length;

    // depth map of frame
    std::vector<std::vector<float>> * depths;

    // voxel grid
    std::vector<std::vector<std::vector<voxel_t>>> voxels;
};

#endif
