#ifndef SDF_H
#define SDF_H

#include <vector>

#include "voxel.h"
#include "point.h"

class sdf_t {
public:
    float distance(point_t p);
    point_t distance_gradient(point_t p);

private:
    voxel_t * voxel_at(point_t p);

    int id;
    point_t size;
    float voxel_length;
    
    std::vector<voxel_t> voxels;
};

#endif
