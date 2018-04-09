#ifndef SDF_H
#define SDF_H

#include <vector>

#include "point.h"

typedef std::vector<std::vector<unsigned char>> * depth_map_t;

class sdf_t {
public:
    // constants 
    static constexpr float delta = 2.0f; //in millimetres

    // constructor
    sdf_t(depth_map_t depths, point_t size, float l);
    ~sdf_t();

    // distance functions
    float distance(point_t p);
    point_t distance_gradient(point_t p);

private:
    // size of voxel grid
    point_t size;

    // side length of each voxel
    float l;

    // depth map of frame
    depth_map_t depths;
};

#endif
