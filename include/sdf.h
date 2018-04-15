#ifndef SDF_H
#define SDF_H

#include <vector>

#include "point.h"
#include "min_params.h"

typedef std::vector<std::vector<unsigned char>> * depth_map_t;

class canon_sdf_t;

class sdf_t {
public:
    // constants 
    static constexpr float delta = 2.0f; //in millimetres

    // constructor
    sdf_t(depth_map_t depths, point_t size, float l);
    ~sdf_t();

    // fusion
    void fuse(canon_sdf_t * canon, sdf_t * previous, min_params_t * ps);

private:
    // size of voxel grid
    point_t size;

    // side length of each voxel
    float l;

    // depth map of frame
    depth_map_t depths;

    // deformation field
    std::vector<point_t> deform_field;
    
    // gradient descent
    void update_rigid(bool * cont, canon_sdf_t * canon, min_params_t * ps);
    void update_nonrigid(bool * cont, canon_sdf_t * canon, min_params_t * ps);

    point_t energy_gradient(int voxel, canon_sdf_t* c, float o_k, float o_s, float gamma, float eps);
    point_t data_energy(point_t p, point_t u, canon_sdf_t * canon);
    point_t killing_energy(point_t p, point_t u, float gamma);
    point_t level_set_energy(point_t p, point_t u, float epsilon);

    //private methods
    point_t voxel_centre(point_t p);
    point_t project(point_t p);
    float distance(point_t p);
    point_t distance_gradient(point_t p);
    point_t voxel_at(int i);
    point_t deformation_at(point_t p);
    int voxel_index(point_t p);
};

#endif
