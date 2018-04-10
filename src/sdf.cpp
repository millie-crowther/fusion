#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

sdf_t::sdf_t(depth_map_t depths, point_t size, float l){
    this->size = size;
    this->l = l;
    this->depths = depths;

    for (int x = 0; x < size.get_x(); x += l){
        deform_field.push_back(std::vector<std::vector<point_t>>());

        for (int y = 0; y < size.get_y(); y += l){
            deform_field[x].push_back(std::vector<point_t>());

            for (int z = 0; z < size.get_z(); z += l){
                deform_field[x][y].push_back(point_t());
            }
        }
    }
}

sdf_t::~sdf_t(){
    delete depths;
}

float
sdf_t::distance(point_t p){
    // true signed distance
    float phi_true = 1;// depths->at(x).at(y) - p.get_z();
    
    // divide by delta
    float phi = phi_true / delta;
    
    // clamp to range [-1..1]
    return phi / std::max(1.0f, std::abs(phi));
}

point_t
sdf_t::distance_gradient(point_t p){

}

point_t
sdf_t::voxel_centre(point_t p){
    point_t p_div = p / l;
    point_t p_floor = point_t((int) p_div.get_x(), (int) p_div.get_y(), (int) p_div.get_z()) * l;
    return p_floor + point_t(l, l, l) / 2.0f; 
}

void 
sdf_t::fuse(canonical_sdf_t * canon, sdf_t * previous, min_params_t * ps){
    // initialise deformation field to that of the previous frame
    std::cout << "Initialising deformation field..." << std::endl;

    for (int x = 0; x < deform_field.size(); x++){
        for (int y = 0; y < deform_field[0].size(); y++){
            for (int z = 0; z < deform_field[0][0].size(); z++){
                deform_field[x][y][z] = previous->deform_field[x][y][z];
            }
        }
    }

    // rigid component
    std::cout << "Performing rigid component of reconstruction..." << std::endl;

    bool should_update;
    do {
        should_update = false;
    } while (should_update);

    // non-rigid component
    std::cout << "Performing non-rigid component of reconstruction..." << std::endl; 

    do {
        should_update = false;
    } while (should_update);

}
