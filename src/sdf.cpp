#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

sdf_t::sdf_t(depth_map_t depths, point_t size, float l){
    this->size = size;
    this->l = l;
    this->depths = depths;

    for (int x = 0; x < size.get_x(); x += l){
        for (int y = 0; y < size.get_y(); y += l){
            for (int z = 0; z < size.get_z(); z += l){
                deform_field.push_back(point_t());
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

    for (auto v : previous->deform_field){
        deform_field.push_back(v);
    }

    // rigid component
    std::cout << "Performing rigid component of reconstruction..." << std::endl;

    bool should_update;
    do {
        should_update = false;
        update_rigid(&should_update, ps);
    } while (should_update);

    // non-rigid component
    std::cout << "Performing non-rigid component of reconstruction..." << std::endl; 

    do {
        should_update = false;
        update_nonrigid(&should_update, ps);
    } while (should_update);

}

void
sdf_t::update_rigid(bool * cont, min_params_t * ps){
    float eta = ps->eta_rigid;
    float threshold = ps->threshold_rigid;

    for (int i = 0; i < deform_field.size(); i++){
        point_t voxel = voxel_at(i);
        point_t e = energy_rigid(voxel);
        
        if ((e * eta).length() > threshold){
            *cont = true;
        }
       
        deform_field[i] -= e * eta;
    }       
}

void
sdf_t::update_nonrigid(bool * cont, min_params_t * ps){
    float eta = ps->eta_nonrigid;
    float threshold = ps->threshold_nonrigid;

    for (int i = 0; i < deform_field.size(); i++){
        point_t voxel = voxel_at(i);
        point_t e = energy_nonrigid(voxel);
        
        if ((e * eta).length() > threshold){
            *cont = true;
        }
       
        deform_field[i] -= e * eta;
    }       

}

point_t
sdf_t::voxel_at(int i){
     point_t dim = size / l;
     int w = dim.get_x();
     int h = dim.get_y();
     int d = dim.get_z();

     int x = i / (h * d);
     int y = (i - x * h * d) / d;
     int z = i - x * h * d - y * d;

     return point_t(x, y, z) * l + point_t(l, l, l) / 2;
}

point_t
sdf_t::energy_rigid(point_t voxel){
     return point_t(); //TODO
}

point_t
sdf_t::energy_nonrigid(point_t voxel){
     return point_t(); //TODO
}
