#include "sdf.h"

#include <cmath>

sdf_t::sdf_t(int id, std::vector<std::vector<float>> * depths, float delta, point_t size, float l){
    this->id = id;
    this->depths = depths;
    this->delta = delta;
    this->size = size;
    this->voxel_length = l;
}

float
sdf_t::distance(point_t p){
    // TODO: apply transformation associated with voxel that p 
    //       occupies before rest of function

    int x = p.get_x();
    int y = p.get_y();

    // check on indices
    if (x < 0 || y < 0 || x >= depths->size() || y >= depths->at(x).size()){
        // TODO
    }

    // true signed distance
    float phi_true = depths->at(x).at(y) - p.get_z();
    
    // divide by delta
    float phi = phi_true / delta;
    
    // clamp to range [-1..1]
    if (fabs(phi) > 1){
        return phi / fabs(phi);
    } else {
        return phi;
    }
}

point_t
sdf_t::distance_gradient(point_t p){

}


