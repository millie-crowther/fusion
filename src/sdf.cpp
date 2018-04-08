#include "sdf.h"

#include <algorithm>
#include <cmath>

sdf_t::sdf_t(float * depths, float delta, point_t size, float l){
    this->depths = depths;
    this->delta = delta;
    this->size = size;
    this->l = l;
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


