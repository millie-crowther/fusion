#include "voxel.h"

#include "sdf.h"

voxel_t::voxel_t(sdf_t * sdf, point_t p){
    this->sdf = sdf;
    this->p = p;
}

bool
voxel_t::update(sdf_t * canon, min_params_t * ps){
    u -= ps->eta() * energy_gradient(canon, ps->omega_k, ps->omega_s, ps->gamma, ps->epsilon);
    return update.length() > ps->threshold;
} 

point_t 
voxel_t::energy_gradient(sdf_t * canon, float omega_k, float omega_s, float gamma, float epsilon){
    return
        data_gradient(canon, u) + 
        killing_gradient(u, gamma) * omega_k + 
        level_set_gradient(u, epsilon) * omega_s;
}

point_t
voxel_t::data_gradient(sdf_t * canon){
    return sdf->distance_gradient(p + u) * (sdf->distance(p + u) - canon->distance(p));
}

point_t
voxel_t::killing_gradient(float gamma){
    return 2;
}

point_t
voxel_t::level_set_gradient(float epsilon){
    point_t g = sdf->distance_gradient(p + u);

    float scale = (g.length() - 1) / (g.length() + epsilon);

    return scale; 
}
