#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

template<class T>
class function_t {
private:
    T function;

public:
    function_t(T f){
        function = f;
    }

    float eval(point_t p){
        return function(p);
    }

    function_t differentiate(int axis){
        float l = 1.0f; //TODO

        point_t axes[3] = {
            point_t(l, 0, 0),
            point_t(0, l, 0),
            point_t(0, 0, l)
        };
        point_t d = axes[axis];
    
        return function_t([function, d](point_t p){
            return (function(p + d) - function(p - d)) / (2 * l);
        });
    }
};

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
sdf_t::fuse(canon_sdf_t * canon, sdf_t * previous, min_params_t * ps){
    // initialise deformation field to that of the previous frame
    for (auto v : previous->deform_field){
        deform_field.push_back(v);
    }

    // rigid component
    bool should_update;
    int i;
    for (i = 0; i == 0 || should_update; i++){
        should_update = false;
        update_rigid(&should_update, canon, ps);
    }
    std::cout << "Rigid transformation update converged in " << i << " iterations." << std::endl;

    // non-rigid component
    for (i = 0; i == 0 || should_update; i++){
        should_update = false;
        update_nonrigid(&should_update, canon, ps);
    }
    std::cout << "Non-rigid transformation update converged in " << i << " iterations." << std::endl;
}

void
sdf_t::update_rigid(bool * cont, canon_sdf_t * canon, min_params_t * ps){
    for (int i = 0; i < deform_field.size(); i++){
        point_t e = data_energy(voxel_at(i), deform_field[i], canon);
        point_t u = e * ps->eta_rigid;
        
        if (u.length() > ps->threshold_rigid){
            *cont = true;
        }
       
        deform_field[i] -= u;
    }       
}

void
sdf_t::update_nonrigid(bool * cont, canon_sdf_t * canon, min_params_t * ps){
    for (int i = 0; i < deform_field.size(); i++){
        point_t e = energy_gradient(i, canon, ps->omega_k, ps->omega_s, ps->gamma, ps->epsilon);
        point_t u = e * ps->eta_nonrigid; 
      
        if (u.length() > ps->threshold_nonrigid){
            *cont = true;
        }
       
        deform_field[i] -= u;
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

     return voxel_centre(point_t(x, y, z) * l);
}

point_t
sdf_t::energy_gradient(int voxel, canon_sdf_t * c, float o_k, float o_s, float gamma, float eps){
     point_t p = voxel_at(voxel);
     point_t u = deform_field[voxel];
     return 
         data_energy(p, u, c) +
         killing_energy(p, u, gamma) * o_k +
         level_set_energy(p, u, eps) * o_s;
}

point_t
sdf_t::data_energy(point_t p, point_t u, canon_sdf_t * canon){
    return point_t(); //TODO
}

point_t
sdf_t::level_set_energy(point_t p, point_t u, float epsilon){
    return point_t(); //TODO
}

point_t
sdf_t::killing_energy(point_t p, point_t u, float gamma){
    return point_t(); //TODO
}
