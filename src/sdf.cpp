#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "matrix.h"
#include "canon_sdf.h"

sdf_t::sdf_t(depth_map_t depths, point_t size, float l){
    this->size = size;
    this->l = l;
    this->depths = depths;

    for (int x = 0; x < size.get(0); x += l){
        for (int y = 0; y < size.get(1); y += l){
            for (int z = 0; z < size.get(2); z += l){
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
sdf_t::deformation_at(point_t p){
    return deform_field[voxel_index(p)];
}

int
sdf_t::voxel_index(point_t p){
    point_t p_grid = p / l;
    int x = p_grid.get(0);
    int y = p_grid.get(1);
    int z = p_grid.get(2);

    point_t grid = size / l;
    int w = grid.get(0);
    int h = grid.get(1);
    int d = grid.get(2);

    return x * h * d + y * d + z;
}

point_t
sdf_t::distance_gradient(point_t p){
    point_t u = deformation_at(p);

    auto phi = [=](point_t x){
        return this->distance(x + u);
    };
    
    function_t<float> f(l, phi);
    function_t<float> g[3] = {
        f.differentiate(0),
        f.differentiate(1),
        f.differentiate(2)
    };

    return point_t(g[0](p), g[1](p), g[2](p));
}

point_t
sdf_t::voxel_centre(point_t p){
    point_t p_div = p / l;
    point_t p_floor = point_t((int) p_div.get(0), (int) p_div.get(1), (int) p_div.get(2)) * l;
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
     int w = dim.get(0);
     int h = dim.get(1);
     int d = dim.get(2);

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
    auto phi = [=](point_t x){
        return distance(x + u);
    };
    
    auto phi_g = [=](point_t x){
        return canon->distance(x + u);
    };
    return point_t(); //TODO
}

point_t
sdf_t::level_set_energy(point_t p, point_t u, float epsilon){
    auto phi = [=](point_t x){
        return distance(x + u);
    };

    matrix_t h = matrix_t::hessian(function_t<float>(l, phi), p);

    point_t g = distance_gradient(p + u);

    float alpha = (g.length() - 1) / (g.length() + epsilon);
    return h * g * alpha;
}

point_t
sdf_t::killing_energy(point_t p, point_t u, float gamma){
    auto psi = [=](point_t p){
        return deformation_at(p);
    };

    matrix_t j = matrix_t::jacobian(function_t<point_t>(l, psi), p);
    std::vector<float> j_v = j.stack();
    std::vector<float> jt_v = j.transpose().stack();

    std::vector<float> v;
    for (int i = 0; i < j_v.size(); i++){
	v.push_back(jt_v[i] + j_v[i] * gamma);
    }

    std::vector<matrix_t> h;
    for (int i = 0; i < 3; i++){
        function_t<float> f(l, [=](point_t p){ return deformation_at(p).get(i);});
        h.push_back(matrix_t::hessian(f, p));
    }

    point_t result;
    for (int i = 0; i < 9; i++){
	result += point_t(
            h[i / 3].get(i % 3, 0) * v[i],
            h[i / 3].get(i % 3, 1) * v[i],
            h[i / 3].get(i % 3, 2) * v[i]
	);
    }

    return result * 2;
}
