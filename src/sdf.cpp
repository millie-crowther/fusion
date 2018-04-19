#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "matrix.h"
#include "canon_sdf.h"

const point_t sdf_t::size = point_t(80);

sdf_t::sdf_t(depth_map_t depths, bool is_multi){
    this->depths = depths;
    this->is_multi = is_multi;
    
    for (int x = 0; x < size.get(0); x += l){
        for (int y = 0; y < size.get(1); y += l){
            for (int z = 0; z < size.get(2); z += l){
                deform_field.push_back(point_t());
            }
        }
    }

    phi = new function_t<float>([=](point_t p){
        return distance(p);
    });

    psi = new function_t<point_t>([=](point_t p){
        return deformation_at(p);
    });

    psi_u = new function_t<float>([=](point_t p){
        return deformation_at(p).get(0);
    });
    
    psi_v = new function_t<float>([=](point_t p){
        return deformation_at(p).get(1);
    });
    
    psi_w = new function_t<float>([=](point_t p){
        return deformation_at(p).get(2);
    });
}

sdf_t::~sdf_t(){
    delete depths;
    delete phi;
    delete psi;
    delete psi_u;
    delete psi_v;
    delete psi_w;
}

float
sdf_t::distance(point_t p){
    point_t x = p + deformation_at(p);

    // true signed distance
    float phi_true = 1;// depths->at(x).at(y) - p.get_z();
    
    // divide by delta
    float phi1 = phi_true / delta;
    
    // clamp to range [-1..1]
    return phi1 / std::max(1.0f, std::abs(phi1));
}

point_t
sdf_t::deformation_at(point_t p){
    int i = voxel_index(p);
    if (i < 0 || i >= deform_field.size()){
	// should only happen if doing gradient check near boundaries...
	return point_t();
    } else {
        return deform_field[i];
    }
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
    auto phi = [=](point_t x){
        return this->distance(x);
    };
   
    function_t<float> f(phi);
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
    for (int i = 0; i < previous->deform_field.size(); i++){
        deform_field[i] = previous->deform_field[i];
    }

    // rigid component
    bool should_update = true;
    int i;
    for (i = 0; should_update; i++){
        should_update = false;
        update_rigid(&should_update, canon, ps);
    }
    std::cout << "Rigid transformation update converged in " << i << " iterations." << std::endl;

    // non-rigid component
    should_update = true;
    for (i = 0; should_update; i++){
        should_update = false;

        update_nonrigid(&should_update, canon, ps);
    }
    std::cout << "Non-rigid transformation update converged in " << i << " iterations." << std::endl;
}

void
sdf_t::update_rigid(bool * cont, canon_sdf_t * canon, min_params_t * ps){
    for (int i = 0; i < deform_field.size(); i++){
        point_t e = data_energy(voxel_at(i), canon);
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
         data_energy(p, c) +
         killing_energy(p, gamma) * o_k +
         level_set_energy(p, eps) * o_s;
}

point_t
sdf_t::data_energy(point_t p, canon_sdf_t * canon){ 
    return distance_gradient(p) * (distance(p) - canon->distance(p));
}

point_t
sdf_t::level_set_energy(point_t p, float epsilon){
    matrix_t h  = matrix_t::hessian(*phi, p);
    point_t g   = distance_gradient(p);
    float alpha = (g.length() - 1) / (g.length() + epsilon);

    return h * g * alpha;
}

point_t
sdf_t::killing_energy(point_t p, float gamma){ 
    matrix_t j = matrix_t::jacobian(*psi, p);
    std::vector<float> j_v = j.stack();
    std::vector<float> jt_v = j.transpose().stack();

    std::vector<float> v;
    for (int i = 0; i < j_v.size(); i++){
	v.push_back(jt_v[i] + j_v[i] * gamma);
    }

    matrix_t h_u = matrix_t::hessian(*psi_u, p);
    matrix_t h_v = matrix_t::hessian(*psi_v, p);
    matrix_t h_w = matrix_t::hessian(*psi_w, p);
    matrix_t h[3] = { h_u, h_v, h_w };    
    
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
