#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "canon_sdf.h"

ctpl::thread_pool sdf_t::pool(8);
std::vector<std::future<void>> sdf_t::futures;
bool sdf_t::is_initialised = false;
sdf_t::deform_field_t sdf_t::deform_field;

void
sdf_t::pool_wait(){
    for (int i = 0; i < futures.size(); i++){
        futures[i].get();
    }
    futures.clear();
}

sdf_t::sdf_t(depth_map_t depths, min_params_t * ps){
    this->depths = depths;
    this->ps = ps;   
    ps->camera_cx = depths->size() / 2;
    ps->camera_cy = depths->at(0).size() / 2;
 
    float l = ps->voxel_length;
    point_t size = ps->size;

    if (!is_initialised){
        std::cout << "Initialising deformation field..." << std::endl; 

        deform_field = std::vector<std::vector<std::vector<point_t>>>(
            size.get(0) / l,
            std::vector<std::vector<point_t>>(
                size.get(1) / l,
                std::vector<point_t>(
                    size.get(2) / l, 
                    point_t()
                )
            )
        );
 
        is_initialised = true;
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
    // deform
    point_t v = p + deformation_at(p);

    // project point
    float x;
    float y;
    project(v, &x, &y);
 
    // in case not in frame
    // TODO: not 100% sure this is the correct way to handle this case
    if (x < 0 || y < 0 || x >= depths->size() || y >= depths->at(0).size()){
	return 1;
    }

    // true signed distance
    int map = depths->at(x).at(y);
    if (map == 0){
        return 1;
    }
    float phi_true = map - v.get(2);

    // divide by delta
    float d = phi_true / ps->delta;
    
    // clamp to range [-1..1]
    float result = d / std::max(1.0f, std::abs(d));
 
    return result;
}

void
sdf_t::project(point_t p, float * x, float * y){
    // centre on origin 
    float rx = p.get(0) - ps->size.get(0) / 2;
    float ry = p.get(1) - ps->size.get(1) / 2;

    // perspective projection
    rx *= ps->camera_fx / (p.get(2) + ps->epsilon);
    ry *= ps->camera_fy / (p.get(2) + ps->epsilon);

    // re-centre in image
    *x = rx + ps->camera_cx;
    *y = ry + ps->camera_cy;
}

point_t
sdf_t::deformation_at(point_t p){
    // align to grid
    point_t v = p / ps->voxel_length;

    // clamp into range of volume
    // deformation lookups only happen outside volume when doing differentation
    // in this case, clamp into volume, reducing derivative to zero.
    int x = std::max(v.get(0), 0.0f);
    int y = std::max(v.get(1), 0.0f);
    int z = std::max(v.get(2), 0.0f);

    x = std::min(x, (int) deform_field.size() - 1);
    y = std::min(y, (int) deform_field[0].size() - 1);
    z = std::min(z, (int) deform_field[0][0].size() - 1);

    return deform_field[x][y][z];
}

point_t
sdf_t::distance_gradient(point_t p){
    return point_t(
        phi->differentiate(0)(p),
        phi->differentiate(1)(p),
        phi->differentiate(2)(p)
    );
}

void 
sdf_t::fuse(canon_sdf_t * canon){
    // rigid component
    bool should_update = true;
    for (int i = 1; should_update; i++){
	std::cout << "Rigid transformation, iteration " << i << "..." << std::endl;

        should_update = false;
        update(true, &should_update, canon);
    }
    std::cout << "Rigid transformation converged." << std::endl;

    // non-rigid component
    should_update = true;
    for (int i = 1; should_update; i++){
        std::cout << "Non-rigid transformation, iteration " << i << "..." << std::endl;

        should_update = false;
        update(false, &should_update, canon);
    }
    std::cout << "Non-rigid transformation converged." << std::endl;
}

void
sdf_t::update(bool is_rigid, bool * cont, canon_sdf_t * canon){
    auto f = [=](int id, int x, int y, int z){
        point_t p = (point_t(x, y, z) + point_t(0.5f)) * ps->voxel_length;
        point_t e = is_rigid ? 
            data_energy(p, canon) :
            energy(p, canon, ps->omega_k, ps->omega_s, ps->gamma, ps->epsilon);
 
        point_t u = e * ps->eta;
        if (u.length() > ps->threshold) {
            *cont = true;
//            std::cout << "  - update length GT threshold: " << u.to_string() << ".length() = " << u.length() << std::endl;
        }
        deform_field[x][y][z] -= u;

        if (!deform_field[x][y][z].is_finite()){
            std::cout << "Error: deformation field has diverged: " << deform_field[x][y][z].to_string() << " at: " << p.to_string() << std::endl;
            throw -1;
        }
    };

    for (int x = 0; x < deform_field.size(); x++){
        for (int y = 0; y < deform_field[0].size(); y++){
            for (int z = 0; z < deform_field[0][0].size(); z++){
                if (ps->mode == fusion_mode::CPU_MULTITHREAD){
                    futures.push_back(pool.push(f, x, y, z));        
                } else {
	            f(0, x, y, z);
                }
            }
	}
    }
     
    // wait for all threads to finish
    pool_wait(); 
}

point_t
sdf_t::energy(point_t v, canon_sdf_t * c, float o_k, float o_s, float gamma, float eps){
     return 
         data_energy(v, c) +
         killing_energy(v, gamma) * o_k +
         level_set_energy(v, eps) * o_s;
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
