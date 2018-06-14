#include "sdf.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "canon_sdf.h"
#include <glm/gtx/string_cast.hpp>
#include <atomic>
#include <thread>
#include <chrono>

ctpl::thread_pool sdf_t::pool(8);
std::vector<std::future<void>> sdf_t::futures;
bool sdf_t::is_initialised = false;
sdf_t::deform_field_t sdf_t::deform_field;

sdf_t::sdf_t(int id, depth_map_t depths, min_params_t * ps){
    this->id = id;
    this->depths = depths;
    this->ps = ps;   

    float l = ps->voxel_length;
    point_t size = ps->size;

    if (!is_initialised){
        std::cout << "Initialising deformation field..." << std::endl; 

        deform_field = std::vector<std::vector<std::vector<point_t>>>(
            size.x / l,
            std::vector<std::vector<point_t>>(
                size.y / l,
                std::vector<point_t>(
                    size.z / l, 
                    point_t()
                )
            )
        );

        is_initialised = true;
    }
}

sdf_t::~sdf_t(){
    delete depths;
}

float
sdf_t::interpolate1D(float a, float b, float alpha){
    return a * (1 - alpha) + b * alpha;
}

float
sdf_t::interpolate2D(float a, float b, float c, float d, float alpha, float beta){
    float s = interpolate1D(a, b, alpha);
    float t = interpolate1D(c, d, alpha);
    return interpolate1D(s, t, beta);
}

float
sdf_t::phi_data(int x, int y, float z){
    x = std::max(0, std::min(x, (int) depths->size() - 1));
    y = std::max(0, std::min(y, (int) depths->at(x).size() - 1));

    int map = depths->at(x).at(y);
    if (map == 0){
        return z + ps->delta;
    } else {
        return map - 500.0f;
    }
}

float
sdf_t::distance_undeformed(point_t p){
    int x = std::max(0, std::min((int) p.x, (int) depths->size() - 1));
    int y = std::max(0, std::min((int) p.y, (int) depths->at(x).size() - 1));
 
    int map = depths->at(x).at(y);
    if (map == 0){
        map = ps->size.z * 2.0f;
    }
    map -= 500.0f;

    return std::abs(map - p.z + ps->sdf_eta / 2.0f) - ps->sdf_eta / 2.0f;
}

float
sdf_t::weight(point_t p){
    p += deformation_at(p);

    int x = std::max(0, std::min((int) p.x, (int) depths->size() - 1));
    int y = std::max(0, std::min((int) p.y, (int) depths->at(x).size() - 1));
    if (depths->at(x).at(y) == 0){
        return 0;
    }       

    float map = phi_data(x, y, p.z);
    if (p.z > map + ps->sdf_eta){
        return 0;
    } else {
        return 1;
    }
}

float
sdf_t::phi_true(point_t v){
    v += deformation_at(v);

    float a = phi_data(v.x,        v.y,        v.z);   
    float b = phi_data(v.x + 1.0f, v.y,        v.z);   
    float c = phi_data(v.x,        v.y + 1.0f, v.z);   
    float d = phi_data(v.x + 1.0f, v.y + 1.0f, v.z);   

    float alpha = v.x - ((int) v.x);
    float beta  = v.y - ((int) v.y);
    
    float depth = interpolate2D(a, b, c, d, alpha, beta); 

    return std::abs(depth - v.z + ps->sdf_eta / 2.0f) - ps->sdf_eta / 2.0f;
}

float
sdf_t::distance(point_t p){
    float d = phi_true(p) / ps->delta;
    if (d > 1.0f) return 1.0f;
    if (d < -1.0f) return -1.0f;
    return d;
}

point_t
sdf_t::deformation_at(point_t p){
    // align to grid
    point_t v = p / ps->voxel_length;

    int x = v.x;
    int y = v.y;
    int z = v.z;
   
    if (
        x < 0 || y < 0 || z < 0 ||
        x >= deform_field.size() ||
        y >= deform_field[0].size() ||
        z >= deform_field[0][0].size()
    ){
        return point_t();
    } else {
        return deform_field[x][y][z];
    }
}

point_t
sdf_t::distance_gradient(point_t p){
    float l = 1.0f;
    auto r = point_t(
        phi_true(p + point_t(l, 0, 0)) - phi_true(p - point_t(l, 0, 0)),
        phi_true(p + point_t(0, l, 0)) - phi_true(p - point_t(0, l, 0)),
        phi_true(p + point_t(0, 0, l)) - phi_true(p - point_t(0, 0, l))
    ) / (2.0f * l);

    return r;
}

void
sdf_t::fuse(bool is_rigid, canon_sdf_t * canon){
    // anonymous function to be evaluated at each voxel
    auto f = [&](int id, int x, int y, int z){
        point_t p = (point_t(x, y, z) + point_t(0.5f)) * ps->voxel_length;

        int px = std::max(0, std::min((int)(p.x + deformation_at(p).x), (int) depths->size() - 1));
        int py = std::max(0, std::min((int)(p.y+deformation_at(p).y), (int) depths->at(0).size()-1));

        if (depths->at(px).at(py) == 0){
            return;
        }

        bool quit = false;
        for (int i = 0; !quit && i < ps->max_iterations; i++){
            point_t e = is_rigid ? 
                data_energy(p, canon) :
                energy(p, canon, ps->omega_k, ps->omega_s, ps->gamma, ps->epsilon);

            deform_field[x][y][z] -= e * ps->eta;
            
            if (glm::length(e) <= ps->threshold) {
                quit = true;
            } 

            // perform check on deformation field to see if it has diverged
	    point_t d = deform_field[x][y][z];
            if (!std::isfinite(d.x) || !std::isfinite(d.y) || !std::isfinite(d.z)){
                std::cout << "Error: deformation field has diverged: " 
          	          << glm::to_string(deform_field[x][y][z])
	  	          << " at: " << glm::to_string(p) << std::endl;
                throw -1;
            }

//	    float max_length = 500.0f;
//	    if (glm::length(d) > max_length){
//	        deform_field[x][y][z] *= max_length / glm::length(d);
//	    }
	}
    };

    // iterate over full volume, in serial or parallel as appropriate
    for (int x = 0; x < deform_field.size(); x++){
        for (int y = 0; y < deform_field[0].size(); y++){
            for (int z = 0; z < deform_field[0][0].size(); z++){
                if (ps->is_multithreaded){
                    futures.push_back(pool.push(f, x, y, z));
                } else {
	            f(0, x, y, z);
                }
            }
	}
    }
     
    // wait for all threads to finish
    for (int i = 0; ps->is_multithreaded && i < futures.size(); i++){
        futures[i].get();
    }
    futures.clear();

}

point_t
sdf_t::energy(point_t v, canon_sdf_t * c, float o_k, float o_s, float gamma, float eps){
    // function that calculates the three components of the energy gradient as 
    // outlined in the killing fusion paper
    return 
         data_energy(v, c) +
         killing_energy(v, gamma) * o_k +
         level_set_energy(v, eps) * o_s;
}

point_t
sdf_t::data_energy(point_t p, canon_sdf_t * canon){
    // data component of energy gradient
    // measures how well aligned the two SDFs are
    return distance_gradient(p) * (distance(p) - canon->distance(p));
}

point_t
sdf_t::level_set_energy(point_t p, float epsilon){
    // level set energy gradient component
    // requires the magnitude of the gradient of the SDF to be unity

    auto phi = [&](point_t v){ return phi_true(v); };
    matrix_t h  = matrix_t::hessian(phi, p, 1.0f);//ps->voxel_length);
    point_t g   = distance_gradient(p);
    float alpha = (glm::length(g) - 1) / (glm::length(g) + epsilon);

    return h * g * alpha;
}

point_t
sdf_t::killing_energy(point_t p, float gamma){
    // killing field energy gradient component
    // requires the deformation field to be a killing vector field

    auto psi = [&](point_t v){ return deformation_at(v); };

    matrix_t j = matrix_t::jacobian(psi, p, ps->voxel_length);
    std::vector<float> j_v = j.stack();
    std::vector<float> jt_v = j.transpose().stack();

    std::vector<float> v;
    for (int i = 0; i < j_v.size(); i++){
	v.push_back(jt_v[i] + j_v[i] * gamma);
    }

    auto psi_u = [&](point_t v){ return deformation_at(p).x; };
    auto psi_v = [&](point_t v){ return deformation_at(p).y; };
    auto psi_w = [&](point_t v){ return deformation_at(p).z; };

    matrix_t h[3] = {
        matrix_t::hessian(psi_u, p, ps->voxel_length),
        matrix_t::hessian(psi_v, p, ps->voxel_length),
        matrix_t::hessian(psi_w, p, ps->voxel_length)
    };

    point_t result;
    for (int i = 0; i < 9; i++){
	result += point_t(
            h[i / 3].get(i % 3, 0) * v[i],
            h[i / 3].get(i % 3, 1) * v[i],
            h[i / 3].get(i % 3, 2) * v[i]
	);
    }

    return result * 2.0f;
}
