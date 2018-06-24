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
sdf_t::interpolate3D(float * vs, point_t alpha){
    float s = interpolate2D(vs[0], vs[1], vs[2], vs[3], alpha.x, alpha.y);
    float t = interpolate2D(vs[4], vs[5], vs[6], vs[7], alpha.x, alpha.y);
    return interpolate1D(s, t, alpha.z);
}

point_t
sdf_t::deformation_at(point_t p){
    point_t v = p / ps->voxel_length;
    int x = v.x;
    int y = v.y;
    int z = v.z;
   
    if (
        x < 0 || y < 0 || z < 0 ||
        x >= deform_field.size() - 1 ||
        y >= deform_field[0].size() - 1 ||
        z >= deform_field[0][0].size() - 1
    ){
        return point_t();
    }

    return deform_field[x][y][z];
}

float
sdf_t::phi_data(int x, int y, float z){
    x = std::max(0, std::min(x, (int) depths->size() - 1));
    y = std::max(0, std::min(y, (int) depths->at(x).size() - 1));

    int map = depths->at(x).at(y);
    return map - ps->near_clip;
}

float
sdf_t::distance_undeformed(point_t p){
    float map = interpolate_sdf(p);
    return std::abs(map - p.z + ps->sdf_eta / 2.0f) - ps->sdf_eta / 2.0f;
}

float
sdf_t::weight(point_t v){
    v += deformation_at(v);
    float d = interpolate_sdf(v) - v.z;
    if (d > -ps->sdf_eta){
	return 1;
    } else {
	return 0;
    }
}

float
sdf_t::interpolate_sdf(point_t v){
    point_t v1 = point_t((int) v.x, (int) v.y, (int) v.z);

    float vs[8];
    for (int i = 0; i < 8; i++){
        point_t v2 = v1;
	if (i & 1) v2.x += 1;
	if (i & 2) v2.y += 1;
	if (i & 4) v2.z += 1;

        vs[i] = phi_data(v2.x, v2.y, v2.z);
    }
   
    point_t alpha = v - v1; 
    return interpolate3D(vs, alpha); 
}

float
sdf_t::phi_true(point_t v){
    v += deformation_at(v);
    float depth = interpolate_sdf(v);
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
sdf_t::distance_gradient(point_t p){
    float l = ps->voxel_length;
    auto r = point_t(
        distance(p + point_t(l, 0, 0)) - distance(p - point_t(l, 0, 0)),
        distance(p + point_t(0, l, 0)) - distance(p - point_t(0, l, 0)),
        distance(p + point_t(0, 0, l)) - distance(p - point_t(0, 0, l))
    ) / (2.0f * l);

    return r * ps->delta;
}

void
sdf_t::fuse(canon_sdf_t * canon){
    // anonymous function to be evaluated at each voxel
    auto f = [&](int id, int x, int y, int z){
        point_t p = (point_t(x, y, z) + point_t(0.5f)) * ps->voxel_length;
        bool quit = false;
        for (int i = 0; !quit && i < ps->max_iterations; i++){
            point_t e = energy(p, canon, ps->omega_k, ps->omega_s, ps->gamma, ps->epsilon);

            deform_field[x][y][z] -= e * ps->eta;
            
            if (glm::length(e) <= ps->threshold) {
                quit = true;
            } 

            // perform check on deformation field to see if it has diverged
	    point_t d = deform_field[x][y][z];
            if (!std::isfinite(d.x) || !std::isfinite(d.y) || !std::isfinite(d.z)){
                std::cout << "Error: deformation field has diverged: " 
          	          << glm::to_string(d) << " at: " << glm::to_string(p) << std::endl;
                throw -1;
            }
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
sdf_t::energy(
    point_t v, canon_sdf_t * c, 
    float o_k, float o_s, float gamma, float eps
){
    // function that calculates the three components of the energy gradient as 
    // outlined in the killing fusion paper
       	
    auto d = data_energy(v, c);
    auto k = killing_energy(v, gamma) * o_k;
    auto ls = level_set_energy(v, eps) * o_s;
    return d + k + ls;

}

point_t
sdf_t::data_energy(point_t p, canon_sdf_t * canon){
    // data component of energy gradient
    // measures how well aligned the two SDFs are
    point_t g = distance_gradient(p);
    if (glm::length(g) > 0){
	g = glm::normalize(g);
    }

    return g * (distance(p) - canon->distance(p));
}

point_t
sdf_t::level_set_energy(point_t p, float epsilon){
    // level set energy gradient component
    // requires the magnitude of the gradient of the SDF to be unity
    point_t g  = distance_gradient(p);
   
    auto phi = [&](point_t v){ return distance(v) * ps->delta; };
    matrix_t h  = matrix_t::hessian(phi, p, ps->voxel_length);

    float alpha = (glm::length(g) - 1) / (glm::length(g) + epsilon);

    return h * g * alpha;
}

point_t
sdf_t::killing_energy(point_t p, float gamma){
    // killing field energy gradient component
    // requires the deformation field to be a killing vector field

    auto psi = [&](point_t v){ return deformation_at(v); };
    auto E = [&](point_t v){
        matrix_t j = matrix_t::jacobian(psi, v, ps->voxel_length);

        std::vector<float> j_v = j.stack();
        std::vector<float> jt_v = j.transpose().stack();
        
	float result = 0;

        for (int i = 0; i < 9; i++){
            result += j_v[i] * j_v[i] + gamma * jt_v[i] * j_v[i];
	}

        // result clamped to upper bound to attempt to prevent divergence
	return std::min(result, 25.0f);
    };



    float l = ps->voxel_length;
    auto r = point_t(
        E(p + point_t(l, 0, 0)) - E(p - point_t(l, 0, 0)),
        E(p + point_t(0, l, 0)) - E(p - point_t(0, l, 0)),
        E(p + point_t(0, 0, l)) - E(p - point_t(0, 0, l))
    ) / (2.0f * l);


    return r;
    /*
      NOTE: previous version as outlined in the killingfusion paper

 //   matrix_t j = matrix_t::jacobian(psi, p, ps->voxel_length);
  //  std::vector<float> j_v = j.stack();
   // std::vector<float> jt_v = j.transpose().stack();
    std::vector<float> v;
    for (int i = 0; i < j_v.size(); i++){
	v.push_back(jt_v[i] + j_v[i] * gamma);
    }

    auto psi_u = [&](point_t v){ return deformation_at(v).x; };
    auto psi_v = [&](point_t v){ return deformation_at(v).y; };
    auto psi_w = [&](point_t v){ return deformation_at(v).z; };

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

    auto r = result * 2.0f;
    return r;*/
}
