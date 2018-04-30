#include "canon_sdf.h"

#include <iostream>
#include <algorithm>

canon_sdf_t::canon_sdf_t(min_params_t * ps){
    voxel_length = ps->voxel_length;
    eta = ps->sdf_eta;    

    weight_sdf_t w; 

    for (int x = 0; x * voxel_length < ps->size.get(0); x++){
        sdf.push_back(std::vector<std::vector<weight_sdf_t>>());

        for (int y = 0; y * voxel_length < ps->size.get(1); y++){
            sdf[x].push_back(std::vector<weight_sdf_t>()); 

            for (int z = 0; z * voxel_length < ps->size.get(2); z++){
                sdf[x][y].push_back(w);
            }
        }
    }
}

canon_sdf_t::~canon_sdf_t(){

}

float 
canon_sdf_t::distance(point_t p){
    int x = p.get(0) / voxel_length;
    int y = p.get(1) / voxel_length;
    int z = p.get(2) / voxel_length;

    // TODO: not 100% these results are correct to return?
    if (x < 0 || y < 0 || z < 0 || x >= sdf.size() || y >= sdf[0].size() || z >= sdf[0][0].size()){
        return 1;
    } else if (sdf[x][y][z].omega == 0.0f){
        return 1;
    } else {
        return sdf[x][y][z].phi / sdf[x][y][z].omega;
    }
}

void
canon_sdf_t::add_sdf(sdf_t * new_sdf){
    for (int x = 0; x < sdf.size(); x++){
        for (int y = 0; y < sdf[0].size(); y++){
            for (int z = 0; z < sdf[0][0].size(); z++){
                point_t p = (point_t(x, y, z) + point_t(0.5f)) * voxel_length;
                float phi = new_sdf->distance(p);
                sdf[x][y][z].phi += phi;
                sdf[x][y][z].omega += weight(phi);
            }
        }
    } 
}

float
canon_sdf_t::weight(float phi_true){
    return phi_true > -eta ? 1 : 0;
}
