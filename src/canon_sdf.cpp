#include "canon_sdf.h"

#include <iostream>
#include <algorithm>

canon_sdf_t::canon_sdf_t(min_params_t * ps){
    voxel_length = ps->voxel_length;
    size = ps->size;
    n = 0;

    weight_sdf_t w; 
    w.phi = 0;
    w.omega = 0;

    for (int x = 0; x * voxel_length < size.get(0); x++){
        sdf.push_back(std::vector<std::vector<weight_sdf_t>>());

        for (int y = 0; y * voxel_length < size.get(1); y++){
            sdf[x].push_back(std::vector<weight_sdf_t>()); 

            for (int z = 0; z * voxel_length < size.get(2); z++){
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

    if (x < 0 || y < 0 || z < 0 || x >= sdf.size() || y >= sdf[0].size() || z >= sdf[0][0].size()){
        return 1;
    } else {
        return sdf[x][y][z].phi;
    }
}

void
canon_sdf_t::add_sdf(sdf_t * new_sdf){
    float a = (float) n / (float)(n + 1);
    float b = 1.0f / (float) n;
    point_t p;
 
    for (int x = 0; x < sdf.size(); x++){
        for (int y = 0; y < sdf[0].size(); y++){
            for (int z = 0; z < sdf[0][0].size(); z++){
                p = (point_t(x, y, z) + point_t(0.5f)) * voxel_length;
                if (n == 0){
                    sdf[x][y][z].phi = new_sdf->distance(p);
                } else {
                    sdf[x][y][z].phi = a * sdf[x][y][z].phi + b * new_sdf->distance(p);  
                }
            }
        }
    } 

    n++;
}
