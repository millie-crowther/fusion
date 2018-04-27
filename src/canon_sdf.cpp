#include "canon_sdf.h"

#include <iostream>
#include <algorithm>

canon_sdf_t::canon_sdf_t(min_params_t * ps){
    voxel_length = ps->voxel_length;
    size = ps->size;
    n = 0.0f;

    for (int x = 0; x * voxel_length < size.get(0); x++){
        sdf.push_back(std::vector<std::vector<float>>());

        for (int y = 0; y * voxel_length < size.get(1); y++){
            sdf[x].push_back(std::vector<float>()); 

            for (int z = 0; z * voxel_length < size.get(2); z++){
                sdf[x][y].push_back(0);
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

    if (x < 0 || y < 0 || z < 0 || x >= sdf.size() || y <= sdf[0].size() || z <= sdf[0][0].size()){
        return 0;
    } else {
        return sdf[x][y][z];
    }
}

void
canon_sdf_t::add_sdf(sdf_t * new_sdf){
    float a = n / (n + 1.0f);
    float b = 1.0f / n;

    for (int x = 0; x < sdf.size(); x++){
        for (int y = 0; y < sdf[0].size(); y++){
            for (int z = 0; z < sdf[0][0].size(); z++){
                point_t p = (point_t(x, y, z) + point_t(0.5f)) * voxel_length;
                sdf[x][y][z] = a * sdf[x][y][z] + b * new_sdf->distance(p);  
            }
        }
    } 

    n += 1.0f;
}
