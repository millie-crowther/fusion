#include "canon_sdf.h"

#include <iostream>
#include <algorithm>

canon_sdf_t::canon_sdf_t(min_params_t * ps){
    voxel_length = ps->voxel_length;
    size = ps->size;

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
   return 0;
}

void
canon_sdf_t::add_sdf(sdf_t * sdf){
}
