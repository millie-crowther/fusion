#include "gpu_sdf.h"

//#include <cuda.h>

/*
__global__ void 
data_energy_kernel(float * x, int w, int h, int d){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= w * h * d){
        return;
    }
}

__global__ void 
killing_energy_kernel(float * x, int w, int h, int d){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= w * h * d){
        return;
    }
}

__global__ void 
level_set_energy_kernel(float * x, int w, int h, int d){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= w * h * d){
        return;
    }
}
*/
gpu_sdf_t::gpu_sdf_t(depth_map_t depths) : sdf_t(depths, true){

}

gpu_sdf_t::~gpu_sdf_t(){

}

void
gpu_sdf_t::update_rigid(bool * cont, canon_sdf_t * canon, min_params_t * ps){

}

void
gpu_sdf_t::update_nonrigid(bool * cont, canon_sdf_t * canon, min_params_t * ps){

}
