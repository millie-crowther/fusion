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
gpu_sdf_t::gpu_sdf_t(depth_map_t depths, min_params_t * ps) : sdf_t(depths, ps){

}

gpu_sdf_t::~gpu_sdf_t(){

}

void
gpu_sdf_t::update(bool is_rigid, bool * cont, canon_sdf_t * canon){

}
