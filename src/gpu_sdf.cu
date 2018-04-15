#include "gpu_sdf.h"

#include <cuda.h>

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

gpu_sdf_t::gpu_sdf_t(depth_map_t depths, point_t size, float l) : sdf_t(depths, size, l, true){

}
