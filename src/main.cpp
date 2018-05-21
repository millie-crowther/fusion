#include "fusion.h"

int main(){
    fusion_t f;

    // declare hyper-parameters for system
    min_params_t ps;
    ps.is_multithreaded = false;
    ps.frames = 550;
    
    // killing fusion paper part 4.2
    ps.epsilon = 0.00005f;
    ps.threshold = 0.1f;

    // killing fusion paper part 5
    ps.eta_nonrigid = 0.1f; // called alpha in the paper
    ps.omega_k = 0.5f;
    ps.omega_s = 0.2f;
    ps.gamma = 0.1f;
   
    // empirically determined
    ps.eta_rigid = 0.2f;
    ps.sdf_eta = 100;
    ps.size = point_t(640, 480, 2048);
    ps.voxel_length = 16; 
    
    // killing fusion paper part 3.1 
    ps.delta = ps.voxel_length * 10;
   
    f.fusion(&ps);
    
    return 0;
}
