#include "fusion.h"

int main(){
    fusion_t f;

    // declare hyper-parameters for system
    min_params_t ps;
    ps.mode = fusion_mode::CPU_MULTITHREAD;
    
    // killing fusion paper part 4.2
    ps.epsilon = 0.00005f;
    ps.threshold = 0.1f;

    // killing fusion paper part 4.3
    ps.size = point_t(80);

    // killing fusion paper part 5
    ps.eta = 0.1f; // called alpha in the paper
    ps.omega_k = 0.5f;
    ps.omega_s = 0.2f;
    ps.voxel_length = 4;  
    ps.gamma = 0.1f;
   
    // killing fusion paper part 3.1 
    ps.delta = ps.voxel_length * 10;
   
    // found this in the dynfu code somewhere 
    ps.camera_fx = 525;
    ps.camera_fy = 525;

    // no citation for this value, just a guess!
    ps.sdf_eta = 0.05f; 
    
    f.fusion(&ps);
    
    return 0;
}
