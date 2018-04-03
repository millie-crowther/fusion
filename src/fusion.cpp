#include "fusion.h"

void
fusion_t::fusion(){
    min_params_t ps;

    // set defaults
    ps.eta = 0.1f;
    ps.omega_k = 0.5f;
    ps.omega_s = 0.2f;
    ps.gamma = 0.1f;
    ps.epsilon = 0.00005f;
    ps.threshold = 0.1f;

    while (true){
        sdf_t * sdf = get_next_sdf();
        
        sdf->fuse(canon, &ps);
    }
}
