#include "fusion.h"

void
fusion_t::fusion(){
    min_params_t ps;

    while (true){
        sdf_t * sdf = get_next_sdf();
        
        sdf->fuse(canon, &ps);
    }
}
