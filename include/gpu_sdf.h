#ifndef GPU_SDF_H
#define GPU_SDF_H

#include "sdf.h"

class gpu_sdf_t : public sdf_t {
public:
    gpu_sdf_t(depth_map_t depths, point_t size, float l);

    void update_rigid(bool * cont, canon_sdf_t * previous, min_params_t * ps) override;
    void update_nonrigid(bool * cont, canon_sdf_t * previous, min_params_t * ps) override;
};

#endif
