#ifndef GPU_SDF_H
#define GPU_SDF_H

#include "sdf.h"

class gpu_sdf_t : public sdf_t {
public:
    gpu_sdf_t(depth_map_t depths, min_params_t * ps);
    ~gpu_sdf_t();

    void update(bool is_rigid, bool * cont, canon_sdf_t * canon) override;
};

#endif
