#ifndef CANONICAL_SDF_H
#define CANONICAL_SDF_H

#include "sdf.h"

typedef std::vector<std::vector<std::vector<float>>> sampled_sdf_t;

class canon_sdf_t {
public:
    // constructor and destructor
    canon_sdf_t(min_params_t * ps);
    ~canon_sdf_t();   

    float distance(point_t p);
    void add_sdf(sdf_t * sdf);

private:
    sampled_sdf_t sdf;
    float voxel_length;
    point_t size;
};

#endif
