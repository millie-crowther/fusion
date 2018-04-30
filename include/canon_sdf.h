#ifndef CANONICAL_SDF_H
#define CANONICAL_SDF_H

#include "sdf.h"


class canon_sdf_t {
public:
    // constructor and destructor
    canon_sdf_t(min_params_t * ps);
    ~canon_sdf_t();   

    float distance(point_t p);
    void add_sdf(sdf_t * new_sdf);

private:
    // types
    struct weight_sdf_t { float phi; float omega; };
    typedef std::vector<std::vector<std::vector<weight_sdf_t>>> sampled_sdf_t;
    
    // private fields
    int n;
    sampled_sdf_t sdf;
    float voxel_length;
    point_t size;

    // private functions
    float weight(float phi_true);
};

#endif
