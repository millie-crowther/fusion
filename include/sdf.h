#ifndef SDF_H
#define SDF_H

#include <vector>
#include "min_params.h"
#include "ctpl_stl.h"

class canon_sdf_t;

class sdf_t {
public:
    // types
    typedef std::vector<std::vector<int>> * depth_map_t;
    typedef std::vector<std::vector<std::vector<point_t>>> deform_field_t;

    // constructor and destructor
    sdf_t(int id, depth_map_t depths, min_params_t * ps);
    ~sdf_t();

    // main fusion method
    void fuse(canon_sdf_t * canon);

    // signed distance function
    float distance(point_t p);
    float distance_undeformed(point_t p);
    float weight(point_t p);

private:
    // singleton thread pool
    static ctpl::thread_pool pool;
    static std::vector<std::future<void>> futures;

    // static deformation field
    static bool is_initialised;
    static deform_field_t deform_field;

    // private fields
    int id;
    depth_map_t depths;
    min_params_t * ps;

    // energy functions 
    point_t energy(point_t p, canon_sdf_t * c, float o_k, float o_s, float gamma, float eps);
    point_t data_energy(point_t p, canon_sdf_t * canon);
    point_t killing_energy(point_t p, float gamma);
    point_t level_set_energy(point_t p, float epsilon);
    
    // other private methods private methods
    float interpolate_sdf(point_t p);
    float phi_data(int x, int y, float z);
    float phi_true(point_t p);
    point_t deformation_at(point_t p);
    point_t distance_gradient(point_t p);
    float interpolate1D(float a, float b, float alpha);
    float interpolate2D(float a, float b, float c, float d, float alpha, float beta);
    float interpolate3D(float * vs, point_t alpha);
};

#endif
