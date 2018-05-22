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
    sdf_t(depth_map_t depths, min_params_t * ps);
    ~sdf_t();

    // main fusion method
    void fuse(bool is_rigid, canon_sdf_t * canon);

    // signed distance function
    float distance(point_t p);
    float phi_true(point_t p);

private:
    // singleton thread pool
    static ctpl::thread_pool pool;
    static std::vector<std::future<void>> futures;

    // static deformation field
    static bool is_initialised;
    static deform_field_t deform_field;

    // private fields
    depth_map_t depths;
    min_params_t * ps;
 
    // differentiable functions (see function.h)
    function_t<float> * phi;
    function_t<point_t> * psi;
    function_t<float> * psi_u;
    function_t<float> * psi_v;
    function_t<float> * psi_w;

    // gradient descent functions

    // energy functions 
    point_t energy(point_t p, canon_sdf_t * c, float o_k, float o_s, float gamma, float eps);
    point_t data_energy(point_t p, canon_sdf_t * canon);
    point_t killing_energy(point_t p, float gamma);
    point_t level_set_energy(point_t p, float epsilon);
    
    // other private methods private methods
    point_t deformation_at(point_t p);
    point_t distance_gradient(point_t p);
    float interpolate1D(float a, float b, float alpha);
    float interpolate2D(float a, float b, float c, float d, float alpha, float beta);
    float interpolate3D(float * xs, float alpha, float beta, float gamma);
};

#endif
