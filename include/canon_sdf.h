#ifndef CANONICAL_SDF_H
#define CANONICAL_SDF_H

#include "sdf.h"
#include <string>
#include <functional>

class canon_sdf_t {
public:
    // constructor and destructor
    canon_sdf_t(min_params_t * ps);
    ~canon_sdf_t();   

    float distance(point_t p);
    void add_sdf(sdf_t * new_sdf);

    static void save_mesh(std::function<float(point_t)> f, std::string model_name, int frame);

private:
    // types
    struct weight_sdf_t { 
        float phi; 
        float omega; 
        weight_sdf_t(){ phi = 0; omega = 0; }
    };
    
    struct triangle_t {
        point_t vertices[3];
        triangle_t(point_t a, point_t b, point_t c);
    };

    struct cell_t {
        float value[8];
        point_t point[8];
        cell_t(point_t p, float l, std::function<float(point_t)> f);
    };

    typedef std::vector<std::vector<std::vector<weight_sdf_t>>> sampled_sdf_t;
    typedef std::vector<triangle_t> mesh_t;
    
    // private fields
    sampled_sdf_t sdf;
    static point_t size;
    static float voxel_length;
    float eta;

    // private functions
    static point_t normal(std::function<float(point_t)> f, point_t p);
    static point_t interpolate(float isolevel, point_t a, point_t b, float alpha, float beta);
    static void create_mesh(std::function<float(point_t)> f, float isolevel, mesh_t * mesh);
    static void create_mesh_for_cell(float isolevel, mesh_t * mesh, cell_t * cell);
};

#endif
