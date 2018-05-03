#ifndef CANONICAL_SDF_H
#define CANONICAL_SDF_H

#include "sdf.h"
#include <string>

class canon_sdf_t {
public:
    // constructor and destructor
    canon_sdf_t(min_params_t * ps);
    ~canon_sdf_t();   

    float distance(point_t p);
    void add_sdf(sdf_t * new_sdf);

    void save_mesh(std::string filename);

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
        float samples[8];
        point_t vertices[8];
        cell_t(point_t p, float l, canon_sdf_t * sdf);
    };

    typedef std::vector<std::vector<std::vector<weight_sdf_t>>> sampled_sdf_t;
    typedef std::vector<triangle_t> mesh_t;
    
    // private fields
    sampled_sdf_t sdf;
    point_t size;
    float voxel_length;
    float eta;

    // private functions
    float weight(float phi_true);
    void create_mesh(mesh_t * mesh);
    void create_mesh_for_cell(mesh_t * mesh, cell_t * cell);
};

#endif
