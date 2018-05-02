#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include "point.h"
#include <vector>
#include "canon_sdf.h"

struct triangle_t {
    point_t vertices[3]; 
    triangle_t(point_t a, point_t b, point_t c);
};

namespace marching {
    typedef std::vector<triangle_t> mesh_t;
    struct cell_t {
        float samples[8];
	point_t vertices[8];
	cell_t(point_t p, float l, canon_sdf_t * sdf);
    };

    void create_mesh(mesh_t * mesh, canon_sdf_t * sdf, point_t size, float voxel_length);
    void create_mesh_for_cell(mesh_t * mesh, cell_t * cell);    
}

#endif
