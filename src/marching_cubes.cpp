#include "marching_cubes.h"

triangle_t::triangle_t(point_t a, point_t b, point_t c){
    vertices[0] = a;
    vertices[1] = b;
    vertices[2] = c;
}

marching::cell_t::cell_t(point_t p, float l, canon_sdf_t * sdf){

}

void
marching::create_mesh(mesh_t * mesh, canon_sdf_t * sdf, point_t size, float voxel_length){
    for (int x = 0; x < size.get(0); x += voxel_length){
        for (int y = 0; y < size.get(1); y += voxel_length){
            for (int z = 0; z < size.get(2); z += voxel_length){
                cell_t cell(point_t(x, y, z), voxel_length, sdf);
        
		create_mesh_for_cell(mesh, &cell);
	    }
	}
    }	
}

void
marching::create_mesh_for_cell(mesh_t * mesh, cell_t * cell){

}
