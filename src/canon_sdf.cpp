#include "canon_sdf.h"

#include <iostream>
#include <algorithm>
#include <fstream>

canon_sdf_t::canon_sdf_t(min_params_t * ps){
    voxel_length = ps->voxel_length;
    eta = ps->sdf_eta;    
    size = ps->size;

    for (int x = 0; x * voxel_length < size.get(0); x++){
        sdf.push_back(std::vector<std::vector<weight_sdf_t>>());

        for (int y = 0; y * voxel_length < size.get(1); y++){
            sdf[x].push_back(std::vector<weight_sdf_t>()); 

            for (int z = 0; z * voxel_length < size.get(2); z++){
                sdf[x][y].push_back(weight_sdf_t());
            }
        }
    }
}

canon_sdf_t::~canon_sdf_t(){

}

float 
canon_sdf_t::distance(point_t p){
    int x = p.get(0) / voxel_length;
    int y = p.get(1) / voxel_length;
    int z = p.get(2) / voxel_length;

    // TODO: not 100% these results are correct to return?
    if (x < 0 || y < 0 || z < 0 || x >= sdf.size() || y >= sdf[0].size() || z >= sdf[0][0].size()){
        return 1;
    } else if (sdf[x][y][z].omega == 0.0f){
        return sdf[x][y][z].phi;
    } else {
        return sdf[x][y][z].phi / sdf[x][y][z].omega;
    }
}

void
canon_sdf_t::add_sdf(sdf_t * new_sdf){
    for (int x = 0; x < sdf.size(); x++){
        for (int y = 0; y < sdf[0].size(); y++){
            for (int z = 0; z < sdf[0][0].size(); z++){
                point_t p = (point_t(x, y, z) + point_t(0.5f)) * voxel_length;
                float phi = new_sdf->distance(p);
                sdf[x][y][z].phi += phi * weight(phi);
                sdf[x][y][z].omega += weight(phi);
            }
        }
    } 
}

float
canon_sdf_t::weight(float phi_true){
    return phi_true > -eta ? 1 : 0;
}

canon_sdf_t::triangle_t::triangle_t(point_t a, point_t b, point_t c){
    vertices[0] = a;
    vertices[1] = b;
    vertices[2] = c;
}

canon_sdf_t::cell_t::cell_t(point_t p, float l, canon_sdf_t * sdf){
    for (int i = 0; i < 8; i++){
        point_t v = p;
        if (i & 1) v += point_t(0, 0, l);
        if (i & 2) v += point_t(0, l, 0);
        if (i & 4) v += point_t(l, 0, 0);

        vertices[i] = v;
        samples[i] = sdf->distance(v);
    }
}

void
canon_sdf_t::create_mesh(mesh_t * mesh){
    for (int x = 0; x < size.get(0) - voxel_length; x += voxel_length){
        for (int y = 0; y < size.get(1) - voxel_length; y += voxel_length){
            for (int z = 0; z < size.get(2) - voxel_length; z += voxel_length){
                point_t p = point_t(x, y, z) + point_t(voxel_length / 2);
                cell_t cell(p, voxel_length, this);

                create_mesh_for_cell(mesh, &cell);
            }
        }
    }
}

point_t
canon_sdf_t::normal(point_t p){
    return point_t(
        //TODO
    );
}

void
canon_sdf_t::create_mesh_for_cell(mesh_t * mesh, cell_t * cell){
    // TODO: marching cubes
}

void
canon_sdf_t::save_mesh(std::string filename){
    std::cout << "Saving SDF to mesh..." << std::endl;   

    // full filename
    std::string full_name = "../data/mesh/" + filename;

    // create triangles
    mesh_t mesh;
    create_mesh(&mesh);

    // create default material file
    std::ofstream mat_file;
    mat_file.open(full_name + ".mtl");
    mat_file << 
    "# Material file for " << filename  << std::endl <<
    "newmtl " << filename << "Material" << std::endl <<
    "Ns 96.078431"                      << std::endl <<
    "Ka 1.000000 1.000000 1.000000"     << std::endl <<
    "Kd 0.640000 0.640000 0.640000"     << std::endl <<
    "Ks 0.500000 0.500000 0.500000"     << std::endl <<
    "Ke 0.000000 0.000000 0.000000"     << std::endl <<
    "Ni 1.000000"                       << std::endl <<
    "d 1.000000"                        << std::endl <<
    "illum 2"                           << std::endl;
    mat_file.close();

    // save to wavefront .obj format
    std::ofstream mesh_file;
    mesh_file.open(full_name + ".obj");
    mesh_file <<
    "# Geometry file for " << filename << std::endl <<
    "mtllib " << filename << ".mtl"    << std::endl <<   
    "o " << filename << "Object"       << std::endl;

    // vertices
    for (auto tri : mesh){
        for (int i = 0; i < 3; i++){
            mesh_file << "v ";
            for (int j = 0; j < 3; j++){
                mesh_file << tri.vertices[i].get(j) << " ";
            }
            mesh_file << std::endl;
        }
    }

    // normals 
    for (auto tri : mesh){
        for (int i = 0; i < 3; i++){
            mesh_file << "vn ";
            for (int j = 0; j < 3; j++){
                mesh_file << normal(tri.vertices[i]).get(j) << " ";
            }
            mesh_file << std::endl;
        }
    }

    mesh_file << 
    "usemtl " << filename << "Material" << std::endl <<
    "s off"                             << std::endl;          

    // faces
    for (int i = 0; i < mesh.size(); i += 3){
        mesh_file << "f " << std::endl;
        for (int j = 0; j < 3; j++){
            mesh_file << i+j << "//" << i+j << " ";
        }
        mesh_file << std::endl;
    }

    mesh_file.close(); 
}
