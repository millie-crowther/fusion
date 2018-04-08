#include "fusion.h"

#include "min_params.h"

#include "png.h"
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

fusion_t::fusion_t(){
    
}

fusion_t::~fusion_t(){

}

sdf_t *
fusion_t::get_sdf(std::string filename){
    std::cout << "loading depth map: " << filename << std::endl; 
 
    int x, y;

    int width, height;
    png_byte color_type;
    png_byte bit_depth;

    png_infop info_ptr;
    int number_of_passes;
    png_bytep * row_pointers;

    char header[8];
    
    // open file
    FILE * fp = fopen(filename.c_str(), "rb");
    if (!fp){
        std::cout << "failed to open file: \"" << filename << "\" for reading!" << std::endl;
        return nullptr;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    return nullptr; 
}

void
fusion_t::load_filenames(std::vector<std::string> * fns){
    for (int i = 0; i < 551; i++){
        std::string padding = i < 10 ? "0" : "";
        if (i < 100){
            padding = "0" + padding;
        }
        fns->push_back("../data/frame-000" + padding + std::to_string(i) + ".depth.png");
    }
}

void
fusion_t::fusion(){
    min_params_t ps;

    // set defaults
    ps.eta = 0.1f;
    ps.omega_k = 0.5f;
    ps.omega_s = 0.2f;
    ps.gamma = 0.1f;
    ps.epsilon = 0.00005f;
    ps.threshold = 0.1f;

    std::vector<std::string> filenames;
    load_filenames(&filenames);
 
    for (int i = 0; i < filenames.size(); i++){
        sdf_t * sdf = get_sdf(filenames[i]);
    }
}
