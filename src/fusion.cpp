#include "fusion.h"

#include "min_params.h"
#include <iostream>
#include "point.h"
#include "gpu_sdf.h"

#include "CImg.h"
using namespace cimg_library;

fusion_t::fusion_t(){
    
}

fusion_t::~fusion_t(){

}

sdf_t *
fusion_t::get_sdf(std::string filename, int mode){
    std::cout << "loading depth map: " << filename << std::endl; 

    CImg<unsigned char> image(filename.c_str());

    depth_map_t depths = new std::vector<std::vector<unsigned char>>(image.width());
    for (int x = 0; x < image.width(); x++){
        depths->push_back(std::vector<unsigned char>(image.height()));
        
        for (int y = 0; y < image.height(); y++){
            depths->at(x).push_back(*image.data(x, y, 0, 0));
        }
    }

    sdf_t * sdf = nullptr;

    if (mode == fusion_mode::GPU){
        sdf = new gpu_sdf_t(depths);
    } else {
        sdf = new sdf_t(depths, mode == fusion_mode::CPU_MULTITHREAD);
    }

    return sdf; 
}

void
fusion_t::load_filenames(std::vector<std::string> * fns){
    for (int i = 0; i < 551; i++){
        std::string padding = i < 10 ? "0" : "";
        if (i < 100){
            padding = "0" + padding;
        }
        fns->push_back("../data/umbrella/depth/frame-000" + padding + std::to_string(i) + ".depth.png");
    }
}

void
fusion_t::fusion(int mode){
    min_params_t ps;

    // set defaults
    ps.eta_rigid = 0.1f;
    ps.eta_nonrigid = 0.1f;
    ps.omega_k = 0.5f;
    ps.omega_s = 0.2f;
    ps.gamma = 0.1f;
    ps.epsilon = 0.00005f;
    ps.threshold_rigid = 0.1f;
    ps.threshold_nonrigid = 0.1f;

    std::vector<std::string> filenames;
    load_filenames(&filenames);
 
    sdf_t * initial = get_sdf(filenames[0], fusion_mode::NIL);
    canon.add_sdf(initial);

    sdf_t * previous = initial;

    for (int i = 1; i < filenames.size(); i++){
        std::cout << "Frame number: " << i << std::endl;     

        sdf_t * sdf = get_sdf(filenames[i], mode);

        sdf->fuse(&canon, previous, &ps);

        canon.add_sdf(sdf);

        previous = sdf;

        std::cout << std::endl;
    }
}
