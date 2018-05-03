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
    delete canon;
}

sdf_t 
fusion_t::get_sdf(std::string filename, min_params_t * ps){
    std::cout << "loading depth map: " << filename << std::endl; 

    CImg<unsigned char> image(filename.c_str());

    sdf_t::depth_map_t depths = new std::vector<std::vector<unsigned char>>(image.width());
    for (int x = 0; x < image.width(); x++){
        depths->push_back(std::vector<unsigned char>(image.height()));
        
        for (int y = 0; y < image.height(); y++){
            depths->at(x).push_back(*image.data(x, y, 0, 0));
        }
    }

    if (ps->mode == fusion_mode::GPU){
        return gpu_sdf_t(depths, ps);
    } else {
        return sdf_t(depths, ps);
    }
}

void
fusion_t::load_filenames(std::vector<std::string> * fns){
    for (int i = 0; i < 25; i++){
        std::string padding = i < 10 ? "0" : "";
        if (i < 100){
            padding = "0" + padding;
        }
        fns->push_back("../data/umbrella/depth/frame-000" + padding + std::to_string(i) + ".depth.png");
    }
}

void
fusion_t::fusion(min_params_t * ps){
    std::vector<std::string> filenames;
    load_filenames(&filenames);

    canon = new canon_sdf_t(ps); 
    sdf_t initial = get_sdf(filenames[0], ps);
    canon->add_sdf(&initial);

    for (int i = 0; i < filenames.size(); i++){
        std::cout << "Frame number: " << i + 1 << std::endl;     

        sdf_t sdf = get_sdf(filenames[i], ps);
        sdf.fuse(canon);
        canon->add_sdf(&sdf);

        std::cout << std::endl;
    }

    canon->save_mesh("test");
}
