#include "fusion.h"

#include "min_params.h"
#include <iostream>
#include <chrono>
#include <unistd.h>
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

    CImg<unsigned short> image(filename.c_str());

    sdf_t::depth_map_t depths = new std::vector<std::vector<int>>();
    for (int x = 0; x < image.width(); x++){
        depths->push_back(std::vector<int>());
        
        for (int y = 0; y < image.height(); y++){
            int d = *image.data(x, y, 0, 0);
            depths->at(x).push_back(d);
        }
    }

    return sdf_t(depths, ps);
}

void
fusion_t::load_filenames(std::vector<std::string> * fns, int frames){
    for (int i = 0; i < frames; i++){
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
    load_filenames(&filenames, ps->frames);

    canon = new canon_sdf_t(ps); 
    sdf_t initial = get_sdf(filenames[0], ps);
    canon->add_sdf(&initial);

    auto start = std::chrono::system_clock::now();
    for (int i = 1; i < filenames.size(); i++){
        std::cout << "Frame number: " << i << std::endl;     

        sdf_t sdf = get_sdf(filenames[i], ps);
        sdf.fuse(canon);
        canon->add_sdf(&sdf);

        if (i % 30 == 0){
            canon->save_mesh("umbrella", i);
        }

        std::cout << std::endl;
    } 
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<float> elapsed_seconds = end - start;    
    float t = elapsed_seconds.count();
    std::cout << "Total time elapsed: " << t << " seconds." << std::endl;
    std::cout << "Average framerate: " << ps->frames / t << " frames per second." << std::endl;

    canon->save_mesh("umbrella", ps->frames);
}
