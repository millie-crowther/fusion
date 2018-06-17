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
fusion_t::get_sdf(int frame, min_params_t * ps){
    std::string filename = get_filename(ps->dataset, frame, true, false);

    std::cout << "loading depth map: " << filename << std::endl; 

    std::string omaskname = get_filename(ps->dataset, frame, true, true);

    CImg<unsigned short> image(filename.c_str());
    CImg<unsigned short> omask(omaskname.c_str());

    sdf_t::depth_map_t depths = new std::vector<std::vector<int>>();
    for (int x = 0; x < image.width(); x++){
        depths->push_back(std::vector<int>());
        
        for (int y = 0; y < image.height(); y++){
            int d = *image.data(x, y, 0, 0);
	    int flag = *omask.data(x, y, 0, 0);

	    if (flag == 255 && d != 0){
                depths->at(x).push_back(d);
	    } else {
                depths->at(x).push_back((int) ps->size.z + 600);
	    }
        }
    }

    return sdf_t(frame, depths, ps);
}

std::string
fusion_t::get_filename(std::string dataset, int frame, bool is_kf_format, bool is_omask){
    std::string pad = frame < 10 ? "0" : "";
    if (frame < 100){
        pad = "0" + pad;
    }

    if (is_omask){
        return "../data/" + dataset + "/omask/omask_000" + pad + std::to_string(frame) + ".png";
    }

    if (is_kf_format){
        return "../data/" + dataset + "/depth/depth_000" + pad + std::to_string(frame) + ".png";
    } else {
        return "../data/" + dataset + "/depth/frame-000" + pad + 
	       std::to_string(frame) + ".depth.png";
    }
}

void
fusion_t::fusion(min_params_t * ps){

    canon = new canon_sdf_t(ps); 
    sdf_t initial = get_sdf(0, ps);
    canon->add_sdf(&initial);

    auto phi_global = [&](point_t p){ return canon->distance(p); };

    canon_sdf_t::save_mesh(phi_global, ps->dataset, 0);

    auto start = std::chrono::system_clock::now();
    for (int i = 1; i < ps->frames - 1; i++){
        std::cout << "Frame number: " << i << std::endl;     

        sdf_t sdf = get_sdf(i, ps);

        // rigid component
//        std::cout << "Calculating rigid deformation..." << std::endl;
  //      sdf.fuse(true, canon);
    //    std::cout << "Rigid deformation converged." << std::endl;
 
        // non-rigid component
        std::cout << "Calculating non-rigid deformation..." << std::endl;
	bool cont = true;
	//for (int i = 1; cont && i < ps->max_iterations; i++){
	  //  std::cout << "Iteration " << i << "..." << std::endl;
	   // cont = false;
            sdf.fuse(false, canon, &cont);
	//}
        std::cout << "Non-rigid deformation converged." << std::endl;
        canon->add_sdf(&sdf);

        if (i % 25 == 0){
            canon_sdf_t::save_mesh(phi_global, ps->dataset, i);
            auto phi_null   = [&](point_t p){ return sdf.distance_undeformed(p); };
            auto phi_n      = [&](point_t p){ return sdf.distance(p); };
            canon_sdf_t::save_mesh(phi_null, ps->dataset, i+1);
            canon_sdf_t::save_mesh(phi_n, ps->dataset, i+2);
        }

        std::cout << std::endl;
    } 
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<float> elapsed_seconds = end - start;    
    float t = elapsed_seconds.count();
    std::cout << "Total time elapsed: " << t << " seconds." << std::endl;
    std::cout << "Average framerate: " << ps->frames / t << " frames per second." << std::endl;

    canon_sdf_t::save_mesh(phi_global, ps->dataset, ps->frames);
}
