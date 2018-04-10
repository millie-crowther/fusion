#include "canon_sdf.h"

#include <iostream>

canon_sdf_t::canon_sdf_t(){

}

canon_sdf_t::~canon_sdf_t(){
    while (!sdfs.empty()){
        delete *sdfs.begin();
        sdfs.pop_front();
    }
}

float
canon_sdf_t::weight(int i){
    return 1;
}

float 
canon_sdf_t::distance(point_t p){
    //unweighted for now
    float r = 0;
    for (auto sdf : sdfs){
//        r += sdf->distance(p);
    }
    return r / sdfs.size();
}

void
canon_sdf_t::add_sdf(sdf_t * sdf){
    sdfs.push_back(sdf);

    if (sdfs.size() > 10){ // TODO: adjust this limit? (or remove entirely...?)
        delete *sdfs.begin();
        sdfs.pop_front();
    }
}
