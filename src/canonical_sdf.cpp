#include "canonical_sdf.h"

canonical_sdf_t::canonical_sdf_t(){

}

canonical_sdf_t::~canonical_sdf_t(){
    while (!sdfs.empty()){
        delete (*sdfs.begin());
        sdfs.pop_front();
    }
}

float
canonical_sdf_t::weight(int i){
    return 1;
}

float 
canonical_sdf_t::distance(point_t p){
    //unweighted for now
    float r = 0;
    for (auto sdf : sdfs){
        r += sdf->distance(p);
    }
    return r / sdfs.size();
}

void
canonical_sdf_t::add_sdf(sdf_t * sdf){
    sdfs.push_back(sdf);

    if (sdfs.size() > 50){ // TODO: adjust this limit? (or remove entirely...?)
        delete (*sdfs.begin());
        sdfs.pop_front();
    }
}
