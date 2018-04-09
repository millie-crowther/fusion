#ifndef CANONICAL_SDF_H
#define CANONICAL_SDF_H

#include "sdf.h"

#include <list>

class canonical_sdf_t {
public:
    // constructor and destructor
    canonical_sdf_t();
    ~canonical_sdf_t();   

    float distance(point_t p);

    void add_sdf(sdf_t * sdf);

private:
    std::list<sdf_t *> sdfs;

    float weight(int i);
};

#endif
