#ifndef FUSION_H
#define FUSION_H

#include "sdf.h"

class fusion_t {
private:
    sdf_t * canon;

    sdf_t * get_next_sdf();

public:
    void fusion();
};

#endif
