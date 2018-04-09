#ifndef FUSION_H
#define FUSION_H

#include "canonical_sdf.h"
#include <string>
#include <vector>

class fusion_t {
private:
    canonical_sdf_t canon;

    void load_filenames(std::vector<std::string> * fns);
    sdf_t * get_sdf(std::string filename);

public:
    // constructors and destructors
    fusion_t();
    ~fusion_t();
    
    // main public method
    void fusion();
};

#endif
