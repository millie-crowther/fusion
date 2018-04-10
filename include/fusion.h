#ifndef FUSION_H
#define FUSION_H

#include "canon_sdf.h"
#include <string>
#include <vector>

class fusion_t {
public:
    // constructors and destructors
    fusion_t();
    ~fusion_t();
    
    // main public method
    void fusion();
private:
    canon_sdf_t canon;

    void load_filenames(std::vector<std::string> * fns);
    sdf_t * get_sdf(std::string filename);
};

#endif
