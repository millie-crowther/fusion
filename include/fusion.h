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
    void fusion(min_params_t * ps);

private:
    canon_sdf_t * canon;

    void load_filenames(std::vector<std::string> * fns, int frames);
    sdf_t get_sdf(int id, std::string filename, min_params_t * ps);
};

#endif
