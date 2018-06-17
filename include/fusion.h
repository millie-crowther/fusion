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

    sdf_t get_sdf(int frame, min_params_t * ps);
    std::string get_filename(std::string dataset, int frame, bool is_kf_format, bool is_omask);
};

#endif
