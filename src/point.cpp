#include "point.h"

#include <math.h>

point_t::point_t() : point_t(0, 0, 0){

}

point_t::point_t(float x) : point_t(x, x, x){

}

point_t::point_t(const point_t& o) : point_t(o.elem[0], o.elem[1], o.elem[2]){

}

point_t::point_t(float x, float y, float z){
    elem[0] = x;
    elem[1] = y;
    elem[2] = z;
}

float
point_t::length(){
    return sqrt(elem[0] * elem[0] + elem[1] * elem[1] + elem[2] * elem[2]); 
}

point_t
point_t::operator*(float scale){
    return point_t(elem[0] * scale, elem[1] * scale, elem[2] * scale);
}

point_t
point_t::operator/(float scale) const {
    return point_t(elem[0] / scale, elem[1] / scale, elem[2] / scale);
}

point_t
point_t::operator+(point_t o){
    return point_t(elem[0] + o.elem[0], elem[1] + o.elem[1], elem[2] + o.elem[2]);
}

void
point_t::operator+=(point_t o){
    for (int i = 0; i < 3; i++){
	elem[i] += o.elem[i];
    }
}

void
point_t::operator-=(point_t o){
    for (int i = 0; i < 3; i++){
	elem[i] -= o.elem[i];
    }
}

point_t
point_t::operator-(point_t o){
    return point_t(elem[0] - o.elem[0], elem[1] - o.elem[1], elem[2] - o.elem[2]);
}

float
point_t::get(int i) const {
    return elem[i];
}

std::string
point_t::to_string(){
    return "point_t(" + 
	    std::to_string(get(0)) + ", " + 
	    std::to_string(get(1)) + ", " + 
	    std::to_string(get(2)) + 
    ")";
}
