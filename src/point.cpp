#include "point.h"

#include <math.h>

point_t::point_t() : point_t(0, 0, 0){

}

point_t::point_t(const point_t& o) : point_t(o.x, o.y, o.z){

}

point_t::point_t(float x, float y, float z){
    this->x = x;
    this->y = y;
    this->z = z;
}

float
point_t::length(){
    return sqrt(x * x + y * y + z * z); 
}

point_t
point_t::operator*(float scale){
    return point_t(x * scale, y * scale, z * scale);
}

point_t
point_t::operator/(float scale){
    return (*this) * (1.0f / scale);
}

point_t
point_t::operator+(point_t o){
    return point_t(x + o.x, y + o.y, z + o.z);
}

void
point_t::operator-=(point_t o){
    x -= o.x;
    y -= o.y;
    z -= o.z;
}

point_t
point_t::operator-(point_t o){
    return point_t(x - o.x, y - o.y, z - o.z);
}

float
point_t::get_x(){
    return x;
}

float
point_t::get_y(){
    return y;
}

float
point_t::get_z(){
    return z;
}

