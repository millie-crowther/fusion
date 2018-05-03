#ifndef POINT_H
#define POINT_H

#include <string>

class point_t {
public:
    // constructors
    point_t();
    point_t(float x);
    point_t(float x, float y, float z);
    point_t(const point_t& p);
    
    // operator overrides
    point_t operator*(float scale);
    point_t operator/(float scale) const;
    point_t operator+(point_t other);
    point_t operator-(point_t other);
    void operator+=(point_t other);
    void operator-=(point_t other);

    // public functions
    float length();
    float get(int i) const;
    bool is_finite();
    std::string to_string();

private:
    // data
    float elem[3];
};

#endif
