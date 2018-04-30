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

    // vector length
    float length();

    point_t hadamard(point_t o);

    // getters
    float get(int i) const;

    std::string to_string();
    bool is_finite();

private:
    float elem[3];
};

#endif
