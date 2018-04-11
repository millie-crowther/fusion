#ifndef FUNCTION_H
#define FUNCTION_H

#include <functional>
#include "point.h"

class function_t {
private:
    float l;
    std::function<float(point_t)> f;

public:
    function_t(float l, std::function<float(point_t)>& f){
        this->f = f;
        this->l = l;
    }

    float operator()(point_t p){
        return f(p);
    }

    function_t differentiate(int axis){
        point_t axes[3] = {
            point_t(l, 0, 0),
            point_t(0, l, 0),
            point_t(0, 0, l)
        };
        point_t u = axes[axis];

        std::function<float(point_t)> f_grad = [=](point_t x){
            return (f(x + u) - f(x - u)) / (2 * l);
        };

        return function_t(l, f_grad);
    }
};

#endif
