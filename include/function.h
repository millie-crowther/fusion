#ifndef FUNCTION_H
#define FUNCTION_H

#include <functional>
#include "point.h"
#include <iostream>

template<class T>
class function_t {
private:
    float l;
    std::function<T(point_t)> f;

public:
    function_t(float l, std::function<T(point_t)> f){
        this->f = f;
        this->l = l;
    }

    T operator()(point_t p){
        return f(p);
    }

    function_t differentiate(int axis){
        point_t axes[3] = {
            point_t(l, 0, 0),
            point_t(0, l, 0),
            point_t(0, 0, l)
        };
        point_t u = axes[axis];

	std::function<T(point_t)> f_grad = [=](point_t x){
            return (this->f(x + u) - this->f(x - u)) / (2 * this->l);
        };

        return function_t(l, f_grad);
    }
};

#endif
