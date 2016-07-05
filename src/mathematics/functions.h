#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include "core/epcDebug.h"

namespace epc {

namespace fun {

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

} // fun

} // epc

#endif