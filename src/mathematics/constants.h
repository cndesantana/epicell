#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <limits>

namespace epc {

namespace mathConstants {
    const double pi = atan(1.0)*4.0;

    const std::vector<double> HexagonAngles = { 0.0, pi/3.0, 2.0*pi/3.0, pi, 4.0*pi/3.0,  5.0*pi/3.0};

    const std::vector<double> SquareAngles = { pi/4.0, 3.0*pi/4.0, 5.0*pi/4.0, 7.0*pi/4.0};

    static constexpr double epsilon = std::numeric_limits<double>::epsilon()*100.0;

    unsigned long SEED = 1;

    const double ForceEpsilon = 1.0e-3;

} // mathConstants

}  // end namespace epc

#endif
