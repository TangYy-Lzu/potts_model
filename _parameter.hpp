#ifndef parameter_hpp
#define parameter_hpp
#include <math.h>

namespace Parameters
{
    // Magnetic coupling constant
    constexpr double J = 1;

    // Quantum number
    constexpr int Q = 2;

    // Type of lattice
    constexpr const char *LATTICE_TYPE = "square";

    // Size of lattice
    constexpr int LATTICE_SIZE = 6; // Should be a multiple of 3 for diced lattice and 2 for union jack lattice

    // Maximum and minimum temperature
    constexpr double T_MAX = 3.0;
    constexpr double T_MIN = 0.0;

    // Critical temperatures
    constexpr double T_CRIT_UP = 1.5;
    constexpr double T_CRIT_DOWN = 1.0;

    // Number of samples and groups
    constexpr int BROKEN_SIZE = 8;
    constexpr int SAMPLE_SIZE = 100;
    constexpr int GROUP_SIZE = 100;

    // Time steps
    constexpr double DELTA_T = 0.1;
    constexpr double DELTA_T_CRIT = 0.01; // For precise calculation near phase transition
}

#endif // PARAMETER_H_INCLUDE
