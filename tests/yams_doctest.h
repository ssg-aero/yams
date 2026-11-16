#pragma once
// Common include of the yams C++ tests (doctest).
#include <doctest/doctest.h>
#include <cmath>
#include <limits>

// Interactive VTK plots in the tests are compiled only with the CMake option
// YAMS_TESTS_PLOT=ON (code under #if YAMS_TESTS_PLOT), so that ctest runs headless.
#ifndef YAMS_TESTS_PLOT
#define YAMS_TESTS_PLOT 0
#endif

namespace yams_tests
{
    // Equivalent of GoogleTest's ASSERT_DOUBLE_EQ: equal within 4 machine epsilons.
    inline doctest::Approx approx_ulp(double value)
    {
        return doctest::Approx(value).epsilon(4 * std::numeric_limits<double>::epsilon());
    }
}
