#pragma once
// Common include of the yams C++ tests (doctest).
#include <doctest/doctest.h>
#include <cmath>
#include <filesystem>
#include <limits>
#include <string>

// Interactive VTK plots in the tests are compiled only with the CMake option
// YAMS_TESTS_PLOT=ON (code under #if YAMS_TESTS_PLOT), so that ctest runs headless.
#ifndef YAMS_TESTS_PLOT
#define YAMS_TESTS_PLOT 0
#endif

namespace yams_tests
{
#ifdef TEST_DIR
    // Directory of the input grids of the tests (tests/data)
    inline std::string data_dir() { return std::string(TEST_DIR) + "/data/"; }

    // True when an input file is missing from tests/data: the test case is skipped
    // (doctest::skip) instead of failing on a file that could not be regenerated.
    inline bool data_missing(const std::string &name) { return !std::filesystem::exists(data_dir() + name); }
#endif

    // Equivalent of GoogleTest's ASSERT_DOUBLE_EQ: equal within 4 machine epsilons.
    inline doctest::Approx approx_ulp(double value)
    {
        return doctest::Approx(value).epsilon(4 * std::numeric_limits<double>::epsilon());
    }
}
