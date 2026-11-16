#include "yams_doctest.h"
#include <correlationscax.h>

// Deviation correlations of axial compressor cascades (Lieblein, as given by Aungier).
TEST_CASE("correlations.m10_naca65")
{
    // m = 0.17 - 0.0333 x + 0.333 x^2, x = beta1 / 100, beta1 in degrees
    REQUIRE(yams::m10_naca65(0.) == doctest::Approx(0.17));
    REQUIRE(yams::m10_naca65(50.) == doctest::Approx(0.17 - 0.0333 * 0.5 + 0.333 * 0.25));
    REQUIRE(yams::m10_naca65(-50.) == doctest::Approx(0.17 + 0.0333 * 0.5 + 0.333 * 0.25));
}
