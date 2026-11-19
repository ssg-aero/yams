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

// Slope of the reference deviation against the camber, NASA SP-36 eq. 269 and 280:
// m = m_sigma=1 / sigma^b, b(beta1) of Fig. 164 (Aungier's fit).
TEST_CASE("correlations.slope_param")
{
    auto b = [](double b1)
    {
        auto x = b1 / 100;
        return 0.9625 - 0.17 * x - 0.85 * x * x * x;
    };
    // Fig. 164: b = 0.965 at 0 deg, 0.555 at 70 deg (chart reading)
    REQUIRE(b(0.) == doctest::Approx(0.965).epsilon(0.01));
    REQUIRE(b(70.) == doctest::Approx(0.555).epsilon(0.01));

    for (auto b1 : {0., 30., 45., 70.})
    {
        auto m10 = yams::m10_naca65(b1);
        // unit solidity: the chart value
        REQUIRE(yams::slope_param(yams::m10_naca65, b1, 1.) == doctest::Approx(m10));
        // the slope decreases with the solidity (it was multiplied by sigma^b)
        REQUIRE(yams::slope_param(yams::m10_naca65, b1, 1.5) == doctest::Approx(m10 / std::pow(1.5, b(b1))));
        REQUIRE(yams::slope_param(yams::m10_naca65, b1, 2.) < m10);
        REQUIRE(yams::slope_param(yams::m10_naca65, b1, 0.5) > m10);
    }
}
