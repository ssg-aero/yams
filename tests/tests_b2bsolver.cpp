#include "yams_doctest.h"
#include <chrono>
#include <iostream>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <baldeToBlade/bladeToBladeCurvatureSolverPost.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <gbs/bscbuild.h>
#include <string>
#include <gas/gasModel.h>
TEST_CASE("tridiagonal.thomas_vec")
{
        using namespace std;
        using T = double;
        vector<T> a{ 0, -1, -1, -1 };
        vector<T> b{ 4,  4,  4,  4 };
        vector<T> c{-1, -1, -1,  0 };
        vector<T> d{ 5,  5, 10, 23 };
        vector<T> r{ 2,  3,  5, 7  };
        thomas_algorithm(a, b, c, d);

        for(int i{}; i < 4 ;i++)
            REQUIRE(std::abs(d[i] - r[i]) <= 1e-9);
}

TEST_CASE("solver.base_novak")
{
    using namespace std;
    using T = double;

    auto test_dir = std::string(TEST_DIR);
    auto f_name = test_dir+std::string("/b2bmsh_medium.vts");
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_coarse.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_medium.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_finest.vts";

    auto [pts, nj]  = yams::getVtkStructuredGridPoints<T>(f_name.c_str());
    auto stream_line = make_shared<gbs::BSCurve<T,2>>( 
        gbs::build_segment<T,2>({0.,1.},{1.,1.})
    );
    // auto stream_line = make_shared<gbs::BSCurve<T,2>>( 
    //     gbs::BSCurve<T,2>(
    //         {
    //             {0.0,1.0},
    //             {0.5,1.1},
    //             {1.0,1.1},
    //         },
    //         {0.,1.},
    //         {3,3},
    //         2
    //     )
    // );
    // auto stream_line = make_shared<gbs::BSCurve<T,2>>( 
    //     gbs::BSCurve<T,2>(
    //         {
    //             {0.0,1.0},
    //             {0.1,1.0},
    //             {0.9,1.1},
    //             {1.0,1.1},
    //         },
    //         {0.,1.},
    //         {4,4},
    //         3
    //     )
    // );

    yams::BladeToBladeCurvatureSolver<T> solver{ 
        pts, 
        nj, 
        stream_line, 
        11,
        4, 15
        // 3, 11
    };

    solver.setRotationSpeed( -210 );

    solver.setFullBladePassageMassFlow(
        // 74.
        // 74.8 
        // 74.9251
        // 75
        // 138
        // 70
        // 60
        92.36282401553996
    );
    auto ni = solver.dimensions()[0];

    auto start = chrono::steady_clock::now();
    solver.computeW();
    // converged by the tolerances (1e-4), not stopped by the iteration limit
    CHECK(solver.converged());
    CHECK(solver.iterations() < solver.maxConvergenceIterations());
    solver.computeW();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() / 1000.;
    cout << elapsed << endl;

    // Converged at the requested mass flow: no reduction (a NaN reduces it), every
    // residual below its tolerance.
    const auto &dat = solver.data();
    REQUIRE(std::ranges::none_of(dat.W, [](T w) { return std::isnan(w); }));
    REQUIRE(solver.fullBladePassageMassFlow() == doctest::Approx(92.36282401553996));
    CHECK(solver.compressibilityResidualAverage().back() < solver.compressibilityRelTol());
    CHECK(solver.geomResidualMax().back() < solver.geomRelTol());
    CHECK(solver.periodicityUpStreamResidual().back() < solver.geomRelTol());
    CHECK(solver.periodicityDownStreamResidual().back() < solver.geomRelTol());
    // every downstream station is corrected (the last one was not: the residual
    // stalled at 2.9e-4)
    CHECK(solver.periodicityDownStreamResidual().back() < 1e-5);
    // Peak relative speed: 309.1 m/s with clang on arm64; 307.7 to 310.6 with the
    // geometric relaxation between 0.02 and 0.005
    auto W_max = *std::ranges::max_element(dat.W);
    CHECK(W_max == doctest::Approx(309.1).epsilon(0.01));
    CHECK(*std::ranges::max_element(dat.MW) < 1.);
    // isentropic flow with the meridional gas model: uniform entropy (it varied by
    // 0.02 J/(kg K) with gamma, R and cp set independently)
    auto [s_min, s_max] = std::ranges::minmax_element(dat.S);
    CHECK(*s_max - *s_min < 1e-8);
    // total enthalpy: cp Tt on the inlet plane, whatever the construction assumed
    for (size_t i = 0; i < ni; i++)
        CHECK(dat.H[i] == doctest::Approx(solver.gas().cp * dat.TT[i]));

    // Blade sides, leading to trailing edge included
    auto n_side = solver.trailingEdgeIndex() - solver.leadingEdgeIndex() + 1;
    auto ps1 = yams::getPsSide1(solver);
    auto ps2 = yams::getPsSide2(solver);
    REQUIRE(ps1.size() == n_side);
    REQUIRE(ps2.size() == n_side);
    REQUIRE(ps1.back() == dat.PS[solver.trailingEdgeIndex() * ni]);
    auto xyz = yams::getXYZSide1(solver);
    REQUIRE(xyz.size() == n_side);
    {
        // theta of the side, not the stream line parameter
        auto id = solver.leadingEdgeIndex() * ni;
        auto th = solver.mesh().TH[id];
        auto r = solver.mesh().R[id];
        CHECK(xyz.front()[0] == doctest::Approx(r * std::cos(th)));
        CHECK(xyz.front()[1] == doctest::Approx(r * std::sin(th)));
    }
    // the nodes stay on the stream line: m within its length
    auto [m1, m2] = gbs::abs_curv(solver.meridionalStreamLine()).bounds();
    CHECK(std::ranges::all_of(solver.mesh().M, [m1, m2](T m) { return m >= m1 && m <= m2; }));
    // a trailing edge plane before the leading edge one is refused, not wrapped around
    CHECK_THROWS_AS(yams::getValueStreamLine(dat.PS, 0, 3, 2, ni), std::invalid_argument);
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    #if YAMS_TESTS_PLOT
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    #endif
    // yams::plot(solver.mesh(), solver.data().Vm,"Vm");
    // yams::plot(solver.mesh(), solver.data().S,"S");
    // yams::plot(solver.mesh(), solver.data().PS,"PS");
    // yams::plot(solver.mesh(), solver.data().TS,"TS");


}
// The blade-to-blade solver takes the meridional solver's gas model
TEST_CASE("solver.gas_model")
{
    yams::PerfectGas<double> air;
    CHECK(air.gamma() == doctest::Approx(1004.5 / (1004.5 - 287.04)));
    CHECK(air.rho(288., 1.01325e5) == doctest::Approx(1.01325e5 / (287.04 * 288.)));
    CHECK(air.entropy(air.Tref, air.Pref) == 0.);
    CHECK(air.sound_speed(288.) == doctest::Approx(std::sqrt(air.gamma() * 287.04 * 288.)));

    yams::BladeToBladeCurvatureSolver<double> solver;
    solver.setGas({.R = 296.8, .cp = 1040.});
    CHECK(solver.gas().R == 296.8);
    CHECK(solver.gas().gamma() == doctest::Approx(1040. / (1040. - 296.8)));
    // gamma is not independent: it sets cp with the current R
    solver.setGamma(1.4);
    CHECK(solver.gas().cp == doctest::Approx(1.4 * 296.8 / 0.4));
    CHECK(solver.gas().gamma() == doctest::Approx(1.4));
    // not a gas: refused, the solver keeps its own
    CHECK_THROWS_AS(solver.setGamma(1.), std::invalid_argument);
    CHECK_THROWS_AS(solver.setCp(200.), std::invalid_argument);
    CHECK_THROWS_AS(solver.setGas({.R = 287.04, .cp = 287.04}), std::invalid_argument);
    CHECK_THROWS_AS(solver.setGasConstant(-1.), std::invalid_argument);
    CHECK(solver.gas().gamma() == doctest::Approx(1.4));
    // the density of the older gas model (euler2d's): p / (r T), it was its square root
    yams::GasPerfectGammaConstant<double> g{287.04, 1.4};
    CHECK(g.rho(288., 1.01325e5) == doctest::Approx(1.01325e5 / (287.04 * 288.)));
}
