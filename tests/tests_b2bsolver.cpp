#include "yams_doctest.h"
#include <chrono>
#include <iostream>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <baldeToBlade/bladeToBladeCurvatureSolverPost.h>
#include <algorithm>
#include <cmath>
#include <gbs/bscbuild.h>
#include <string>
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
    // Peak relative speed: 308 m/s with the geometric tolerance at 1e-4 (relaxation
    // 0.02 or 0.01); about 300 at the default 1e-3.
    auto W_max = *std::ranges::max_element(dat.W);
    CHECK(W_max == doctest::Approx(308.).epsilon(0.05));
    CHECK(*std::ranges::max_element(dat.MW) < 1.);

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
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    #if YAMS_TESTS_PLOT
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    #endif
    // yams::plot(solver.mesh(), solver.data().Vm,"Vm");
    // yams::plot(solver.mesh(), solver.data().S,"S");
    // yams::plot(solver.mesh(), solver.data().PS,"PS");
    // yams::plot(solver.mesh(), solver.data().TS,"TS");


}