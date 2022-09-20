#include <gtest/gtest.h>
#include <chrono>
#include <iostream>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <gbs/bscbuild.h>

TEST(solver, base_novak)
{
    using namespace std;
    using T = double;

    auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh.vts";

    auto [pts, nj]  = yams::getVtkStructuredGridPoints<T>(f_name);
    // auto stream_line = make_shared<gbs::BSCurve<T,2>>( 
    //     gbs::build_segment<T,2>({0.,1.},{1.,1.})
    // );
    auto stream_line = make_shared<gbs::BSCurve<T,2>>( 
        gbs::BSCurve<T,2>(
            {
                {0.0,1.0},
                {0.5,1.1},
                {1.0,1.1},
            },
            {0.,1.},
            {3,3},
            2
        )
    );
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

    yams::BladeToBladeCurvatureSolver<T> solver{ pts, nj, stream_line, 11 ,5,16};

    solver.setRotationSpeed( -150 );

    solver.setFullBladePassageMassFlow(
        // 74.
        // 74.8 
        // 74.9251
        // 75
        // 138
        // 70
        60
    );

    auto start = chrono::steady_clock::now();
    solver.computeW();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() / 1000.;
    cout << elapsed << endl;

    // yams::plot(solver.mesh(), nj ,solver.data().M,"Absolute Mach Number");
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    // yams::plot(solver.mesh(), solver.data().PT,"TotalPressure");
    // yams::plot(solver.mesh(), solver.data().TT,"TotalPressure");


}