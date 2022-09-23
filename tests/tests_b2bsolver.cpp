#include <gtest/gtest.h>
#include <chrono>
#include <iostream>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <gbs/bscbuild.h>
TEST(tridiagonal,thomas_vec)
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
            ASSERT_NEAR(d[i], r[i], 1e-9);
}

TEST(solver, base_novak)
{
    using namespace std;
    using T = double;

    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_coarse.vts";
    auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_medium.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_finest.vts";

    auto [pts, nj]  = yams::getVtkStructuredGridPoints<T>(f_name);
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
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    // yams::plot(solver.mesh(), solver.data().Vm,"Vm");
    // yams::plot(solver.mesh(), solver.data().S,"S");
    // yams::plot(solver.mesh(), solver.data().PS,"PS");
    // yams::plot(solver.mesh(), solver.data().TS,"TS");


}