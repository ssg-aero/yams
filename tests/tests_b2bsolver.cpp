#include <gtest/gtest.h>
#include <chrono>
#include <iostream>

#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <baldeToBlade/gridReader.h>
#include <gbs/bscbuild.h>
#include <matrix/tridiagonal.h>

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

    auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_coarse.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_medium.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh.vts";
    // auto f_name = "C:/Users/sebastien/workspace/yams/tests/b2bmsh_finest.vts";

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

    yams::BladeToBladeCurvatureSolver<T> solver{ 
        pts, 
        nj, 
        stream_line, 
        11,
        // 4,
        // 15
        3,11
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

    // solver.applyDeltaStagnationLineDownStream({ 0.000, 0.00, 0.00});
    // solver.applyDeltaStagnationLineUpStream({ 0.006,-0.002, 0.00});


    auto start = chrono::steady_clock::now();
    solver.computeW();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() / 1000.;
    cout << elapsed << endl;
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    // yams::plot(solver.mesh(), solver.data().Vm,"Vm");
    // yams::plot(solver.mesh(), solver.data().S,"S");
    // yams::plot(solver.mesh(), solver.data().PS,"PS");
    // yams::plot(solver.mesh(), solver.data().TS,"TS");
    // return;

    const auto &TH = solver.mesh().TH;
    const auto &W = solver.data().W;
    auto jLe = solver.leadingEdgeIndex();
    auto jTe = solver.trailingEdgeIndex();
    auto j1 = jLe - 1;
    auto j2 = jLe - 2;
    // // auto j1 = jTe + 1;
    // // auto j2 = jTe + 2;
    // vector<size_t> j_stations{jLe - 1,jLe - 2};

    int nit = 2;
    int nit_sub = 3;
    for(int i{}; i < nit * nit_sub; i++)
    {
        size_t i_mid = (ni - 1) / 2 + 1;
        auto f_A = [i_mid, ni, &VM = solver.data().Vm, &R = solver.mesh().R](size_t j)
        {
            auto rVm = R[i_mid + ni * j]*VM[i_mid + ni * j];
            return 2 * rVm * rVm;
        };
        auto f_B = [i_mid, ni, &VM = solver.data().Vm, &DG3 = solver.data().G3](size_t j)
        {
            return 2 * VM[i_mid + ni * j] * DG3[i_mid + ni * j];
        };

        auto f_a = [i_mid, ni, &m = solver.mesh().M](size_t j, vector<T> &A, vector<T> &B)
        {
            auto id1 = i_mid + ni * (j - 1);
            auto id2 = i_mid + ni * (j);
            auto id3 = i_mid + ni * (j + 1);
            return 2 * A[j] / (m[id3] - m[id1]) / (m[id2] - m[id1]) + B[j] / 2 / (m[id2] - m[id1]);
        };

        auto f_c = [i_mid, ni, &m = solver.mesh().M](size_t j, vector<T> &A, vector<T> &B)
        {
            auto id1 = i_mid + ni * (j - 1);
            auto id2 = i_mid + ni * (j);
            auto id3 = i_mid + ni * (j + 1);
            return 2 * A[j] / (m[id3] - m[id1]) / (m[id3] - m[id2]) + B[j] / 2 / (m[id3] - m[id2]);
        };

        auto f_b = [i_mid, ni, &m = solver.mesh().M](size_t j, vector<T> &A, vector<T> &B)
        {
            auto id1 = i_mid + ni * (j - 1);
            auto id2 = i_mid + ni * (j);
            auto id3 = i_mid + ni * (j + 1);
            return -2 * A[j] / (m[id3] - m[id1]) * (1 / (m[id3] - m[id2]) + 1 / (m[id2] - m[id1])) - B[j] / 2 * (1 / (m[id3] - m[id2]) + 1 / (m[id2] - m[id1]));
        };

        auto f_d = [ ni,&W=solver.data().W, &TH=solver.mesh().TH](size_t j)
        {
            auto id1 = ni * j;
            auto id2 = ni - 1 + ni * j; 
            return (W[id2] * W[id2] - W[id1] * W[id1]) / (TH[id2] - TH[id1]);
        };

        vector<T> A(nj), B(nj);
        A[j1] = f_A(j1);
        A[j2] = f_A(j2);
        B[j1] = f_B(j1);
        B[j2] = f_B(j2);

        auto b1 = f_b(j1, A, B);
        auto c1 = f_c(j1, A, B);
        auto a2 = f_a(j2, A, B);
        auto b2 = f_b(j2, A, B);
        auto c2 = f_c(j2, A, B);


        auto d1 = f_d(j1);
        auto d2 = f_d(j2);

        vector<T> a{0.,a2};
        vector<T> b{b1,b2+c2};
        vector<T> c{c1,0.};
        vector<T> d{-d1/nit_sub,-d2/nit_sub};
        vector<T> delta(2);
        thomas_algorithm(a, b, c, d);

        

        // auto DM = b1 * (b2 + c2) - a2 * c1;
        // auto delta1 = (d2 * (b2 + c2) - d2 * c1) / DM;
        // auto delta2 = (b1 * d2 - a2 * d1) / DM;

        // solver.applyDeltaStagnationLineDownStream({-delta1 / nit,-delta2 / nit});
        // solver.applyDeltaStagnationLineUpStream({-delta1/nit,-delta2/nit});
        if(j2>j1)
            solver.applyDeltaStagnationLineDownStream(d);
        else if(j1>j2)
            solver.applyDeltaStagnationLineUpStream(d);
        else
            break;
        solver.computeMeshData();
        solver.computeW();
        for( int j = j1 ; j >= 0; j--)
        // for( int j = j1 ; j < nj; j++)
            cout << W[ni - 1 + ni * j] - W[0 + ni * j] << " ";
        cout << endl;

    }

    // for(int i=0; i < 10; i++)
    // {
    //     solver.applyDeltaStagnationLineUpStream({-0.002});
    //     solver.computeMeshData();
    //     solver.computeW();
    //     for( int j = j1 ; j >= 0; j--)
    //         cout << W[ni - 1 + ni * j] - W[0 + ni * j] << " ";
    //     cout << endl;
    // }

    // yams::plot(solver.mesh(), nj ,solver.data().M,"Absolute Mach Number");
    // yams::plot(solver.mesh(), solver.data().W,"Relative Speed");
    yams::plot(solver.mesh(), solver.data().Q,"Relative Speed");
    yams::plot(solver.mesh(), solver.data().MW,"Relative Mach Number");
    // yams::plot(solver.mesh(), solver.data().PT,"TotalPressure");
    // yams::plot(solver.mesh(), solver.data().TT,"TotalPressure");

    // cout << solver.data().W[0 + ni * (2)] << " ;" << solver.data().W[ni-1 + ni * (2)] << endl;


}