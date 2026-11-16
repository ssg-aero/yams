#include "yams_doctest.h"
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>
#include <gbs/execution.h>
#include <algorithm>
#include <chrono>

using namespace yams;

const double PI = acos(-1.);

TEST_CASE("tests_datastorage.Array2d")
{
    Array2d<double> a(3,5);
    REQUIRE(a(0,0) == yams_tests::approx_ulp(0.));
    std::for_each(
        a.begin(0),
        a.end(0),
        [](const auto &v_){ REQUIRE(v_ == yams_tests::approx_ulp(0.)); }
    );
    
    Array2d<double> b(3,5,1.);
    REQUIRE(b(0,0) == yams_tests::approx_ulp(1.));
    std::for_each(
        b.begin(1),
        b.end(1),
        [](const auto &v_){ REQUIRE(v_ == yams_tests::approx_ulp(1.)); }
    );
}

// TEST_CASE("tests_datastorage.ArrayX2d")
// {
//     ArrayX2d<double> a(3,5);// on contrary of std::vector, values are not initialized to 0.
//     REQUIRE(std::abs(a(0,0) - 0.) <= 1e-30);
//     std::for_each(
//         a.begin(0),
//         a.end(0),
//         [](const auto &v_){ REQUIRE(std::abs(v_ - 0.) <= 1e-30); }
//     );
    
//     ArrayX2d<double> b(3,5,1.);
//     REQUIRE(b(0,0) == yams_tests::approx_ulp(1.));
//     std::for_each(
//         b.begin(1),
//         b.end(1),
//         [](const auto &v_){ REQUIRE(v_ == yams_tests::approx_ulp(1.)); }
//     );

//     std::for_each(
//         GBS_PAR_EXEC
//         b.begin(),
//         b.end(),
//         [](const auto &v_){ REQUIRE(v_ == yams_tests::approx_ulp(1.)); }
//     );
// }


template <typename T>
auto m = [](MeridionalGridPoint<T> &gp){return gp.m;};
template <typename T>
auto l = [](MeridionalGridPoint<T> &gp){return gp.l;};
template <typename T>
auto z = [](MeridionalGridPoint<T> &gp){return gp.x;};
template <typename T>
auto r = [](MeridionalGridPoint<T> &gp){return gp.y;};
template <typename T>
auto phi = [](MeridionalGridPoint<T> &gp){return gp.phi;};
TEST_CASE("tests_datastorage.Grid")
{
    size_t ni = 210;
    size_t nj = 130;
    MeridionalGrid<double> g(ni,nj);
    auto r1 =  2.;
    auto r2 =  1.;
    auto r3 =  3.;
    auto z3 =  r1;
    auto t1 =  PI;
    auto t2 =  2*PI;
    make_circular_grid(r1,r2,t1,t2,{z3,r3},g);
    compute_abscissas(g);

    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for(int c = 0 ; c < 1000 ; c++)
        for (auto i = 0; i < ni; i++)
        {
            REQUIRE(std::abs((l<double>(g(i,nj-1))) - (r1-r2)) <= 1e-7);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr  << std::fixed << "access with lambda "
                  << " took " << ms.count() << " ms\n";
    }

    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for(int c = 0 ; c < 1000 ; c++)
        for (auto i = 0; i < ni; i++)
        {
            REQUIRE(std::abs(g(i, nj - 1).l - (r1 - r2)) <= 1e-7);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr  << "direct access "
                  << " took " << ms.count() << " ms\n";
    }

    REQUIRE(g(0,0).x == yams_tests::approx_ulp(0.));
    REQUIRE(g(0,0).y == yams_tests::approx_ulp(r3));
    REQUIRE(g(0,nj-1).x == yams_tests::approx_ulp(r1-r2));
    REQUIRE(g(0,nj-1).y == yams_tests::approx_ulp(r3));
    REQUIRE(g(ni-1,0).x == yams_tests::approx_ulp(2*z3));
    REQUIRE(g(ni-1,0).y == yams_tests::approx_ulp(r3));
    std::cerr << g(ni/2,0).y << ' ' << g(ni/2,0).x << std::endl;
    MeridionalGridPoint<double> O;
    O.x = z3;
    O.y = r3;
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            REQUIRE(std::abs(distance(g(i, j) , O) - r) <= 1e-7);
        }
    }

}

TEST_CASE("tests_datastorage.Grid_Diff")
{
    size_t ni = 100;
    size_t nj = 30;
    MeridionalGrid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    
    make_uniform_grid(2.,1.,g);
    for (auto j = 0; j < nj; j++)
    {
        REQUIRE(g(ni - 1, j).x == yams_tests::approx_ulp(2.));
    }
        for (auto i = 0; i < ni; i++)
    {
        REQUIRE(g(i, nj - 1).y == yams_tests::approx_ulp(1.));
    }
    compute_abscissas(g);
    for (auto j = 0; j < nj; j++)
    {
        REQUIRE(g(ni - 1, j).m == yams_tests::approx_ulp(2.));
    }
        for (auto i = 0; i < ni; i++)
    {
        REQUIRE(g(i, nj - 1).l == yams_tests::approx_ulp(1.));
    }

    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            double drqdm,dzqdm,drqdl,dzqdl;
            drqdm = D1_O2_i(g,i,j,r<double>,m<double>);
            dzqdm = D1_O2_i(g,i,j,z<double>,m<double>);
            REQUIRE(drqdm == yams_tests::approx_ulp(0.));
            REQUIRE(dzqdm == yams_tests::approx_ulp(1.));
            drqdl = D1_O2_j(g,i,j,r<double>,l<double>);
            dzqdl = D1_O2_j(g,i,j,z<double>,l<double>);
            REQUIRE(drqdl == yams_tests::approx_ulp(1.));
            REQUIRE(dzqdl == yams_tests::approx_ulp(0.));
        }
    }

    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    for (auto j = 0; j < nj; j++)
    {
        // REQUIRE(g(ni - 1, j).m == yams_tests::approx_ulp(2.));
    }
        for (auto i = 0; i < ni; i++)
    {
        // REQUIRE(g(i, nj - 1).l == yams_tests::approx_ulp(1.));
    }
    compute_angles(g);
    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            REQUIRE(fmod(-PI /2  + PI * i / (ni - 1.)-g(i,j).phi,2*PI) <= 1e-5);
            // gam = atan2( dzqdl,drqdl ); // Span line angle
            // std::cerr << gam << ' ' << dzqdl << ' ' << drqdl << std::endl;
            // REQUIRE(std::abs(gam - gam_) <= 1e-5);
        }
    }
    compute_curvature(g);
    double DphiDm;
    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            // DphiDm = D1_O2_i(g, i, j, phi<double>, m<double>);
            DphiDm = g(i,j).cur;
            REQUIRE(std::abs((1. / DphiDm) - (r1 + (r2 - r1) * j / (nj - 1.))) <= 1e-3);
        }
    }
}

TEST_CASE("tests_datastorage.GridStd_perfo")
{
    size_t ni = 2100;
    size_t nj = 1300;
    
    Array2d<yams::MeridionalGridPoint<double>> g{ni,nj};
    auto r1 =  1.;
    auto r2 =  2.;
    auto r3 =  3.;
    auto z3 =  3.;
    auto t1 = -PI/2.;
    auto t2 =  PI/2.;
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            g(i, j).y = r3 - r * sin(th);
            g(i, j).x = z3 - r * cos(th);
        }
    }
}

#ifdef USE_XT
TEST_CASE("tests_datastorage.GridX_perfo")
{
    size_t ni = 2100;
    size_t nj = 1300;
    
    ArrayX2d<yams::MeridionalGridPoint<double>> g{ni,nj};
    auto r1 =  1.;
    auto r2 =  2.;
    auto r3 =  3.;
    auto z3 =  3.;
    auto t1 = -PI/2.;
    auto t2 =  PI/2.;
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            g(i, j).y = r3 - r * sin(th);
            g(i, j).x = z3 - r * cos(th);
        }
    }
}
TEST_CASE("tests_datastorage.convertion")
{
    MeridionalGridPoint<double> gpd;
    MeridionalGridPoint<float> gpf;
    gpf.x=1.2;
    gpf.y=3.1;

    size_t n =  sizeof(gpf) / sizeof(float);
    yams::copy(gpf,gpd,n);
    REQUIRE(std::abs(gpd.x - gpd.x) <= 1e-30);
    REQUIRE(std::abs(gpd.y - gpd.y) <= 1e-30);

    size_t ni = 21;
    size_t nj = 13;
    
    ArrayX2d<yams::MeridionalGridPoint<double>> g{ni,nj};
    auto r1 =  1.;
    auto r2 =  2.;
    auto r3 =  3.;
    auto z3 =  3.;
    auto t1 = -PI/2.;
    auto t2 =  PI/2.;
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            g(i, j).y = r3 - r * sin(th);
            g(i, j).x = z3 - r * cos(th);
        }
    }

    Array2d<yams::MeridionalGridPoint<float>> g1{ni,nj};

    yams::copy(g,g1);

    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            REQUIRE(std::abs(g1(i, j).y - (r3 - r * sin(th))) <= 1e-6);
            REQUIRE(std::abs(g1(i, j).x - (z3 - r * cos(th))) <= 1e-6);
        }
    }

}

#include <cstddef>
#include <vector>
#include "xsimd/xsimd.hpp"
#include "xsimd/stl/algorithms.hpp"

const size_t size_v = 1000000;
const size_t n_it = 1000;

TEST_CASE("tests_datastorage.xsimd_perfo")
{
    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        using vector_type = std::vector<double, xsimd::aligned_allocator<double, XSIMD_DEFAULT_ALIGNMENT>>;
        vector_type a(size_v);
        vector_type b(size_v);
        vector_type c(size_v);
        for (auto i = 0; i < n_it; i++)
            xsimd::transform(
                // GBS_PAR_EXEC // not working
                a.begin(), a.end(), b.begin(), c.begin(),
                [](const auto &x, const auto &y) { return (x + y) / 2.; });
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr << std::fixed << "xtensor container "
                  << " took " << ms.count() << " ms\n";
    }
    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::vector<double> a(size_v);
        std::vector<double> b(size_v);
        std::vector<double> c(size_v);
        for (auto i = 0; i < n_it; i++)
            std::transform(
                a.begin(), a.end(), b.begin(), c.begin(),
                [](const auto &x, const auto &y) { return (x + y) / 2.; });
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr << std::fixed << "stl container "
                  << " took " << ms.count() << " ms\n";
    }
        {
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::vector<double> a(size_v);
        std::vector<double> b(size_v);
        std::vector<double> c(size_v);
        for (auto i = 0; i < n_it; i++)
            std::transform(
                GBS_PAR_EXEC
                a.begin(), a.end(), b.begin(), c.begin(),
                [](const auto &x, const auto &y) { return (x + y) / 2.; });
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr << std::fixed << "stl container  execution::par"
                  << " took " << ms.count() << " ms\n";
    }
}
#endif