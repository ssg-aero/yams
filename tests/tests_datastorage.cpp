#include <gtest/gtest.h>
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <execution>
#include <algorithm>
#include <chrono>

using namespace quiss;

const double PI = acos(-1.);

TEST(tests_datastorage, Array2d)
{
    Array2d<double> a(3,5);
    ASSERT_DOUBLE_EQ(a(0,0),0.);
    std::for_each(
        a.begin(0),
        a.end(0),
        [](const auto &v_){ ASSERT_DOUBLE_EQ(v_,0.); }
    );
    
    Array2d<double> b(3,5,1.);
    ASSERT_DOUBLE_EQ(b(0,0),1.);
    std::for_each(
        b.begin(1),
        b.end(1),
        [](const auto &v_){ ASSERT_DOUBLE_EQ(v_,1.); }
    );
}

TEST(tests_datastorage, ArrayX2d)
{
    ArrayX2d<double> a(3,5);// on contrary of std::vector, values are not initialized to 0.
    ASSERT_NEAR(a(0,0),0.,1e-30);
    std::for_each(
        a.begin(0),
        a.end(0),
        [](const auto &v_){ ASSERT_NEAR(v_,0.,1e-30); }
    );
    
    ArrayX2d<double> b(3,5,1.);
    ASSERT_DOUBLE_EQ(b(0,0),1.);
    std::for_each(
        b.begin(1),
        b.end(1),
        [](const auto &v_){ ASSERT_DOUBLE_EQ(v_,1.); }
    );

    std::for_each(
        std::execution::par,
        b.begin(),
        b.end(),
        [](const auto &v_){ ASSERT_DOUBLE_EQ(v_,1.); }
    );
}

template <typename T>
inline auto distance(const GridPoint<T> &gp1,const GridPoint<T> &gp2) -> T
{
    return sqrt((gp1.y - gp2.y) * (gp1.y - gp2.y) + (gp1.x - gp2.x) * (gp1.x - gp2.x));
}

template <typename T>
inline auto compute_abscissas(Grid<T> &g)
{
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            g(i, j).m = i == 0 ? 0. : g(i - 1, j).m + distance(g(i, j), g(i - 1, j));
            g(i, j).l = j == 0 ? 0. : g(i, j - 1).l + distance(g(i, j), g(i, j - 1));
        }
    }
}

template <typename T>
auto m = [](GridPoint<T> &gp){return gp.m;};
template <typename T>
auto l = [](GridPoint<T> &gp){return gp.l;};
template <typename T>
auto z = [](GridPoint<T> &gp){return gp.x;};
template <typename T>
auto r = [](GridPoint<T> &gp){return gp.y;};
TEST(tests_datastorage, Grid)
{
    size_t ni = 210;
    size_t nj = 130;
    Grid<double> g(ni,nj);
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
            ASSERT_NEAR(l<double>(g(i,nj-1)),r1-r2,1e-7);
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
            ASSERT_NEAR(g(i, nj - 1).l, r1 - r2, 1e-7);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cerr  << "direct access "
                  << " took " << ms.count() << " ms\n";
    }

    ASSERT_DOUBLE_EQ(g(0,0).x,0.);
    ASSERT_DOUBLE_EQ(g(0,0).y,r3);
    ASSERT_DOUBLE_EQ(g(0,nj-1).x,r1-r2);
    ASSERT_DOUBLE_EQ(g(0,nj-1).y,r3);
    ASSERT_DOUBLE_EQ(g(ni-1,0).x,2*z3);
    ASSERT_DOUBLE_EQ(g(ni-1,0).y,r3);
    std::cerr << g(ni/2,0).y << ' ' << g(ni/2,0).x << std::endl;
    GridPoint<double> O;
    O.x = z3;
    O.y = r3;
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            ASSERT_NEAR(distance(g(i, j) , O),r,1e-7);
        }
    }

}

TEST(tests_datastorage, Grid_Diff)
{
    size_t ni = 100;
    size_t nj = 30;
    Grid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    
    make_uniform_grid(2.,1.,g);
    for (auto j = 0; j < nj; j++)
    {
        ASSERT_DOUBLE_EQ(g(ni - 1, j).x, 2.);
    }
        for (auto i = 0; i < ni; i++)
    {
        ASSERT_DOUBLE_EQ(g(i, nj - 1).y, 1.);
    }
    compute_abscissas(g);
    for (auto j = 0; j < nj; j++)
    {
        ASSERT_DOUBLE_EQ(g(ni - 1, j).m, 2.);
    }
        for (auto i = 0; i < ni; i++)
    {
        ASSERT_DOUBLE_EQ(g(i, nj - 1).l, 1.);
    }

    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            double drqdm,dzqdm,drqdl,dzqdl;
            drqdm = D1_O2_i(g,i,j,r<double>,m<double>);
            dzqdm = D1_O2_i(g,i,j,z<double>,m<double>);
            ASSERT_DOUBLE_EQ(drqdm,0.);
            ASSERT_DOUBLE_EQ(dzqdm,1.);
            drqdl = D1_O2_j(g,i,j,r<double>,l<double>);
            dzqdl = D1_O2_j(g,i,j,z<double>,l<double>);
            ASSERT_DOUBLE_EQ(drqdl,1.);
            ASSERT_DOUBLE_EQ(dzqdl,0.);
        }
    }

    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    double drqdm, dzqdm, phi,drqdl,dzqdl,gam;
    for (auto i = 0; i < ni; i++)
    {
        auto phi_ = -PI /2  + PI * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto gam_ =  t1  + (t2-t1) * j / (nj - 1.);
            drqdm = D1_O2_i(g,i,j,r<double>,m<double>);
            dzqdm = D1_O2_i(g,i,j,z<double>,m<double>);
            phi = atan2( drqdm,dzqdm ); // Stream line angle
            ASSERT_NEAR(phi,phi_,1e-5);
            drqdl = D1_O2_j(g,i,j,r<double>,l<double>);
            dzqdl = D1_O2_j(g,i,j,z<double>,l<double>);
            gam = atan2( dzqdl,drqdl ); // Span line angle
            // std::cerr << gam << ' ' << dzqdl << ' ' << drqdl << std::endl;
            // ASSERT_NEAR(gam,gam_,1e-5);
        }
    }
}

TEST(tests_datastorage, GridStd_perfo)
{
    size_t ni = 21000;
    size_t nj = 13000;
    
    Array2d<quiss::GridPoint<double>> g{ni,nj};
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

TEST(tests_datastorage, GridX_perfo)
{
    size_t ni = 21000;
    size_t nj = 13000;
    
    ArrayX2d<quiss::GridPoint<double>> g{ni,nj};
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

TEST(tests_datastorage, convertion)
{
    GridPoint<double> gpd;
    GridPoint<float> gpf;
    gpf.x=1.2;
    gpf.y=3.1;

    size_t n =  sizeof(gpf) / sizeof(float);
    quiss::copy(gpf,gpd,n);
    ASSERT_NEAR(gpd.x, gpd.x, 1e-30);
    ASSERT_NEAR(gpd.y, gpd.y, 1e-30);

    size_t ni = 21;
    size_t nj = 13;
    
    ArrayX2d<quiss::GridPoint<double>> g{ni,nj};
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

    Array2d<quiss::GridPoint<float>> g1{ni,nj};

    quiss::copy(g,g1);

    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            ASSERT_NEAR(g1(i, j).y , r3 - r * sin(th), 1e-6);
            ASSERT_NEAR(g1(i, j).x , z3 - r * cos(th), 1e-6);
        }
    }

}