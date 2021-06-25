#pragma once
#include <gtest/gtest.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridreader.h>

using namespace quiss;

const double PI = acos(-1.);

TEST(tests_diff, D1_O2)
{
    struct gp
    {
        double x;
        double y;
        double v;
    };

    Array2d<gp>   g(100,100);
    make_uniform_grid(2,1.,g);

    auto f     = [](auto & gp){gp.v = gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            auto di = D1_O2_i(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O2_j(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-3);
        }
    }

    make_uniform_clustered_grid(2,1.,g);

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            auto di = D1_O2_i(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O2_j(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-3);
        }
    }

    for (auto i = 2; i < ni; i++)
    {
        for (auto j = 2; j < nj; j++)
        {
            auto di = D1_O2_i_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O2_j_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-3);
        }
    }

    for (auto i = 0; i < ni-2; i++)
    {
        for (auto j = 0; j < nj-2; j++)
        {
            auto di = D1_O2_i_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O2_j_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-3);
        }
    }

    make_uniform_clustered_grid(2,1.,g,0.,PI / 8.);

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            auto di = D1_O2_i(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            // auto dj = D1_O2_j(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            // std:: cerr<< di << ' ' << dj << std::endl;
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            // ASSERT_NEAR(
            //     dj,
            //     dfqdy(g(i,j)),
            //     1e-3);
        }
    }
}

TEST(tests_diff, D1_O2_X)
{
    struct gp
    {
        double x;
        double y;
        double v;
    };

    ArrayX2d<gp>   g(100,100);
    make_uniform_grid(2,1.,g);

    auto f     = [](auto & gp){return  gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            g(i,j). v = f(g(i,j));
        }
    }

    for (auto i = 0; i < ni; i++)
    {
        for (auto j = 0; j < nj; j++)
        {
            auto di = D1_O2_i(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O2_j(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-3);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-3);
        }
    }
}

TEST(tests_diff, D1_O1)
{
    struct gp
    {
        double x;
        double y;
        double v;
    };

    Array2d<gp>   g(100,100);
    make_uniform_grid(2,1.,g);

    auto f     = [](auto & gp){gp.v = gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 1; i < ni; i++)
    {
        for (auto j = 1; j < nj; j++)
        {
            auto di = D1_O1_i_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O1_j_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-1);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-1);
        }
    }

    for (auto i = 0; i < ni-1; i++)
    {
        for (auto j = 0; j < nj-1; j++)
        {
            auto di = D1_O1_i_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O1_j_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-1);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-1);
        }
    }

    make_uniform_clustered_grid(2,1.,g);

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    for (auto i = 1; i < ni; i++)
    {
        for (auto j = 1; j < nj; j++)
        {
            auto di = D1_O1_i_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O1_j_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-1);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-1);
        }
    }

    for (auto i = 0; i < ni-1; i++)
    {
        for (auto j = 0; j < nj-1; j++)
        {
            auto di = D1_O1_i_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
            auto dj = D1_O1_j_fw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
            ASSERT_NEAR(
                di,
                dfqdx(g(i,j)),
                1e-1);
            ASSERT_NEAR(
                dj,
                dfqdy(g(i,j)),
                1e-1);
        }
    }

    // make_uniform_clustered_grid(2,1.,g,0.,PI / 8.);

    // std::for_each(
    //     g.begin(),
    //     g.end(),
    //     f
    // );

    // for (auto i = 1; i < ni; i++)
    // {
    //     for (auto j = 1; j < nj; j++)
    //     {
    //         auto di = D1_O1_j_bw(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.x;});
    //         // auto dj = D1_O2_j(g,i,j,[](const auto &g){return g.v;},[](const auto &g){return g.y;});
    //         // std:: cerr<< di << ' ' << dj << std::endl;
    //         ASSERT_NEAR(
    //             di,
    //             dfqdx(g(i,j)),
    //             1e-1);
    //         // ASSERT_NEAR(
    //         //     dj,
    //         //     dfqdy(g(i,j)),
    //         //     1e-3);
    //     }
    // }
}

TEST(tests_diff, D1_O2_test_001)
{
    struct gp
    {
        double x;
        double y;
        double v;
    };

    struct  gp_metrics
    {
        double x_ksi;
        double x_eth;
        double y_ksi;
        double y_eth;
        double J;
    };
    

    Array2d<gp>   g;
    quiss::read_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_001_250x21.vts");

    auto f     = [](auto & gp){gp.v = gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    Array2d<gp_metrics>   gp_metrics(ni,nj);
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);
    double err_max_x = 0.;
    double err_max_y = 0.;
    auto fx = [&g](const auto &gp) { return gp.x; };
    auto fy = [&g](const auto &gp) { return gp.y; };
    auto fv = [&g](const auto &gp) { return gp.v; };
    quiss::compute_metrics(g,fx,fy,gp_metrics);

    auto f_dv_dx = [&g, &gp_metrics, ksi, eth](size_t i, size_t j,const auto &fv_)
    {
        auto v_ksi = quiss::D1_O2_ksi(g, i, j, fv_, ksi);
        auto v_eth = quiss::D1_O2_eth(g, i, j, fv_, eth);
        auto gp_m = gp_metrics(i, j);
        return gp_m.J * gp_m.y_eth * v_ksi - gp_m.J * gp_m.y_ksi * v_eth;
    };
    auto f_dv_dy = [&g, &gp_metrics, ksi, eth](size_t i, size_t j,const auto &fv_)
    {
        auto v_ksi = quiss::D1_O2_ksi(g, i, j, fv_, ksi);
        auto v_eth = quiss::D1_O2_eth(g, i, j, fv_, eth);
        auto gp_m = gp_metrics(i, j);
        return gp_m.J * gp_m.x_eth * v_ksi + gp_m.J * gp_m.x_ksi * v_eth;
    };
    for (auto j = 0; j < nj; j++)
    {
        for (auto i = 0; i < ni; i++)
        {
            auto v_x = f_dv_dx(i, j, fv);
            auto v_y = f_dv_dy(i, j, fv);
            err_max_x = fmax(v_x - dfqdx(g(i, j)), err_max_x);
            err_max_y = fmax(v_y - dfqdy(g(i, j)), err_max_y);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                2e-3);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                2e-3);
        }
    }

    std::cout << "err_max: " << err_max_x << " " << err_max_y << std::endl;
}