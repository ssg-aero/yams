#include "gridsbuilders.h"

#include <yams/diffop.h>
#include <yams/gridreader.h>
#include <yams/gridmetrics.h>
#include <yams/gridrender.h>
#include <yams/meshtools.h>

#include <gtest/gtest.h>

#include <gbs-render/vtkgridrender.h>

using namespace yams;

const double PI = acos(-1.);
const bool PLOT_ON = true;
#ifdef _WIN32
    const std::string test_files_path = "C:/Users/sebastien/workspace/tbslib/tests/";
#else
    const std::string test_files_path = "../../../tbslib/tests/";
#endif
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



    Array2d<gp>   g;
    vtkNew<vtkXMLStructuredGridReader> reader;
    reader->SetFileName((test_files_path+"/out/test_001_250x21.vts").c_str());
    // reader->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_250x21.vts");
    reader->Update();
    auto sgrid = reader->GetOutput();
    auto points= sgrid->GetPoints();
    auto dims  =sgrid->GetDimensions();
    size_t ni = dims[0];
    size_t nj = dims[1];
    ASSERT_TRUE(nj != 0);
    ASSERT_TRUE(ni != 0);
    g.resize(ni,nj);
    vtkIdType id {};
    for(size_t j {} ; j < nj ; j++)
    {
        for(size_t i {} ; i < ni ; i++)
        {
            auto pt = points->GetPoint(id);
            g(i,j).x = pt[0];
            g(i,j).y = pt[1];
            id++;
        }
    }
    // yams::Array2d<double>(g,"../../../tbslib/tests/out/test_001_250x21.vts");
        // yams::read_vtk_grid(g,"../../../tbslib/tests/out/test_002.vts");

    auto f     = [](auto & gp){gp.v = gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    ni = g.nRows();
    nj = g.nCols();
    Array2d<yams::Grid2dMetricsPoint<double>>   gp_metrics(ni,nj);
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);
    double err_max_x = 0.;
    double err_max_y = 0.;
    auto fx = [&g](const auto &gp) { return gp.x; };
    auto fy = [&g](const auto &gp) { return gp.y; };
    auto fv = [&g](const auto &gp) { return gp.v; };
    yams::compute_metrics(g,fx,fy,gp_metrics);

    for (auto j = 0; j < nj; j++)
    {
        for (auto i = 0; i < ni; i++)
        {
            auto v_x = yams::D1_O2_dx1(g,gp_metrics,i,j,ksi,eth,fv);
            auto v_y = yams::D1_O2_dx2(g,gp_metrics,i,j,ksi,eth,fv);
            err_max_x = fmax(v_x - dfqdx(g(i, j)), err_max_x);
            err_max_y = fmax(v_y - dfqdy(g(i, j)), err_max_y);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                1e-2);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                1e-2);
        }
    }

    std::cout << "err_max: " << err_max_x << " " << err_max_y << std::endl;
}
/*
TEST(tests_diff, D1_O2_test_001_ghost)
{
    struct gp
    {
        double x;
        double y;
        double v;
    };



    Array2d<gp,1,1>   g;
    vtkNew<vtkXMLStructuredGridReader> reader;
    reader->SetFileName((test_files_path+"/out/test_001_250x21.vts").c_str());
    // reader->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_250x21.vts");
    reader->Update();
    auto sgrid = reader->GetOutput();
    auto points= sgrid->GetPoints();
    auto dims  =sgrid->GetDimensions();
    size_t ni = dims[0];
    size_t nj = dims[1];
    ASSERT_TRUE(nj != 0);
    ASSERT_TRUE(ni != 0);
    g.resize(ni-2,nj-2);
    vtkIdType id {};
    for(std::intmax_t j {-1} ; j <= nj ; j++)
    {
        for(std::intmax_t i {-1} ; i <= ni ; i++)
        {
            auto pt = points->GetPoint(id);
            g(i,j).x = pt[0];
            g(i,j).y = pt[1];
            id++;
        }
    }
    // yams::Array2d<double>(g,"../../../tbslib/tests/out/test_001_250x21.vts");
        // yams::read_vtk_grid(g,"../../../tbslib/tests/out/test_002.vts");

    auto f     = [](auto & gp){gp.v = gp.x * gp.y + gp.y * sin(gp.x);};
    auto dfqdx = [](auto & gp){return gp.y + gp.y * cos(gp.x);};
    auto dfqdy = [](auto & gp){return gp.x + sin(gp.x);};

    std::for_each(
        g.begin(),
        g.end(),
        f
    );

    ni = g.nRows();
    nj = g.nCols();
    Array2d<yams::Grid2dMetricsPoint<double>>   gp_metrics(ni,nj);
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);
    double err_max_x = 0.;
    double err_max_y = 0.;
    auto fx = [&g](const auto &gp) { return gp.x; };
    auto fy = [&g](const auto &gp) { return gp.y; };
    auto fv = [&g](const auto &gp) { return gp.v; };
    yams::compute_metrics(g,fx,fy,gp_metrics);

    for (auto j = 0; j < nj; j++)
    {
        for (auto i = 0; i < ni; i++)
        {
            auto v_x = yams::D1_O2_dx1(g,gp_metrics,i,j,ksi,eth,fv);
            auto v_y = yams::D1_O2_dx2(g,gp_metrics,i,j,ksi,eth,fv);
            err_max_x = fmax(v_x - dfqdx(g(i, j)), err_max_x);
            err_max_y = fmax(v_y - dfqdy(g(i, j)), err_max_y);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                1e-2);
            ASSERT_NEAR(
                v_x,
                dfqdx(g(i, j)),
                1e-2);
        }
    }

    std::cout << "err_max: " << err_max_x << " " << err_max_y << std::endl;
}
*/

TEST(tests_diff, D1_O2_cir)
{
    using T = double;
    size_t ni = 5000;
    size_t nj = 1500;
    T d_ksi = 1. / (ni - 1.);
    T d_eth = 1. / (nj - 1.);
    MeridionalGrid<T> g(ni,nj);
    Array2d<Grid2dMetricsPoint<T>> g_metrics(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto t1 = std::numbers::pi_v<T>;
    auto t2 = std::numbers::pi_v<T> * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    compute_metrics(g,[](const auto & gp){return gp.m;},[](const auto & gp){return gp.l;},g_metrics);
    compute_angles(g);

    auto  f_Vm    = [](const auto &gp){return sin(gp.m)*gp.l + cos(gp.l)*gp.m;};
    auto df_Vm_dm = [](const auto &gp){return  cos(gp.m)*gp.l + cos(gp.l);};
    auto df_Vm_dl = [](const auto &gp){return -sin(gp.l)*gp.m + sin(gp.m);};
    std::for_each(
        g.begin(), g.end(),
        [&f_Vm](auto &gp){gp.Vm = f_Vm(gp);}
    );

    for(size_t i{}; i < ni; i++)
    {
        for(size_t j{}; j < nj; j++)
        {
            auto df_Vm_dm_ = D1_O2_dx1(g, g_metrics, i, j, d_ksi, d_eth, [](const auto &gp)
                                       { return gp.Vm; });
            auto df_Vm_dl_ = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, [](const auto &gp)
                                      { return gp.Vm; }); 
            if (PLOT_ON)
            {
                g(i, j).Vu = df_Vm_dm_; // use variable only for storage
                g(i, j).H  = df_Vm_dl_;
            }
            ASSERT_NEAR(df_Vm_dm_, df_Vm_dm(g(i, j)), 1e-5);
            ASSERT_NEAR(df_Vm_dl_, df_Vm_dl(g(i, j)), 1e-5);
        }
    }
    if (PLOT_ON)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "m", [](const auto &gp)
                  { return gp.m; });
        add_value(g, structuredGrid, "l", [](const auto &gp)
                  { return gp.l; });
        add_value(g, structuredGrid, "Vm", [](const auto &gp)
                  { return gp.Vm; });
        add_value(g, structuredGrid, "dVm_dm", [&df_Vm_dm](const auto &gp)
                  { return df_Vm_dm(gp); });
        add_value(g, structuredGrid, "dVm_dl", [&df_Vm_dl](const auto &gp)
                  { return df_Vm_dl(gp); });
        add_value(g, structuredGrid, "dVm_dm_err", [&](const auto &gp)
                  { return df_Vm_dm(gp) - gp.Vu; });
        add_value(g, structuredGrid, "dVm_dl_err", [&](const auto &gp)
                  { return df_Vm_dl(gp) - gp.H; });
        plot_vtkStructuredGrid(structuredGrid, "Vm", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dm", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dl", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dm_err", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dl_err", false);
    }
}

TEST(tests_diff, D1_O2_bump)
{
    using T = double;
    std::vector<std::array<T,2>> poles1{
            {0.0,0.0},
            {1.0*std::numbers::pi_v<T>,0.0},
        };
    std::vector<std::array<T,2>> poles3{
            {0.0,0.5*std::numbers::pi_v<T>},
            {1.0*std::numbers::pi_v<T>,0.5*std::numbers::pi_v<T>},
        };
    std::vector<std::array<T,2>> poles2{
            {0.0,1.0*std::numbers::pi_v<T>},
            {0.1*std::numbers::pi_v<T>,1.0*std::numbers::pi_v<T>},
            {0.4*std::numbers::pi_v<T>,0.7*std::numbers::pi_v<T>},
            {0.6*std::numbers::pi_v<T>,0.7*std::numbers::pi_v<T>},
            {0.9*std::numbers::pi_v<T>,1.0*std::numbers::pi_v<T>},
            {1.0*std::numbers::pi_v<T>,1.0*std::numbers::pi_v<T>},
        };
    std::vector<T> knots1{0.,1.};
    std::vector<size_t> mult1{2,2};
    size_t deg1{1};
    std::vector<T> knots2{0.,0.5,1.};
    std::vector<size_t> mult2{5,1,5};
    size_t deg2{4};

    auto crv1 = std::make_shared<gbs::BSCurve<T,2>>(
        poles1,
        knots1,
        mult1,
        deg1
    );
    auto crv3 = std::make_shared<gbs::BSCurve<T,2>>(
        poles3,
        knots1,
        mult1,
        deg1
    );
    auto crv2 = std::make_shared<gbs::BSCurve<T,2>>(
        poles2,
        knots2,
        mult2,
        deg2
    );
    size_t nu{300};
    size_t nv{150};
    yams::crv_vector<T> crv_lst{crv1, crv3, crv2};
    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = mesh_channel<T>(crv_lst, knots1, nv, nu);
    T d_ksi = 1. / (ni - 1.);
    T d_eth = 1. / (nj - 1.);
    auto sgrid = gbs::make_structuredgrid(pts, ni, nj);
    SolverCase<T> solver_case;
    solver_case.gi = make_grid_info<T>(sgrid);

    auto &g = *(solver_case.gi->g);
    auto &g_metrics = *(solver_case.gi->g_metrics);

    auto  f_Vm    = [](const auto &gp){return sin(gp.m)*gp.l + cos(gp.l)*gp.m;};
    auto df_Vm_dm = [](const auto &gp){return  cos(gp.m)*gp.l + cos(gp.l);};
    auto df_Vm_dl = [](const auto &gp){return -sin(gp.l)*gp.m + sin(gp.m);};
    std::for_each(
        g.begin(), g.end(),
        [&f_Vm](auto &gp){gp.Vm = f_Vm(gp);}
    );

    for(size_t i{}; i < ni; i++)
    {
        for(size_t j{}; j < nj; j++)
        {
            auto df_Vm_dm_ = D1_O2_dx1(g, g_metrics, i, j, d_ksi, d_eth, [](const auto &gp)
                                       { return gp.Vm; });
            auto df_Vm_dl_ = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, [](const auto &gp)
                                      { return gp.Vm; }); 
            if (PLOT_ON)
            {
                g(i, j).Vu = df_Vm_dm_; // use variable only for storage
                g(i, j).H  = df_Vm_dl_;
            }
            ASSERT_NEAR(df_Vm_dm_, df_Vm_dm(g(i, j)), 1e-3);
            ASSERT_NEAR(df_Vm_dl_, df_Vm_dl(g(i, j)), 1e-3);
        }
    }

    if(PLOT_ON)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "m", [](const auto &gp)
                  { return gp.m; });
        add_value(g, structuredGrid, "l", [](const auto &gp)
                  { return gp.l; });
        add_value(g, structuredGrid, "Vm", [](const auto &gp)
                  { return gp.Vm; });
        add_value(g, structuredGrid, "dVm_dm", [&df_Vm_dm](const auto &gp)
                  { return df_Vm_dm(gp); });
        add_value(g, structuredGrid, "dVm_dl", [&df_Vm_dl](const auto &gp)
                  { return df_Vm_dl(gp); });
        add_value(g, structuredGrid, "dVm_dm_err", [&](const auto &gp)
                  { return df_Vm_dm(gp) - gp.Vu; });
        add_value(g, structuredGrid, "dVm_dl_err", [&](const auto &gp)
                  { return df_Vm_dl(gp) - gp.H; });
        plot_vtkStructuredGrid(structuredGrid, "Vm", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dm", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dl", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dm_err", false);
        plot_vtkStructuredGrid(structuredGrid, "dVm_dl_err", false);
        auto sgrid_actor = gbs::make_structuredgrid_actor(pts, ni, nj);
        gbs::plot(
            crv_lst
            , sgrid_actor
        );
    }

}