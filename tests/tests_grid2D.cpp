#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-render/vtkcurvesrender.h>
#include <gbs-render/vtkgridrender.h>
#include <meshtools.h>
#include <meridionalsolvercase.h>
#include <gridmetrics.h>
#include <gridrender.h>
#include <gridsbuilders.h>
#include <numbers>
using T = double;
using namespace yams;
const bool PLOT_ON = false;
const T tol_ang{static_cast<T>(1e-3)};
TEST(grid_2d, grid_angle_phi)
{
    const T phi{std::numbers::pi_v<T>/3.};
    std::vector<std::array<T,2>> poles1{
            {0.0,0.0},
            {1.0,std::tan(phi)},
        };
    std::vector<std::array<T,2>> poles2{
            {0.0,1.0},
            {1.0,1.0 + std::tan(phi)},
        };
    std::vector<T> knots{0.,1.};
    std::vector<size_t> mult{2,2};
    size_t deg{1};

    auto crv1 = std::make_shared<gbs::BSCurve<T,2>>(
        poles1,
        knots,
        mult,
        deg
    );
    auto crv2 = std::make_shared<gbs::BSCurve<T,2>>(
        poles2,
        knots,
        mult,
        deg
    );
    size_t nu{30};
    size_t nv{15};
    yams::crv_vector<T> crv_lst{crv1, crv2};
    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = mesh_channel<T>(crv_lst, knots, nv, nu);
    auto sgrid = gbs::make_structuredgrid(pts, ni, nj);
    SolverCase<T> solver_case;
    solver_case.gi = make_grid_info<T>(sgrid);

    auto &g = *(solver_case.gi->g);
    compute_abscissas(g);
    compute_angles(g);

    for(size_t i{}; i < ni; i++)
    {
        for(size_t j{}; j < nj; j++)
        {
            ASSERT_NEAR(phi * 180 / std::numbers::pi_v<T> ,g(i, j).phi * 180 / std::numbers::pi_v<T>, tol_ang);
            ASSERT_NEAR( 0. ,g(i, j).gam * 180 / std::numbers::pi_v<T>, tol_ang);
        }
    }


    if(PLOT_ON)
    {
        sgrid = make_vtk_grid(g);
        plot_vtkStructuredGrid(sgrid,"phi");
        auto sgrid_actor = gbs::make_structuredgrid_actor(pts, ni, nj);
        gbs::plot(
            crv_lst, 
            sgrid_actor
        );
    }
}

TEST(grid_2d, grid_angle_gam)
{
    const T gam{std::numbers::pi_v<T>/3.};
    std::vector<std::array<T,2>> poles1{
            {0.0,0.0},
            {1.0,0.0},
        };
    std::vector<std::array<T,2>> poles2{
            {0.0 + std::tan(gam),1.0},
            {1.0 + std::tan(gam),1.0},
        };
    std::vector<T> knots{0.,1.};
    std::vector<size_t> mult{2,2};
    size_t deg{1};

    auto crv1 = std::make_shared<gbs::BSCurve<T,2>>(
        poles1,
        knots,
        mult,
        deg
    );
    auto crv2 = std::make_shared<gbs::BSCurve<T,2>>(
        poles2,
        knots,
        mult,
        deg
    );
    size_t nu{30};
    size_t nv{15};
    yams::crv_vector<T> crv_lst{crv1, crv2};
    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = mesh_channel<T>(crv_lst, knots, nv, nu);
    auto sgrid = gbs::make_structuredgrid(pts, ni, nj);
    SolverCase<T> solver_case;
    solver_case.gi = make_grid_info<T>(sgrid);

    auto &g = *(solver_case.gi->g);
    compute_abscissas(g);
    compute_angles(g);

    for(size_t i{}; i < ni; i++)
    {
        for(size_t j{}; j < nj; j++)
        {
            ASSERT_NEAR(gam * 180 / std::numbers::pi_v<T> ,g(i, j).gam * 180 / std::numbers::pi_v<T>, tol_ang);
        }
    }


    if(PLOT_ON)
    {
        sgrid = make_vtk_grid(g);
        plot_vtkStructuredGrid(sgrid,"gam");
        auto sgrid_actor = gbs::make_structuredgrid_actor(pts, ni, nj);
        gbs::plot(
            crv_lst, 
            sgrid_actor
        );
    }
}

TEST(grid_2d, grid_angle_gam_phi)
{

    size_t ni = 5000;
    size_t nj = 1500;
    MeridionalGrid<T> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto t1 = std::numbers::pi_v<T>;
    auto t2 = std::numbers::pi_v<T> * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    compute_angles(g);

    for (size_t i{}; i < ni ; i++)
    {
        auto phi = -std::numbers::pi_v<T> / 2 + std::numbers::pi_v<T> * i / (ni - 1.);
        auto gam = 3 * std::numbers::pi_v<T> / 2 - std::numbers::pi_v<T> * i / (ni - 1.);
        for (size_t j{}; j < nj ; j++)
        {
            const auto &gp = g(i, j);
            ASSERT_NEAR(fmod(gam * 180 / std::numbers::pi_v<T> + 360, 360), fmod(gp.gam * 180 / std::numbers::pi_v<T> + 360, 360), 0.1);
            ASSERT_NEAR(fmod(phi * 180 / std::numbers::pi_v<T> + 360, 360), fmod(gp.phi * 180 / std::numbers::pi_v<T> + 360, 360), 0.1);
            ASSERT_NEAR(0., gp.sgp, 0.01);
            ASSERT_NEAR(1., gp.cgp, 0.01);
        }
    }
    if (PLOT_ON)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "phi", [](const auto &gp)
                  { return gp.phi * 180 / std::numbers::pi_v<T>; });
        add_value(g, structuredGrid, "gam", [](const auto &gp)
                  { return fmod(gp.gam * 180 / std::numbers::pi_v<T> + 360, 360); });
        add_value(g, structuredGrid, "cgp", [](const auto &gp)
                  { return gp.cgp; });
        add_value(g, structuredGrid, "sgp", [](const auto &gp)
                  { return gp.sgp; });
        plot_vtkStructuredGrid(structuredGrid, "phi");
        plot_vtkStructuredGrid(structuredGrid, "gam");
        plot_vtkStructuredGrid(structuredGrid, "cgp");
        plot_vtkStructuredGrid(structuredGrid, "sgp");
    }
}