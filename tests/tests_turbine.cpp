#include <gtest/gtest.h>
#include <curvaturesolver.h>
#include <plots.h>
#include <gridreader.h>
#include <gridsbuilders.h>
#include <gridrender.h>
#include <datastorage.h>
#include <vtkXMLStructuredGridWriter.h>


#ifdef _WIN32
    const std::string test_files_path = "../tests/";
#else
    const std::string test_files_path = "../../../tbslib/tests/";
#endif

#include <chrono>
using namespace yams;
const bool TESTS_USE_PLOT = true;

TEST(tests_turbine, vtk_no_blades)
{
    using T = double;
    using namespace std::chrono;

    auto g = read_vtk_grid<T>(test_files_path+"turbine_grid.vts");

    auto Vm = 104.;
    // auto dH = 1004. * 10.;
    auto ga = 1.32;
    auto Pt = 3.8e5;
    auto Tt = 1100.;
    auto rg = 287.04;
    auto cp = rg * ga /( ga - 1 );
    auto Ts = Tt - Vm*Vm / 2. / cp;
    auto a  = std::sqrt(ga*rg*Ts);
    auto M  = Vm / a;
    auto Ps = Pt / std::pow(1+0.5*(ga-1)*M*M, ga / ( ga - 1 ) );
    size_t max_geom=5;

    // init values
    std::for_each(g.begin(), g.end(), [Vm, cp, ga, rg, Tt, Pt, Ts, Ps](auto &gp) {gp.Vm=Vm; gp.Cp=cp, gp.ga=ga; gp.Vu=0.; gp.H=cp*Tt; gp.Pt=Pt; gp.Tt=Tt; gp.rho=Ps/rg/Ts;});
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    Array2d<Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    compute_grid_metrics(g,g_metrics,f_m,f_l);// TODO run in //

    GridInfo<T> gi{
        .g = std::make_shared< MeridionalGrid<T> >( g ),
        .g_metrics = std::make_shared< Grid2dMetrics<T> >( g_metrics ),
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.05,
    };

    SolverCase<T> solver_case{
        .gi = std::make_shared< GridInfo<T> >( gi ),
        .inlet = InletBC<T>{
            .mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [Ts](auto l_rel){return Ts;},
            .Vm_moy=Vm
        },
        .max_geom = max_geom,
        // .relocate=false,
        .mf_ref_span = ni / 2,
    };
    
    {
        // shadowing
        auto &gi= *(solver_case.gi);
        auto &g = *(solver_case.gi->g);
        if(solver_case.use_meridional_grad)
        {
            // solver_case.gi->RF /= 10.;
            // solver_case.max_geom = 1000;
        }
        auto start = high_resolution_clock::now();
        curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = make_vtk_grid(g);
        if (TESTS_USE_PLOT)
        {
            plot_vtkStructuredGrid(structuredGrid,"rho", true);
            // plot_residual(solver_case.log);
        }
    }

}

#include <array>
#include <iostream>
#include <list>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

void print(auto const rem, auto const& range)
{
    for (std::cout << rem; auto const& elem : range)
        std::cout << elem << ' ';
    std::cout << '\n';
}

TEST(tests_turbine, zip)
{
    auto x = std::vector{ 1, 2, 3, 4 };
    auto y = std::list<std::string>{ "alpha", "beta", "gamma", "delat", "epsilon" };
    auto z = std::array{ 'A', 'B', 'C', 'D', 'E', 'F' };

    print("Source views:", "");
    print("x: ", x);
    print("y: ", y);
    print("z: ", z);

    print("\nzip(x,y,z):", "");

    for (auto elem : std::views::zip(x, y, z))
    {
        auto &[x_, y_, z_] = elem;
        std::cout << x_ << ' ' << y_ << ' ' << z_ << std::endl;

        std::get<char&>(elem) += ('a' - 'A'); // modifies the element of z
    }

    print("\nAfter modification, z: ", z);
}