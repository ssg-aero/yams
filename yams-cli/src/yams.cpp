#include <yams/datastorage.h>
#include <yams/gridreader.h>
#include <yams/curvaturesolver.h>
#include <yams/plots.h>
#include <yams/gridrender.h>

#include <string>


int main(int argc, char *argv[])
{
    using T = double;
    auto Vm = 100.;

    // std::string fname { "../../../cosapp-tubomachine/tests/data/cmp_ax_1_stage_cad" };
    // std::string fname { "../../../tbslib/tests/in/test_009" };
    std::string fname { argv[1] };
    std::cout << "opening " << fname << std::endl;
    auto g = yams::read_vtk_grid<T>( (fname+".vts").c_str() );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = std::make_shared< yams::MeridionalGrid<T> >( g ),
        .g_metrics = std::make_shared< yams::Grid2dMetrics<T> >( g_metrics ),
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.01,
    };

    yams::SolverCase<T> solver_case{
        .gi = std::make_shared< yams::GridInfo<T> >( gi ),
        .inlet = yams::InletBC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            // .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            // .Vu   = [Vm](auto l_rel){return  Vm * 0.5 * (1.-l_rel) + Vm * 1.5 * l_rel;},
            .Vm_moy=Vm
        },
    };
    std::cout << "opening blade info" << (fname+"_bld.json").c_str() << std::endl;
    yams::read_blade_info( (fname+"_bld.json").c_str(), solver_case );
    std::cout << "Done.\n";
    solver_case.max_geom = 300;
    // solver_case.tol_rel_mf=0.05;
    // solver_case.tol_rel_pos=0.05;

    std::for_each(g.begin(), g.end(), [&Vm](auto &gp)
                  {
                      if (gp.iB != -1)
                          gp.omg_ = 0.01;
                  });
    std::cout << "launch computation.\n";
    using namespace std::chrono;
    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,(fname+"_design.vts").c_str());
        // if (TESTS_USE_PLOT)
        {
            // yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"Pt", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"Tt", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"V", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"rho", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"bet", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"alf", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"s", true);
            // yams::plot_vtkStructuredGrid(structuredGrid,"H", true);
            yams::plot_residual(solver_case.log);
        }
    }
}