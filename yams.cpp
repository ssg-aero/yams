#include <string>
#include <datastorage.h>
#include <gridreader.h>
#include <curvaturesolver.h>
#include <plots.h>
#include <gridrender.h>

int main(int argc, char *argv[])
{
    using T = double;
    auto Vm = 100.;

    // std::string fname { "C:/Users/sebastien/workspace/tbslib/tests/in/test_009" };
    std::string fname { argv[1] };
    auto g = quiss::read_vtk_grid<T>( (fname+".vts").c_str() );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    quiss::Array2d<quiss::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    quiss::compute_grid_metrics(g,g_metrics,quiss::f_m,quiss::f_l);// TODO run in //

    quiss::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.01,
    };

    quiss::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = quiss::Inlet_BC<T>{
            .mode = quiss::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            // .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            // .Vu   = [Vm](auto l_rel){return  Vm * 0.5 * (1.-l_rel) + Vm * 1.5 * l_rel;},
            .Vm_moy=Vm
        },
    };

    quiss::read_blade_info( (fname+"_bld.json").c_str(), solver_case );
    solver_case.max_geom = 300;
    // solver_case.tol_rel_mf=0.05;
    // solver_case.tol_rel_pos=0.05;

    std::for_each(g.begin(), g.end(), [&Vm](auto &gp)
                  {
                      if (gp.iB != -1)
                          gp.omg_ = 0.01;
                  });

    using namespace std::chrono;
    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = quiss::write_vtk_grid(g,(fname+"_design.vts").c_str());
        // if (TESTS_USE_PLOT)
        {
            // quiss::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"Pt", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"Tt", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"V", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"rho", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"bet", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"alf", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"s", true);
            // quiss::plot_vtkStructuredGrid(structuredGrid,"H", true);
            quiss::plot_residual(solver_case.log);
        }
    }
}