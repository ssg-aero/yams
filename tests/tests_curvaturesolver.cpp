#include <gtest/gtest.h>
#include <curvaturesolver.h>
#include <plots.h>
#include <gridreader.h>
#include <gridsbuilders.h>
#include <gridrender.h>
#include <datastorage.h>
#include <vtkXMLStructuredGridWriter.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/writer.h>
#include <tbs/io/json2.h>
#include <tbs/io/iges.h>
#include <tbs/ui/plots.h>

// #include <CoolPropLib.h>

#include <chrono>
const bool TESTS_USE_PLOT = true;
TEST(tests_curvature_solver, vtk_no_blades)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/out/test_001.vts");
    auto Vm = 30.;
    auto dH = 1004. * 10.;
    auto Pt = 133337.02;
    auto Tt = 300.;
    auto Ps = Pt - Vm*Vm / 2. / g(0,0).rho;
    auto Ts = Tt - Vm*Vm / 2. / g(0,0).Cp;
    size_t max_geom=500;
    // init values
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; });
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .RF = 0.05
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [Ts](auto l_rel){return Ts;},
            .Vm_moy=Vm
        },
        // .max_geom = 2
    };
    
    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_001_Vm_rho_cst.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_residual(solver_case.log);
        }
    }
// return;
    {
        
        gi.rho_cst = false;
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_001_Vm_rho_var.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_residual(solver_case.log);
        }
    }



    {
        gi.rho_cst = true;
        solver_case.inlet.Vu = [Vm](auto l_rel){return 0.5*Vm;};
        solver_case.inlet.Vm_moy = 0.5* Vm;

        std::for_each(g.begin(), g.end(), [&Vm,&gi](auto &gp) {gp.rho= 1.6432411e5 / (gi.R) / 300.; });
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_001_Vm_swirl_rho_cst.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }

    {
        gi.rho_cst = false;
        // std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm*0.5;gp.Vu=Vm*0.5;gp.H=gp.Cp*gp.Tt;gp.Pt=1.6432411e5; });
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_001_Vm_swirl_rho_var.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }

    {
        gi.rho_cst = false;
        // solver_case.max_geom=1;
        std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm*0.5;gp.Vu=Vm*0.5;gp.H=gp.Cp*gp.Tt;gp.Pt=1.6432411e5; });
        solver_case.inlet.Vu = [Vm](auto l_rel){return 0.5*Vm;};
        solver_case.inlet.Ps = [](auto l_rel){return 1.6432411e5;};
        solver_case.inlet.Ts = [](auto l_rel){return 300. * (1. - l_rel) + 310 * l_rel;};
        // auto r1 = g(0,0).y;
        // auto r2 = g(0,nj-1).y;
        // std::for_each(g.begin(0), g.end(0), [r1,r2](auto &gp) {gp.Tt = 300. * (gp.y - r2) /(r1 -r2) - 310. * (gp.y - r1) /(r1 -r2);gp.H = gp.Tt * gp.Cp - 1004. * 288.;});

        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_001_Vm_swirl_rho_var_Tt_ramp.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_residual(solver_case.log);
        }
    }
}

TEST(tests_curvature_solver, vtk_static_blades1)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_002.vts");
    auto Vm = 30.;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; });
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .bld_info_lst {yams::BladeInfo<T>{
            .mode=yams::MeridionalBladeMode::DIRECT,
            }},
        .max_geom=1,
        };

    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        size_t i_bf = 89;

        auto ri = g(i_bf,0).y;
        auto Vmi= g(i_bf,0).Vm;
        auto alp= g(i_bf,0).bet;

        for(size_t j {} ; j < nj ; j++)
        {
            auto r = g(i_bf,j).y;
            auto Vm= g(i_bf,j).Vm;
            auto res = Vm / Vmi;
            auto res_analytic = std::pow( ri / r, std::pow( std::sin(alp), 2 ) );
            auto err_pc = (res - res_analytic) / res_analytic * 100.;
            ASSERT_LT(err_pc,1.);
            std::cout << res << " " << res_analytic << " " << err_pc << "%" << std::endl;
        }

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_002.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_static_blades2)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_003.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    // auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.1; });
    // std::cout << PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]") << std::endl;

    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.05,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .bld_info_lst {yams::BladeInfo<T>{
            .mode=yams::MeridionalBladeMode::DIRECT,
            }},
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [](auto l_rel){return 300. + 50. * std::sin( l_rel * std::numbers::pi );},
            .Vm_moy=Vm
        },
        .max_geom = 1000
    };

    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_003.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }
}

TEST(tests_curvature_solver, vtk_static_blades3)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_004.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.1; });
    // std::cout << PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]") << std::endl;

    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.06,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .bld_info_lst {yams::BladeInfo<T>{
            .mode=yams::MeridionalBladeMode::DIRECT,
            }},
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [](auto l_rel){return 300. + 50. * std::sin( l_rel * std::numbers::pi );},
            .Vm_moy=Vm
        },
        .max_geom = 1000,
    };

    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_004.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Pt", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_stream_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_005.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [Ps,Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Ps=Ps; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.0; });

    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.025,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            .Vm_moy=Vm
        },
        .max_geom = 2000
    };

    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_005.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_span_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_006.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.0; });

    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.05,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            .Ps   = [Ps](auto l_rel){return Ps;},
            .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            .Vm_moy=Vm
        },
        .max_geom = 5000
    };


    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_006.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_residual(solver_case.log);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_mixed_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = yams::read_vtk_grid<T>("../../../tbslib/tests/in/test_007.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    // auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.0; });

    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.025,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            .Vu   = [Vm](auto l_rel){return  Vm * 0.8 * (1.-l_rel) + Vm * 1.2 * l_rel;},
            .Vm_moy=Vm
        },
        .max_geom = 5000,
        .tol_rel_pos = 2e-5,
    };


    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_007.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"rho", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"s", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"H", true);
            yams::plot_residual(solver_case.log);
        }
    }

}

TEST(tests_curvature_solver, vtk_read_blade_info)
{
    using T = double;
    auto Vm = 30.;

    std::string fname { "../../../tbslib/tests/in/test_004" };
    auto g = yams::read_vtk_grid<T>( (fname+".vts").c_str() );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.025,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            .Vu   = [Vm](auto l_rel){return  Vm * 0.5 * (1.-l_rel) + Vm * 1.5 * l_rel;},
            .Vm_moy=Vm
        },
    };

    yams::read_blade_info( (fname+"_bld.json").c_str(), solver_case );


    using namespace std::chrono;
    {
        auto start = high_resolution_clock::now();
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,"../../../tbslib/tests/out/test_004_design.vts");
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"V", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"rho", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"bet", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"s", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"H", true);
            yams::plot_residual(solver_case.log);
        }
    }

}


TEST(tests_curvature_solver, vtk_fan_design)
{
    using T = double;
    auto Vm = 100.;

    std::string fname { "../../../tbslib/tests/in/test_008" };
    auto g = yams::read_vtk_grid<T>( (fname+".vts").c_str() );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.025,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            // .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            // .Vu   = [Vm](auto l_rel){return  Vm * 0.5 * (1.-l_rel) + Vm * 1.5 * l_rel;},
            .Vm_moy=Vm
        },
    };

    yams::read_blade_info( (fname+"_bld.json").c_str(), solver_case );
    solver_case.max_geom = 5000;
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
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,(fname+"_design.vts").c_str());
        if (TESTS_USE_PLOT)
        {
            yams::plot_vtkStructuredGrid(structuredGrid,"Ps", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Pt", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Tt", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"V", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"rho", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"bet", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"s", true);
            yams::plot_vtkStructuredGrid(structuredGrid,"H", true);
            yams::plot_residual(solver_case.log);
        }
    }

    rapidjson::Document document;
    gbs::parse_file((fname+".json").c_str(), document);

    {
        tbs::AeroMachine<double> b_set;
        tbs::read_and_store_aero_machine(b_set, document);
        tbs::plot_3d(b_set);
    }

    auto i1 = solver_case.bld_info_lst[0].i1;
    auto i2 = solver_case.bld_info_lst[0].i2;
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][0]["k"][0] = g(i1,0).bet;
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][0]["k"][1] = g(i2,0).bet;
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][0]["gauge"]= 0.5 *( g(i1,0).bet + g(i2,0).bet);
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][1]["k"][0] = g(i1,nj-1).bet;
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][1]["k"][1] = g(i2,nj-1).bet;
    document["channel_sets"][0]["channels"][1]["blades"][0]["sections"][1]["gauge"]= 0.5 *( g(i1,nj-1).bet + g(i2,nj-1).bet);

    FILE *fp = fopen((fname+"_design.json").c_str(), "wb"); // non-Windows use "w"

    char writeBuffer[65536];
    rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));

    rapidjson::Writer<rapidjson::FileWriteStream> writer(os);
    document.Accept(writer);

    fclose(fp);

    {
        tbs::AeroMachine<double> b_set;
        tbs::read_and_store_aero_machine(b_set, document);
        tbs::plot_3d(b_set);
        tbs::export_iges(b_set,(fname+"_design.iges").c_str());
    }
}

TEST(tests_curvature_solver, vtk_fan_ogv_design)
{
    using T = double;
    auto Vm = 100.;

    std::string fname { "../../../tbslib/tests/in/test_009" };
    auto g = yams::read_vtk_grid<T>( (fname+".vts").c_str() );

    size_t ni = g.nRows();
    size_t nj = g.nCols();
    
    double ksi = 1. / (ni-1.);
    double eth = 1. / (nj-1.);

    yams::Array2d<yams::Grid2dMetricsPoint<T>>   g_metrics(ni,nj);
    yams::compute_grid_metrics(g,g_metrics,yams::f_m,yams::f_l);// TODO run in //

    yams::GridInfo<T> gi{
        .g = g,
        .g_metrics = g_metrics,
        .d_ksi = ksi,
        .d_eth = eth,
        .ni = ni,
        .nj = nj,
        .rho_cst=false,
        .RF = 0.01,
    };

    yams::SolverCase<T> solver_case{
        .gi = gi,
        .inlet = yams::Inlet_BC<T>{
            .mode = yams::MeridionalBC::INLET_VmMoy_Ts_Ps_Vu,
            // .Ps   = [Ps](auto l_rel){return 1e5+0.1e5*l_rel;}, //TODO investigate this mess
            // .Ts   = [](auto l_rel){return 300. +50. * std::sin( l_rel * std::numbers::pi );},
            // .Vu   = [Vm](auto l_rel){return  Vm * 0.5 * (1.-l_rel) + Vm * 1.5 * l_rel;},
            .Vm_moy=Vm
        },
    };

    yams::read_blade_info( (fname+"_bld.json").c_str(), solver_case );
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
        yams::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = yams::write_vtk_grid(g,(fname+"_design.vts").c_str());
        if (TESTS_USE_PLOT)
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