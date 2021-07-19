#include <gtest/gtest.h>
#include <curvaturesolver.h>
#include <plots.h>
#include <gridreader.h>
#include <gridsbuilders.h>
#include <gridrender.h>
#include <datastorage.h>
#include <vtkXMLStructuredGridWriter.h>
// #include <CoolPropLib.h>

#include <chrono>
const bool TESTS_USE_PLOT = true;
TEST(tests_curvature_solver, vtk_no_blades)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/out/test_001.vts");
    auto Vm = 30.;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; });
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
        // .RF = 0.01
    };

    quiss::SolverCase<T> solver_case{
        .gi = gi,
    };

    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        if (TESTS_USE_PLOT)
        {
            auto structuredGrid = quiss::make_vtkStructuredGrid(g);
            add_value(g, structuredGrid, "Vm", [](const auto &gp)
                      { return gp.Vm; });
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_residual(solver_case.log);

            vtkNew<vtkXMLStructuredGridWriter> writer;
            writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm_rho_cst.vts");
            writer->SetInputData(structuredGrid);
            writer->Write();
        }
    }

    {
        
        gi.rho_cst = false;
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        if (TESTS_USE_PLOT)
        {
            auto structuredGrid = quiss::make_vtkStructuredGrid(g);
            add_value(g, structuredGrid, "Vm", [](const auto &gp)
                      { return gp.Vm; });
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);

            vtkNew<vtkXMLStructuredGridWriter> writer;
            writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm_rho_var.vts");
            writer->SetInputData(structuredGrid);
            writer->Write();
        }
    }


    {
        gi.rho_cst = true;
        // std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm*0.5;gp.Vu=Vm*0.5;gp.H=gp.Cp*gp.Tt;gp.Pt=1.6432411e5; });
        std::for_each(g.begin(0), g.end(0), [Vm, &gi](auto &gp)
            {
                gp.Vm = Vm * 0.5;
                gp.Vu = Vm * 0.5;
                // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
                // gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
                gp.Ts = 300.;
                gp.Ps = 1.6432411e5;
                gp.rho = gp.Ps / (gi.R) / gp.Ts;
                // gp.Tt = 300.;
                // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
                gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / gp.Cp;
                gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
                gp.H = gp.Tt * gp.Cp;
                gp.s = std::log(pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04));
            });
        std::for_each(g.begin(), g.end(), [&Vm,&gi](auto &gp) {gp.rho= 1.6432411e5 / (gi.R) / 300.; });
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        if (TESTS_USE_PLOT)
        {
            auto structuredGrid = quiss::make_vtkStructuredGrid(g);
            add_value(g, structuredGrid, "Vm", [](const auto &gp)
                      { return gp.Vm; });
            add_value(g, structuredGrid, "Vu", [](const auto &gp)
                      { return gp.Vu; });
            add_value(g, structuredGrid, "rho", [](const auto &gp)
                      { return gp.rho; });
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);

            vtkNew<vtkXMLStructuredGridWriter> writer;
            writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm_swirl_rho_cst.vts");
            writer->SetInputData(structuredGrid);
            writer->Write();
        }
    }

    // {
    //     gi.rho_cst = false;
    //     // std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm*0.5;gp.Vu=Vm*0.5;gp.H=gp.Cp*gp.Tt;gp.Pt=1.6432411e5; });
    //     auto start = high_resolution_clock::now();
    //     quiss::curvature_solver(solver_case);
    //     auto stop = high_resolution_clock::now();
    //     auto duration = duration_cast<microseconds>(stop - start);
    //     cout << "Time taken by meridian computation: "
    //          << duration.count() << " microseconds" << endl;

    //     if (TESTS_USE_PLOT)
    //     {
    //         auto structuredGrid = quiss::make_vtkStructuredGrid(g);
    //         add_value(g, structuredGrid, "Vm", [](const auto &gp)
    //                   { return gp.Vm; });
    //         add_value(g, structuredGrid, "Vu", [](const auto &gp)
    //                   { return gp.Vu; });
    //         structuredGrid->GetPointData()->SetActiveScalars("Vu");
    //         quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);

    //         vtkNew<vtkXMLStructuredGridWriter> writer;
    //         writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm_swirl_rho_var.vts");
    //         writer->SetInputData(structuredGrid);
    //         writer->Write();
    //     }
    // }

    // {
    //     gi.rho_cst = false;
    //     solver_case.max_geom=1;
    //     std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm*0.5;gp.Vu=Vm*0.5;gp.H=gp.Cp*gp.Tt;gp.Pt=1.6432411e5; });
    //     auto r1 = g(0,0).y;
    //     auto r2 = g(0,nj-1).y;
    //     std::for_each(g.begin(0), g.end(0), [r1,r2](auto &gp) {gp.Tt = 300. * (gp.y - r2) /(r1 -r2) - 310. * (gp.y - r1) /(r1 -r2);gp.H = gp.Tt * gp.Cp - 1004. * 288.;});
    //     auto start = high_resolution_clock::now();
    //     quiss::curvature_solver(solver_case);
    //     auto stop = high_resolution_clock::now();
    //     auto duration = duration_cast<microseconds>(stop - start);
    //     cout << "Time taken by meridian computation: "
    //          << duration.count() << " microseconds" << endl;

    //     if (TESTS_USE_PLOT)
    //     {
    //         auto structuredGrid = quiss::make_vtkStructuredGrid(g);
    //         add_value(g, structuredGrid, "Vm", [](const auto &gp)
    //                   { return gp.Vm; });
    //         add_value(g, structuredGrid, "Vu", [](const auto &gp)
    //                   { return gp.Vu; });
    //         add_value(g, structuredGrid, "Tt", [](const auto &gp)
    //                   { return gp.Tt; });
    //         add_value(g, structuredGrid, "Ts", [](const auto &gp)
    //                   { return gp.Ts; });
    //         add_value(g, structuredGrid, "rho", [](const auto &gp)
    //                   { return gp.rho; });
    //         quiss::plot_vtkStructuredGrid(structuredGrid,"rho", true);

    //         vtkNew<vtkXMLStructuredGridWriter> writer;
    //         writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm_swirl_rho_var_Tt_ramp.vts");
    //         writer->SetInputData(structuredGrid);
    //         writer->Write();
    //     }
    // }
}

TEST(tests_curvature_solver, vtk_static_blades1)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_002.vts");
    auto Vm = 30.;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; });
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
        .nj = nj
    };

    quiss::SolverCase<T> solver_case{
        .gi = gi,
        .max_geom=1};

    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
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

        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_002.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_static_blades2)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_003.vts");
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
        .RF = 0.06,
    };

    quiss::SolverCase<T> solver_case{
        .gi = gi,
        .max_geom = 1000
    };

    // Inlet
    auto r1 = g(0, 0).y;
    auto r2 = g(0, nj - 1).y;
    // std::for_each(g.begin(0), g.end(0), [r1,r2](auto &gp) 
    //     {
    //         // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
    //         gp.Tt = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
    //         // gp.Tt = 300.;
    //         // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
    //         gp.Ts = gp.Tt - (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
    //         gp.Ps = gp.Pt - (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
    //         gp.H = gp.Tt * gp.Cp;
    //         gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
    //     }
    // );
    std::for_each(g.begin(0), g.end(0), [r1,r2,Ps,&gi](auto &gp) 
        {
            // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
            gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Ps = Ps;
            gp.rho= gp.Ps / (gi.R) / gp.Ts;
            // gp.Tt = 300.;
            // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
            gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
            gp.H = gp.Tt * gp.Cp;
            gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
        }
    );


    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_003.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_vtkStructuredGrid(structuredGrid,"r", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_static_blades3)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_004.vts");
    auto Vm = 30.;
    auto Ps = 1.2e5;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt;gp.Pt=133337.02; /*gp.iB=-1;*/  if(gp.iB!=-1) gp.omg_=0.0; });
    // std::cout << PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]") << std::endl;

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
        .RF = 0.06,
    };

    quiss::SolverCase<T> solver_case{
        .gi = gi,
        .max_geom = 1000
    };

    // Inlet
    auto r1 = g(0, 0).y;
    auto r2 = g(0, nj - 1).y;
    // std::for_each(g.begin(0), g.end(0), [r1,r2](auto &gp) 
    //     {
    //         // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
    //         gp.Tt = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
    //         // gp.Tt = 300.;
    //         // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
    //         gp.Ts = gp.Tt - (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
    //         gp.Ps = gp.Pt - (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
    //         gp.H = gp.Tt * gp.Cp;
    //         gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
    //     }
    // );
    std::for_each(g.begin(0), g.end(0), [r1,r2,Ps,&gi](auto &gp) 
        {
            // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
            gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            // gp.Ts = 300.;
            gp.Ps = Ps;
            gp.rho= gp.Ps / (gi.R) / gp.Ts;
            // gp.Tt = 300.;
            // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
            gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
            gp.H = gp.Tt * gp.Cp;
            gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
        }
    );


    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;

        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_004.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_vtkStructuredGrid(structuredGrid,"r", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_stream_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_005.vts");
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
        .max_geom = 5000
    };

    // Inlet
    auto r1 = g(0, 0).y;
    auto r2 = g(0, nj - 1).y;
    std::for_each(g.begin(0), g.end(0), [r1,r2,Ps,Vm,&gi](auto &gp) 
        {
            gp.Vm=Vm;
            gp.Vu=Vm;
            // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
            gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            // gp.Ts = 300.;
            gp.Ps = Ps;
            gp.rho= gp.Ps / (gi.R) / gp.Ts;
            // gp.Tt = 300.;
            // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
            gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
            gp.H = gp.Tt * gp.Cp;
            gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
        }
    );


    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_005.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_vtkStructuredGrid(structuredGrid,"r", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_span_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_006.vts");
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
        .max_geom = 5000
    };

    // Inlet
    auto r1 = g(0, 0).y;
    auto r2 = g(0, nj - 1).y;
    std::for_each(g.begin(0), g.end(0), [r1,r2,Ps,Vm,&gi](auto &gp) 
        {
            gp.Vm=Vm;
            gp.Vu=Vm;
            // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
            gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            // gp.Ts = 300.;
            gp.Ps = Ps;
            gp.rho= gp.Ps / (gi.R) / gp.Ts;
            // gp.Tt = 300.;
            // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
            gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
            gp.H = gp.Tt * gp.Cp;
            gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
        }
    );


    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_006.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_vtkStructuredGrid(structuredGrid,"r", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
        }
    }
}

TEST(tests_curvature_solver, vtk_non_otrtho_channel_mixed_dir)
{
    using T = double;
    using namespace std::chrono;

    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/in/test_007.vts");
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
        .max_geom = 5000
    };

    // Inlet
    auto r1 = g(0, 0).y;
    auto r2 = g(0, nj - 1).y;
    std::for_each(g.begin(0), g.end(0), [r1,r2,Ps,Vm,&gi](auto &gp) 
        {
            gp.Vm=Vm;
            gp.Vu= Vm * 0.8 * (gp.y - r2) /(r1 -r2) - Vm * 1.2 * (gp.y - r1) /(r1 -r2);
            // gp.Tt = 350. * (gp.y - r2) /(r1 -r2) - 300. * (gp.y - r1) /(r1 -r2);
            gp.Ts = 300. + 50. * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            // gp.Ts = 300.;
            gp.Ps = Ps;
            gp.rho= gp.Ps / (gi.R) / gp.Ts;
            // gp.Tt = 300.;
            // gp.Pt = 133337.02 + 0.2e5 * std::sin( (gp.y - r2) /(r1 -r2) * std::numbers::pi );
            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / 1004.;
            gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
            gp.H = gp.Tt * gp.Cp;
            gp.s = std::log( pow(gp.Ts / 288., gp.Cp) / std::pow(gp.Ps / 1.01325e5, 287.04) );
        }
    );


    {
        auto start = high_resolution_clock::now();
        quiss::curvature_solver(solver_case);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by meridian computation: "
             << duration.count() << " microseconds" << endl;


        auto structuredGrid = quiss::write_vtk_grid(g,"C:/Users/sebastien/workspace/tbslib/tests/out/test_007.vts");
        if (TESTS_USE_PLOT)
        {
            quiss::plot_residual(solver_case.log);
            quiss::plot_vtkStructuredGrid(structuredGrid,"r", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Ts", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vm", true);
            quiss::plot_vtkStructuredGrid(structuredGrid,"Vu", true);
        }
    }
}