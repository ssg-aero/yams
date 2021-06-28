#include <gtest/gtest.h>
#include <curvaturesolver.h>
#include <gridreader.h>
#include <gridsbuilders.h>
#include <gridrender.h>
#include <datastorage.h>
#include <vtkXMLStructuredGridWriter.h>
const bool TESTS_USE_PLOT = true;
TEST(tests_curvature_solver, vtk_no_blades)
{
    using T = double;
    // auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/out/alpx001.vts");
    auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/out/test_001.vts");
    // auto g = quiss::read_vtk_grid<T>("C:/Users/sebastien/workspace/tbslib/tests/out/test_002.vts");
    auto Vm = 30.;
    auto dH = 1004. * 10.;
    size_t max_geom=500;
    // init values
    std::for_each(g.begin(), g.end(), [&Vm](auto &gp) {gp.Vm=Vm;gp.Vu=0.;gp.H=gp.Cp*gp.Tt; });
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
    quiss::curvature_solver(gi);
    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = quiss::make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        // add_value(g, structuredGrid, "beta_metal_deg", [](const auto &gp) { return gp.bet * 180 / PI; });
        // structuredGrid->GetPointData()->SetActiveScalars("beta_metal_deg");
        quiss::plot_vtkStructuredGrid(structuredGrid, true);

        vtkNew<vtkXMLStructuredGridWriter> writer;
        writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_001_Vm.vts");
        // writer->SetFileName("C:/Users/sebastien/workspace/tbslib/tests/out/test_002_Vm.vts");
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
}