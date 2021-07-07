#include <gtest/gtest.h>
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>
#include <gridrender.h>

const double PI = acos(-1.);

using namespace quiss;

TEST(tests_visu, render_curvature)
{
    size_t ni = 50;
    size_t nj = 15;
    MeridionalGrid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);

    auto structuredGrid = make_vtkStructuredGrid(g);
    add_value(g,structuredGrid,"Curvature",[](const auto &gp){return gp.cur;});
    plot_vtkStructuredGrid(structuredGrid,"Curvature");

}