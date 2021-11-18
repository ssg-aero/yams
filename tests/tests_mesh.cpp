#include <gtest/gtest.h>
#include <meshtools.h>
#include <gbs-render/vtkcurvesrender.h>
#include <gbs-render/vtkgridrender.h>
using T = double;
const T tol{static_cast<T>(1e-6)};
const bool PLOT_ON = true;

TEST(tests_mesh, channel)
{

    std::vector<std::array<T,2>> poles1{
            {0.0,0.0},
            {0.2,0.0},
            {0.8,0.3},
            {1.0,0.3},
        };
    std::vector<std::array<T,2>> poles2{
            {0.0,0.5},
            {0.4,0.5},
            {0.8,1.0},
            {1.0,1.0},
        };
    std::vector<T> knots{0.,0.5,1.};
    std::vector<size_t> mult{3,1,3};
    size_t deg{2};

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

    yams::crv_vector<T> crv_lst{crv1, crv2};
    auto [ crv_lst_opp , v] = yams::build_iso_ksi_curves<T>(crv_lst,knots);

    ASSERT_LT(
        gbs::distance(
            crv1->begin(),
            crv_lst_opp.front()->begin()),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv1->end(),
            crv_lst_opp.back()->begin()),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv2->begin(),
            crv_lst_opp.front()->end()),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv2->end(),
            crv_lst_opp.back()->end()),
        tol);

    size_t nu{30};
    size_t nv{15};

    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = yams::mesh_channel<T>(crv_lst, knots, nu, nv);

    ASSERT_LT(
        gbs::distance(
            crv1->begin(),
            pts[0]),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv2->begin(),
            pts[ni - 1]),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv1->end(),
            pts[ni * (nj - 1)]),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv2->end(),
            pts[ni-1 + ni * (nj - 1)]),
        tol);

    if(PLOT_ON)
    {
        auto sgrid_actor = gbs::make_structuredgrid_actor(pts, ni, nj);
        gbs::plot(crv_lst, crv_lst_opp, sgrid_actor);
    }
}