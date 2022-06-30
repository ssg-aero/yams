#include <yams/meshtools.h>

#include <gtest/gtest.h>

#include <gbs-render/vtkcurvesrender.h>
#include <gbs-render/vtkgridrender.h>

using T = double;
const T tol{static_cast<T>(1e-6)};
const bool PLOT_ON = true;

TEST(tests_mesh, channel)
{

    std::vector<std::array<T,2>> poles1{
            {0.0,0.0},
            {0.1,0.0},
            {0.5,0.2},
            {0.8,0.3},
            {1.0,0.3},
        };
    std::vector<std::array<T,2>> poles2{
            {0.0,0.5},
            {0.1,0.5},
            {0.5,0.6},
            {0.8,1.0},
            {1.0,1.0},
        };
    std::vector<T> knots{0.,0.33,0.66,1.};
    std::vector<size_t> mult{3,1,1,3};
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

    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = yams::mesh_channel<T>(crv_lst, knots, nv, nu);

    ASSERT_LT(
        gbs::distance(
            crv1->begin(),
            pts[0]),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv1->end(),
            pts[ni - 1]),
        tol);
    ASSERT_LT(
        gbs::distance(
            crv2->begin(),
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
        gbs::plot(
            crv_lst, 
            crv_lst_opp, 
            sgrid_actor
        );
    }
}

TEST(tests_mesh, channel_arbitrary)
{
    size_t p = 3;
    using T = double;
    using namespace gbs;
    // Set 2 curves with diferent parametrization
    points_vector<T,2> crv1_pts{
        {0.00045, 0.00000}, // i0
        {0.00601, 0.00000},
        {0.02664, 0.01632},
        {0.03090, 0.02746},
        {0.03591, 0.03414},
        {0.03900, 0.03600}, //i1
        {0.05082, 0.04138},
        {0.06000, 0.04550}, //i2
        {0.06557, 0.04677},
        {0.07500, 0.04700}, //i3
        {0.07999, 0.04700},
        {0.10500, 0.04700}, //i4
        {0.11014, 0.04699},
        {0.13825, 0.04493},
        {0.16003, 0.03984}  // i5
    };
    auto crv1 = interpolate(
        crv1_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );

    points_vector<T,2> crv2_pts{
        {0.00010, 0.08550}, // i0
        {0.03900, 0.08550}, // i1
        {0.06000, 0.08550}, // i2
        {0.08000, 0.08550}, // i3
        {0.11000, 0.08550}, // i4
        {0.16000, 0.08550}  // i5
    };
    auto crv2 = interpolate(
        crv2_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );
    // set hard points
    std::vector<size_t> i_crv_1 {0, 5, 7, 9, 11, crv1_pts.size() - 1 };
    std::vector<size_t> i_crv_2 {0, 1, 2, 3, 4 , crv2_pts.size() - 1};
    std::vector<points_vector<T,2>> hard_points(2);
    auto n_hard_pts = i_crv_1.size();
    for(size_t i{}; i < n_hard_pts; i++)
    {
        hard_points[0].push_back(crv1_pts[i_crv_1[i]]);
        hard_points[1].push_back(crv2_pts[i_crv_2[i]]);
    }
    //
    size_t nu{120};
    size_t nv{21};
    // size_t nu{40};
    // size_t nv{7};
    auto [pts, ni, nj, n_iso_eth, n_iso_ksi] = yams::mesh_channel<T>(
        {
            std::make_shared<BSCurve<T,2>>(crv1),
            std::make_shared<BSCurve<T,2>>(crv2)
        },
        hard_points,
        nv,
        nu,
        1e-6
    );

    yams::smooth_mesh(pts,nj, n_iso_eth);

    if(PLOT_ON)
    {
        auto sgrid_actor = gbs::make_structuredgrid_actor(pts, ni, nj);
        gbs::plot(
            sgrid_actor,
            crv1,
            crv2
        );
    }
}