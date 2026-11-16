#include "yams_doctest.h"
#include <meridionalapi.h>

#include <numbers>

namespace
{
    using yams::api::Curve2d;

    Curve2d straight(double z1, double r1, double z2, double r2)
    {
        return Curve2d{{{z1, r1}, {z2, r2}}, {0., 0., 1., 1.}, 1};
    }

    // Annulus with one rotor in DIRECT mode
    yams::api::Problem annulus_rotor()
    {
        constexpr double deg = std::numbers::pi / 180.;
        yams::api::Problem pb;
        pb.hub = straight(0., 0.5, 1., 0.5);
        pb.shroud = straight(0., 1.0, 1., 1.0);
        pb.rows.push_back({
            .name = "R1",
            .z_le = 0.4,
            .z_te = 0.6,
            .rpm = 500.,
            .blade_count = 20,
            .metal_angle = [deg](double m, double) { return ((1 - m) * -50. + m * -20.) * deg; },
        });
        pb.inlet.vm_mean = 30.;
        pb.inlet.Tt = [](double l) { return 288.15 + 10. * l; };
        pb.options.n_stations = 41;
        return pb;
    }
}

TEST_CASE("meridionalapi.rotor")
{
    auto pb = annulus_rotor();
    size_t calls{};
    auto res = yams::api::solve(pb, [&calls](size_t, double, double) { calls++; return true; });

    REQUIRE(res.status == yams::api::Status::Converged);
    REQUIRE(res.iterations > 0);
    REQUIRE(calls == res.iterations);
    REQUIRE(res.rows.size() == 1);

    const auto &g = *res.solver_case.gi->g;
    const auto [i_le, i_te] = res.rows[0];
    const size_t nj = g.nCols();
    REQUIRE(nj == pb.options.n_streams);
    for (size_t j = 0; j < nj; j++)
    {
        // edges on computing stations
        REQUIRE(g(i_le, j).x == doctest::Approx(0.4));
        REQUIRE(g(i_te, j).x == doctest::Approx(0.6));
        // flow leaves at the metal angle, rothalpy conserved, positive work
        REQUIRE(g(i_te, j).bet == doctest::Approx(-20. * std::numbers::pi / 180.));
        REQUIRE(g(i_te, j).I == doctest::Approx(g(i_le - 1, j).I).epsilon(1e-12));
        REQUIRE(g(i_te, j).H > g(i_le - 1, j).H);
    }
    // the inlet values are the totals, given at the initial (evenly spread) stream lines
    for (size_t j = 0; j < nj; j++)
    {
        const double l = double(j) / (nj - 1);
        REQUIRE(g(0, j).Tt == doctest::Approx(288.15 + 10. * l));
        REQUIRE(g(0, j).Pt == doctest::Approx(1.01325e5));
    }
}

TEST_CASE("meridionalapi.cancel")
{
    auto res = yams::api::solve(annulus_rotor(), [](size_t it, double, double) { return it < 3; });
    REQUIRE(res.status == yams::api::Status::Cancelled);
    REQUIRE(res.iterations == 3);
}

TEST_CASE("meridionalapi.invalid_problems")
{
    auto pb = annulus_rotor();
    pb.rows[0].z_te = 0.3; // trailing edge before leading edge
    REQUIRE_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);

    pb = annulus_rotor();
    pb.rows[0].z_le = -0.1; // outside the channel
    REQUIRE_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);

    pb = annulus_rotor();
    pb.rows[0].metal_angle = nullptr;
    REQUIRE_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);

    pb = annulus_rotor();
    pb.hub = straight(1., 0.5, 0., 0.5); // z decreasing along the hub
    REQUIRE_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);
}

TEST_CASE("meridionalapi.stations_on_hard_positions")
{
    // Bounds for which z_in + (z_out - z_in) * n / n rounds one ulp above z_out:
    // the outlet station must be z_out itself, not outside the channel.
    yams::api::Problem pb;
    pb.hub = straight(0.12, 0.5, 1.43, 0.5);
    pb.shroud = straight(0.12, 1.0, 1.43, 1.0);
    pb.options.n_stations = 21;
    auto sc = yams::api::make_case(pb);
    const auto &g = *sc.gi->g;
    REQUIRE(g.nRows() == 21);
    REQUIRE(g(0, 0).x == 0.12);
    REQUIRE(g(g.nRows() - 1, 0).x == 1.43);

    // row edges are exactly on their stations
    std::vector<yams::api::RowStations> rows;
    sc = yams::api::make_case(annulus_rotor(), &rows);
    REQUIRE((*sc.gi->g)(rows[0].i_le, 0).x == 0.4);
    REQUIRE((*sc.gi->g)(rows[0].i_te, 0).x == 0.6);
}

TEST_CASE("meridionalapi.empty_grid")
{
    // an empty grid is an explicit error, not a division by zero
    REQUIRE_THROWS_AS(yams::make_grid_info<double>(yams::MeridionalGrid<double>{}), std::invalid_argument);
    REQUIRE_THROWS_AS(
        yams::make_solver_case<double>(yams::MeridionalGrid<double>{}, std::vector<yams::BladeInfo<double>>{}),
        std::invalid_argument);
}
