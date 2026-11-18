#include "yams_doctest.h"
#include <meridionalapi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>
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

namespace
{
    // Distance from (z, r) to the polyline of the curve sampled at n + 1 points.
    double distance_to(const Curve2d &c, double z, double r, size_t n = 2000)
    {
        auto [u1, u2] = c.bounds();
        double d = std::numeric_limits<double>::max();
        auto a = c(u1);
        for (size_t k = 1; k <= n; k++)
        {
            auto b = c(u1 + (u2 - u1) * k / n);
            const double ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
            const double t = l2 > 0. ? std::clamp(((z - a[0]) * ex + (r - a[1]) * ey) / l2, 0., 1.) : 0.;
            d = std::min(d, std::hypot(a[0] + t * ex - z, a[1] + t * ey - r));
            a = b;
        }
        return d;
    }

    // The rotor with curved edges: a bowed leading edge, a leaning trailing edge.
    yams::api::Problem annulus_rotor_traces()
    {
        auto pb = annulus_rotor();
        pb.rows[0].le_trace = Curve2d{{{0.38, 0.5}, {0.34, 0.75}, {0.42, 1.0}}, {0., 0., 0., 1., 1., 1.}, 2};
        pb.rows[0].te_trace = straight(0.58, 0.5, 0.64, 1.0);
        pb.options.n_stations = 41;
        return pb;
    }
}

TEST_CASE("meridionalapi.edge_traces")
{
    const auto pb = annulus_rotor_traces();
    std::vector<yams::api::RowStations> rows;
    {
        // The edges' stations follow their traces, their ends on the walls; the
        // stations in between blend them.
        auto sc = yams::api::make_case(pb, &rows);
        const auto &g = *sc.gi->g;
        const size_t nj = g.nCols();
        const auto [i_le, i_te] = rows[0];
        REQUIRE(g(i_le, 0).x == doctest::Approx(0.38));
        REQUIRE(g(i_le, 0).y == doctest::Approx(0.5));
        REQUIRE(g(i_le, nj - 1).x == doctest::Approx(0.42));
        REQUIRE(g(i_le, nj - 1).y == doctest::Approx(1.0));
        for (size_t j = 0; j < nj; j++)
        {
            REQUIRE(distance_to(*pb.rows[0].le_trace, g(i_le, j).x, g(i_le, j).y) < 1e-4);
            REQUIRE(distance_to(*pb.rows[0].te_trace, g(i_te, j).x, g(i_te, j).y) < 1e-4);
            for (size_t i = 1; i < g.nRows(); i++)
                REQUIRE(g(i, j).x > g(i - 1, j).x);
        }
        // the bow: the middle of the leading edge upstream of its ends
        REQUIRE(g(i_le, nj / 2).x < 0.38);
        // the inlet and outlet stay at constant z
        for (size_t j = 0; j < nj; j++)
        {
            REQUIRE(g(0, j).x == 0.);
            REQUIRE(g(g.nRows() - 1, j).x == 1.);
        }
    }

    // A trace given from shroud to hub: the same grid.
    {
        auto pb2 = pb;
        auto poles = pb2.rows[0].le_trace->poles();
        std::reverse(poles.begin(), poles.end());
        pb2.rows[0].le_trace = Curve2d{poles, {0., 0., 0., 1., 1., 1.}, 2};
        auto a = yams::api::make_case(pb);
        auto b = yams::api::make_case(pb2);
        const auto &ga = *a.gi->g, &gb = *b.gi->g;
        for (size_t i = 0; i < ga.nRows(); i++)
            for (size_t j = 0; j < ga.nCols(); j++)
            {
                REQUIRE(gb(i, j).x == doctest::Approx(ga(i, j).x));
                REQUIRE(gb(i, j).y == doctest::Approx(ga(i, j).y));
            }
    }

    // Solved: the stream lines move along the stations, which keep their shape; the
    // flow leaves at the metal angle, rothalpy conserved.
    {
        auto res = yams::api::solve(pb);
        REQUIRE(res.status == yams::api::Status::Converged);
        const auto &g = *res.solver_case.gi->g;
        const auto [i_le, i_te] = res.rows[0];
        for (size_t j = 0; j < g.nCols(); j++)
        {
            // as close as before solving: the stations keep their shapes
            REQUIRE(distance_to(*pb.rows[0].le_trace, g(i_le, j).x, g(i_le, j).y) < 1e-4);
            REQUIRE(distance_to(*pb.rows[0].te_trace, g(i_te, j).x, g(i_te, j).y) < 1e-6);
            REQUIRE(g(i_te, j).bet == doctest::Approx(-20. * std::numbers::pi / 180.));
            REQUIRE(g(i_te, j).I == doctest::Approx(g(i_le - 1, j).I).epsilon(1e-12));
        }
    }

    // Edges that cross: refused.
    {
        auto bad = pb;
        bad.rows[0].te_trace = straight(0.6, 0.5, 0.36, 1.0);
        REQUIRE_THROWS_AS(yams::api::make_case(bad), std::invalid_argument);
    }
}

TEST_CASE("meridionalapi.edge_traces_off_the_walls")
{
    // The same edges, their traces going beyond the walls or stopping short of them
    // (a tip gap): the stations are the same, on the part between hub and shroud.
    const auto pb = annulus_rotor_traces();
    const auto ref = yams::api::make_case(pb);
    const auto &g_ref = *ref.gi->g;
    // the trailing edge's line: z = 0.58 + 0.12 (r - 0.5)
    auto te = [](double r1, double r2) { return straight(0.58 + 0.12 * (r1 - 0.5), r1, 0.58 + 0.12 * (r2 - 0.5), r2); };
    for (const auto &trace : {te(0.4, 1.1), te(0.52, 0.97), te(1.1, 0.52)})
    {
        auto pb2 = pb;
        pb2.rows[0].te_trace = trace;
        const auto sc = yams::api::make_case(pb2);
        const auto &g = *sc.gi->g;
        REQUIRE(g.nRows() == g_ref.nRows());
        for (size_t i = 0; i < g.nRows(); i++)
            for (size_t j = 0; j < g.nCols(); j++)
            {
                REQUIRE(g(i, j).x == doctest::Approx(g_ref(i, j).x).epsilon(1e-9));
                REQUIRE(g(i, j).y == doctest::Approx(g_ref(i, j).y).epsilon(1e-9));
            }
    }

    // a trace that is not across the channel
    auto bad = pb;
    bad.rows[0].te_trace = straight(0.58, 1.2, 0.64, 1.5);
    REQUIRE_THROWS_AS(yams::api::make_case(bad), std::invalid_argument);
}

TEST_CASE("meridionalapi.edge_trace_and_constant_z")
{
    // A row with one edge along its trace and the other at constant z
    auto pb = annulus_rotor_traces();
    pb.rows[0].te_trace.reset();
    std::vector<yams::api::RowStations> rows;
    const auto sc = yams::api::make_case(pb, &rows);
    const auto &g = *sc.gi->g;
    for (size_t j = 0; j < g.nCols(); j++)
    {
        REQUIRE(distance_to(*pb.rows[0].le_trace, g(rows[0].i_le, j).x, g(rows[0].i_le, j).y) < 1e-4);
        REQUIRE(g(rows[0].i_te, j).x == 0.6);
    }
}

TEST_CASE("meridionalapi.stations_inside_the_channel")
{
    // A hub with a bump and a leading edge leaning back over it: the line between the
    // trace's ends goes through the hub downstream of the bump. The station is the
    // part of the trace from the shroud down to the hub, and every station stays
    // between the walls.
    auto pb = annulus_rotor();
    pb.hub = Curve2d{{{0., 0.5}, {0.5, 0.9}, {1., 0.5}}, {0., 0., 0., 1., 1., 1.}, 2};
    pb.rows[0].le_trace = straight(0.9, 0.572, 0.1, 1.0);
    pb.rows[0].z_te = 0.95;
    std::vector<yams::api::RowStations> rows;
    const auto sc = yams::api::make_case(pb, &rows);
    const auto &g = *sc.gi->g;
    const yams::api::detail::Wall hub{pb.hub, "hub"};
    REQUIRE(g(rows[0].i_le, 0).x < 0.8); // cut where the trace first meets the hub
    for (size_t i = 0; i < g.nRows(); i++)
        for (size_t j = 1; j < g.nCols(); j++)
        {
            REQUIRE(g(i, j).y > hub.r(g(i, j).x));
            REQUIRE(g(i, j).y <= 1.0);
        }
}

TEST_CASE("meridionalapi.station_spacing")
{
    // 21 stations on a channel of length 1 and span 0.5: 0.05 / 0.5
    auto pb = annulus_rotor();
    pb.options.n_stations = 21;
    pb.options.max_geom = 1;
    REQUIRE(yams::api::solve(pb).station_spacing == doctest::Approx(0.1));
}

TEST_CASE("parallel.par_for")
{
    // every call is made, and the failure of the lowest index is the one rethrown
    std::vector<int> calls(8, 0);
    auto f = [&calls](size_t i)
    {
        calls[i]++;
        if (i == 5 || i == 2)
            throw std::runtime_error(std::to_string(i));
    };
    REQUIRE_THROWS_WITH_AS(yams::detail::par_for(calls.size(), f), "2", std::runtime_error);
    REQUIRE(calls == std::vector<int>(8, 1));
    REQUIRE_NOTHROW(yams::detail::par_for(0, f));
}

TEST_CASE("meridionalapi.diverged")
{
    // Stations much closer than the span, with a fixed relaxation: the stream
    // lines become unstable and cross. The run ends as diverged, with the reason,
    // rather than by an exception from the solver's middle, and gives the grid of
    // its last complete iteration.
    auto pb = annulus_rotor();
    pb.options.n_stations = 121;
    pb.options.local_relaxation = false;
    const auto res = yams::api::solve(pb);
    REQUIRE(res.status == yams::api::Status::Diverged);
    REQUIRE(res.solver_case.diverged);
    REQUIRE_FALSE(res.solver_case.converged);
    REQUIRE(res.divergence.find("cross") != std::string::npos);
    REQUIRE(res.iterations > 0);
    REQUIRE(res.iterations < pb.options.max_geom);
    REQUIRE(res.residual_max.size() == res.iterations);
    for (double r : res.residual_max)
        REQUIRE(std::isfinite(r));

    const auto &g = *res.solver_case.gi->g;
    for (size_t i = 0; i < g.nRows(); i++)
        for (size_t j = 1; j < g.nCols(); j++)
        {
            // the stream lines in order from hub to shroud
            REQUIRE(g(i, j).y > g(i, j - 1).y);
            // the abscissas are those of these positions
            const double dl = std::hypot(g(i, j).x - g(i, j - 1).x, g(i, j).y - g(i, j - 1).y);
            REQUIRE(g(i, j).l - g(i, j - 1).l == doctest::Approx(dl).epsilon(1e-9));
        }

    // the same case solves again from scratch: the status is this run's
    auto pb2 = annulus_rotor();
    pb2.options.n_stations = 21;
    REQUIRE(yams::api::solve(pb2).status == yams::api::Status::Converged);
}

TEST_CASE("meridionalapi.local_relaxation")
{
    // Stations at 0.0125 of the span: a fixed relaxation of 0.05 diverges, the
    // relaxation at each station from its aspect ratio converges, close to the
    // solution of a coarser grid (within its discretisation error). With the
    // most stable curvature, for the fewest iterations.
    auto pb = annulus_rotor();
    pb.options.n_stations = 81;
    pb.options.curvature = yams::CurvatureMethod::DoubleDifferentiation;
    pb.options.local_relaxation = false;
    // unstable: the stream lines cross or keep oscillating, whichever comes first
    REQUIRE(yams::api::solve(pb).status != yams::api::Status::Converged);

    pb.options.local_relaxation = true;
    const auto fine = yams::api::solve(pb);
    REQUIRE(fine.status == yams::api::Status::Converged);
    pb.options.n_stations = 41;
    const auto coarse = yams::api::solve(pb);
    REQUIRE(coarse.status == yams::api::Status::Converged);
    // fewer iterations on the coarser grid: its stations need less relaxation
    REQUIRE(coarse.iterations < fine.iterations);

    // same flow at the outlet, on each stream line
    const auto &gf = *fine.solver_case.gi->g, &gc = *coarse.solver_case.gi->g;
    const size_t nj = gf.nCols();
    for (size_t j = 0; j < nj; j++)
    {
        const auto &pf = gf(gf.nRows() - 1, j), &pc = gc(gc.nRows() - 1, j);
        CHECK(pf.y == doctest::Approx(pc.y).epsilon(2e-3));
        CHECK(pf.Vm == doctest::Approx(pc.Vm).epsilon(1e-2));
        CHECK(pf.Tt == doctest::Approx(pc.Tt).epsilon(1e-4));
    }

    // the factor of a station: Wilkinson's optimum from its spacing and Mach number
    auto gi = *fine.solver_case.gi;
    const size_t i = gf.nRows() / 2;
    double dm = std::numeric_limits<double>::max(), Mm2 = 0.;
    for (size_t j = 0; j < nj; j++)
    {
        const auto &gp = gf(i, j);
        dm = std::min({dm, gp.m - gf(i - 1, j).m, gf(i + 1, j).m - gp.m});
        Mm2 += gp.Vm * gp.Vm / (gp.ga * gi.R * gp.Ts) / nj;
    }
    const double A = gf(i, nj - 1).l / dm;
    // the curvature method's constant, by default
    const double B = yams::api::default_relax_constant(yams::CurvatureMethod::DoubleDifferentiation);
    REQUIRE(gi.RF_B == B);
    REQUIRE(Mm2 < gi.RF_Mm2_max);
    CHECK(yams::station_relaxation(gi, i) == doctest::Approx(1. / (1. + (1. - Mm2) * A * A / B)));
    // the Mach number limited, a fixed factor when not local
    gi.RF_Mm2_max = Mm2 / 2;
    CHECK(yams::station_relaxation(gi, i) == doctest::Approx(1. / (1. + (1. - Mm2 / 2) * A * A / B)));
    // no Mach number effect for a constant density
    gi.rho_cst = true;
    CHECK(yams::station_relaxation(gi, i) == doctest::Approx(1. / (1. + A * A / B)));
    gi.rho_cst = false;
    // a constant or a Mach limit that make no factor are refused
    gi.RF_B = 0.;
    CHECK_THROWS_AS(yams::station_relaxation(gi, i), std::invalid_argument);
    gi.RF_B = B;
    gi.RF_Mm2_max = 1.;
    CHECK_THROWS_AS(yams::station_relaxation(gi, i), std::invalid_argument);
    gi.RF_Mm2_max = 0.6;
    pb.options.relax_constant = -1.;
    CHECK_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);
    pb.options.relax_constant.reset();
    pb.options.relax_mach2_max = 1.;
    CHECK_THROWS_AS(yams::api::make_case(pb), std::invalid_argument);
    // two stations that touch on a stream line: no factor, said as a divergence
    // rather than the fixed factor (on a copy of the grid)
    gi.g = std::make_shared<yams::MeridionalGrid<double>>(gf);
    (*gi.g)(i, 0).m = (*gi.g)(i - 1, 0).m;
    CHECK_THROWS_AS(yams::station_relaxation(gi, i), yams::SolverDivergence);
    gi.RF_local = false;
    CHECK(yams::station_relaxation(gi, i) == gi.RF);
}

TEST_CASE("meridionalapi.local_relaxation_incompressible")
{
    // a constant density at a high velocity converges on the fine grid (the
    // factor then has no Mach number term: checked in local_relaxation)
    auto pb = annulus_rotor();
    pb.options.n_stations = 81;
    pb.options.curvature = yams::CurvatureMethod::DoubleDifferentiation;
    pb.options.compressible = false;
    pb.inlet.vm_mean = 200.;
    const auto res = yams::api::solve(pb);
    REQUIRE(res.status == yams::api::Status::Converged);
}

namespace
{
    // The largest zigzag of the stream lines' slope from one station to the next:
    // half the jump of its change where that change alternates in sign.
    double zigzag(const yams::MeridionalGrid<double> &g)
    {
        double worst = 0.;
        for (size_t j = 1; j + 1 < g.nCols(); j++)
            for (size_t i = 1; i + 2 < g.nRows(); i++)
            {
                auto slope = [&](size_t k) { return (g(k + 1, j).y - g(k, j).y) / (g(k + 1, j).x - g(k, j).x); };
                const double d1 = slope(i) - slope(i - 1), d2 = slope(i + 1) - slope(i);
                if (d1 * d2 < 0.)
                    worst = std::max(worst, std::abs(d1 - d2) / 2);
            }
        return worst;
    }
}

TEST_CASE("meridionalapi.curvature_methods")
{
    using yams::CurvatureMethod;
    // On a circle of radius R, all methods give 1 / R, positive turning
    // counterclockwise in (z, r); on a zigzag from one station to the next, the
    // double differentiation sees no curvature at all, the others do.
    const double R = 2.;
    const size_t ni = 21, nj = 3;
    std::vector<std::array<double, 2>> arc, zig;
    for (size_t j = 0; j < nj; j++)
        for (size_t i = 0; i < ni; i++)
        {
            const double t = 0.5 * static_cast<double>(i) / (ni - 1), rr = R + 0.1 * j;
            arc.push_back({rr * std::sin(t), 3. - rr * std::cos(t)});
            zig.push_back({0.05 * i, 1. + 0.1 * j + (i % 2 ? 0.002 : -0.002)});
        }
    for (auto m : {CurvatureMethod::DoubleDifferentiation, CurvatureMethod::Blended, CurvatureMethod::ThreePoint})
    {
        auto g = yams::make_meridional_grid<double>(arc, ni, nj);
        auto g_metrics = yams::Grid2dMetrics<double>{ni, nj};
        yams::compute_grid_metrics(g, g_metrics, yams::fm, yams::fl, m);
        for (size_t i = 1; i + 1 < ni; i++)
            CHECK(g(i, 0).cur == doctest::Approx(1. / R).epsilon(1e-3));

        auto z = yams::make_meridional_grid<double>(zig, ni, nj);
        yams::compute_grid_metrics(z, g_metrics, yams::fm, yams::fl, m);
        const double k = std::abs(z(ni / 2, 1).cur);
        if (m == CurvatureMethod::DoubleDifferentiation)
            CHECK(k < 1e-9);
        else
            CHECK(k > 1.);
    }
    // the less stable the method, the smaller its relaxation constant
    CHECK(yams::api::default_relax_constant(CurvatureMethod::Blended) <
          yams::api::default_relax_constant(CurvatureMethod::DoubleDifferentiation));
    CHECK(yams::api::default_relax_constant(CurvatureMethod::ThreePoint) <
          yams::api::default_relax_constant(CurvatureMethod::Blended));
}

TEST_CASE("meridionalapi.curvature_waves")
{
    // The rotor's leading edge leaves a zigzag on the stream lines with the
    // double differentiation; the default (Blended) curvature damps it, for the
    // same flow, the case taking its method's relaxation constant.
    auto pb = annulus_rotor_traces();
    pb.options.curvature = yams::CurvatureMethod::DoubleDifferentiation;
    const auto dd = yams::api::solve(pb);
    REQUIRE(dd.status == yams::api::Status::Converged);
    pb.options.curvature = yams::api::Options{}.curvature;
    REQUIRE(pb.options.curvature == yams::CurvatureMethod::Blended);
    const auto bl = yams::api::solve(pb);
    REQUIRE(bl.status == yams::api::Status::Converged);
    CHECK(bl.solver_case.gi->curvature == yams::CurvatureMethod::Blended);
    CHECK(bl.solver_case.gi->RF_B == yams::api::default_relax_constant(yams::CurvatureMethod::Blended));

    const auto &gd = *dd.solver_case.gi->g, &gb = *bl.solver_case.gi->g;
    std::printf("zigzag: double differentiation %g, blended %g\n", zigzag(gd), zigzag(gb));
    CHECK(zigzag(gb) < zigzag(gd) / 3);
    for (size_t j = 0; j < gd.nCols(); j++)
    {
        const auto &pd = gd(gd.nRows() - 1, j), &pb_ = gb(gb.nRows() - 1, j);
        CHECK(pb_.Tt == doctest::Approx(pd.Tt).epsilon(1e-4));
        CHECK(pb_.Vm == doctest::Approx(pd.Vm).epsilon(1e-2));
    }
}
