#include "yams_doctest.h"
#include "rotor37_case.h"
#include <streamtube.h>

#include <cmath>
#include <numbers>

namespace
{
    constexpr double deg = std::numbers::pi / 180.;

    // NASA Rotor 37 at 55 % of its speed and of its design flow: the same flow
    // angles, a subsonic relative flow (the blade-to-blade solver is subsonic).
    yams::api::Problem rotor37_part_speed()
    {
        auto pb = yams::tests::rotor37(0.55 * 150.);
        pb.rows[0].rpm *= 0.55;
        return pb;
    }

    // A section on the tube from the row's metal angle at relative span l: camber
    // dtheta/dm = tan(kappa) / r, parabolic thickness of t_c times the chord,
    // normal to the camber line.
    yams::api::BladeSection section(const yams::api::Problem &pb, const yams::api::StreamTube &tube, double l,
                                    double t_c = 0.06, size_t n = 41)
    {
        const auto &rw = tube.rows[0];
        const auto c_m = rw.m_te - rw.m_le;
        std::vector<double> m(n), th(n), kap(n);
        for (size_t k = 0; k < n; k++)
        {
            m[k] = rw.m_le + c_m * k / (n - 1);
            kap[k] = pb.rows[0].metal_angle(double(k) / (n - 1), l);
        }
        auto slope = [&](size_t q) { return std::tan(kap[q]) / tube.at(m[q]).r; };
        th[0] = 0.;
        double chord = 0.;
        for (size_t k = 1; k < n; k++)
        {
            th[k] = th[k - 1] + 0.5 * (slope(k - 1) + slope(k)) * (m[k] - m[k - 1]);
            chord += (m[k] - m[k - 1]) / std::cos(0.5 * (kap[k - 1] + kap[k]));
        }
        gbs::points_vector<double, 2> s1, s2;
        for (size_t k = 0; k < n; k++)
        {
            const auto x = double(k) / (n - 1);
            const auto h = 0.5 * t_c * chord * 4. * x * (1. - x);
            const auto dth = h / (tube.at(m[k]).r * std::cos(kap[k]));
            s1.push_back({m[k], th[k] + dth});
            s2.push_back({m[k], th[k] - dth});
        }
        return {gbs::interpolate(s1, m, 3), gbs::interpolate(s2, m, 3)};
    }

    struct Rotor37Tube
    {
        yams::api::Problem pb = rotor37_part_speed();
        yams::api::Result res = yams::api::solve(pb);
        size_t j{};
        yams::api::StreamTube tube;
        double l{}; // relative span of the stream line at mid chord
        Rotor37Tube()
        {
            REQUIRE(res.status == yams::api::Status::Converged);
            const auto &g = *res.solver_case.gi->g;
            j = g.nCols() / 2;
            tube = yams::api::stream_tube(pb, res, j);
            const auto i_mid = (res.rows[0].i_le + res.rows[0].i_te) / 2;
            l = g(i_mid, j).l / g(i_mid, g.nCols() - 1).l;
        }
    };
}

TEST_CASE("streamtube.rotor37")
{
    const Rotor37Tube r;
    const auto &g = *r.res.solver_case.gi->g;
    const auto &tube = r.tube;
    const size_t n_sl = g.nCols(), n_st = g.nRows(), j = r.j;
    REQUIRE(tube.stations.size() == n_st);
    // halfway in mass flow to the neighbours: a tenth of the flow with 11 stream
    // lines (within the meridional solver's tolerance on the stream lines' flow)
    CHECK(tube.mass_flow == doctest::Approx(g(0, n_sl - 1).q / (n_sl - 1)).epsilon(1e-3));
    for (size_t i = 0; i < n_st; i++)
    {
        const auto &st = tube.stations[i];
        CHECK(st.b > 0.);
        // the bounding surfaces lie between the stream line and its neighbours
        CHECK(tube.lower(st.m)[1] > g(i, j - 1).y);
        CHECK(tube.lower(st.m)[1] < st.r);
        CHECK(tube.upper(st.m)[1] > st.r);
        CHECK(tube.upper(st.m)[1] < g(i, j + 1).y);
        // the mean curve passes through the stream line's nodes, parameter m
        CHECK(tube.mean(st.m)[0] == doctest::Approx(st.z));
        CHECK(tube.mean(st.m)[1] == doctest::Approx(st.r));
    }
    const auto mid = tube.at(0.5 * (tube.stations[3].m + tube.stations[4].m));
    CHECK(mid.r == doctest::Approx(0.5 * (tube.stations[3].r + tube.stations[4].r)));
    REQUIRE(tube.rows.size() == 1);
    CHECK(tube.rows[0].omega == doctest::Approx(0.55 * 17188.7 * std::numbers::pi / 30.));
    CHECK(tube.rows[0].blade_count == 36);
    CHECK(tube.rows[0].m_le == g(r.res.rows[0].i_le, j).m);
    CHECK(tube.rows[0].m_te == g(r.res.rows[0].i_te, j).m);

    CHECK_THROWS_AS(yams::api::stream_tube(r.pb, r.res, n_sl), std::invalid_argument);
}

TEST_CASE("streamtube.b2b_rotor37")
{
    const Rotor37Tube r;
    const auto &tube = r.tube;
    const auto sec = section(r.pb, tube, r.l);

    const auto c = yams::api::make_b2b_case(tube, 0, sec);
    CHECK(c.points.size() == 11 * c.n_planes);
    CHECK(c.j_le > 1);
    CHECK(c.j_te > c.j_le + 1);
    CHECK(c.j_te + 2 < c.n_planes);
    CHECK(c.m[c.j_le] == doctest::Approx(tube.rows[0].m_le));
    CHECK(c.m[c.j_te] == doctest::Approx(tube.rows[0].m_te));
    // relative inlet flow of a rotor without inlet swirl
    const auto in = tube.at(c.m.front());
    CHECK(c.beta_in == doctest::Approx(std::atan2(-c.omega * in.r, in.Vm)));
    for (auto dr : c.dr)
        CHECK(dr > 0.);

    auto solver = yams::api::make_b2b_solver(c);
    solver.computeW();
    REQUIRE(solver.converged());
    // 696 iterations with clang on arm64
    CHECK(solver.iterations() < 1000);
    CHECK(solver.fullBladePassageMassFlow() == doctest::Approx(tube.mass_flow));

    const auto &dat = solver.data();
    const auto &msh = solver.mesh();
    const auto [ni, nj] = solver.dimensions();
    CHECK(*std::ranges::max_element(dat.MW) < 1.);
    // the inlet keeps its flow angle on the stagnation line, within the solver's
    // relocation of the passage's inner nodes
    CHECK(msh.BT[0] == doctest::Approx(c.beta_in).epsilon(0.01));
    CHECK(msh.BT[ni / 2] == doctest::Approx(c.beta_in).epsilon(0.05));
    // Outlet: the meridional run imposes the metal angle at the trailing edge
    // (DIRECT mode); the blade-to-blade flow leaves it with a deviation and does
    // less work (-38.0 deg against -34.8, Tt rise 25.3 K against 29.5 with clang on arm64)
    const auto out = tube.at(c.m.back());
    const auto beta_merid = std::atan2(out.Vu - c.omega * out.r, out.Vm);
    const auto id_out = ni / 2 + ni * (nj - 1);
    CHECK(msh.BT[id_out] < beta_merid - 1. * deg);
    CHECK(msh.BT[id_out] > beta_merid - 8. * deg);
    CHECK(dat.TT[id_out] < out.Tt);
    CHECK(dat.TT[id_out] > out.Tt - 0.25 * (out.Tt - in.Tt));
}

TEST_CASE("streamtube.b2b_errors")
{
    const Rotor37Tube r;
    const auto &tube = r.tube;
    const auto sec = section(r.pb, tube, r.l);
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 1, sec), std::invalid_argument);
    yams::api::B2BOptions opt;
    opt.plane_spacing = 0.;
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, sec, opt), std::invalid_argument);
    opt = {};
    opt.n_pitch = 2;
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, sec, opt), std::invalid_argument);
    // a side running backwards in m
    auto back = sec;
    back.side1.reverse();
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, back), std::invalid_argument);
    // a section beyond the tube's outlet
    const auto shift = tube.stations.back().m;
    gbs::points_vector<double, 2> p1{{shift + 0.01, 0.}, {shift + 0.02, 0.01}}, p2{{shift + 0.01, 0.}, {shift + 0.02, 0.}};
    yams::api::BladeSection outside{gbs::interpolate(p1, std::vector<double>{0., 1.}, 1),
                                    gbs::interpolate(p2, std::vector<double>{0., 1.}, 1)};
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, outside), std::invalid_argument);
}
