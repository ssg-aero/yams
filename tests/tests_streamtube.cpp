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
    CHECK(tube.rows[0].m_le == tube.stations[r.res.rows[0].i_le].m);
    CHECK(tube.rows[0].m_te == tube.stations[r.res.rows[0].i_te].m);
    // m is the mean curve's arc length, the abscissa the blade-to-blade solver
    // computes from the parameter, close to the meridional grid's polyline length
    for (size_t i = 1; i < n_st; i++)
    {
        const auto &st = tube.stations[i];
        CHECK(gbs::length(tube.mean, tube.stations.front().m, st.m) == doctest::Approx(st.m - tube.stations.front().m).epsilon(1e-6));
        CHECK(st.m == doctest::Approx(g(i, j).m).epsilon(1e-3));
    }

    // the tubes at the hub and at the shroud, and the channel's flow shared
    double flow = 0.;
    for (size_t k = 0; k < n_sl; k++)
    {
        const auto t = yams::api::stream_tube(r.pb, r.res, k);
        flow += t.mass_flow;
        for (const auto &st : t.stations)
            CHECK(st.b > 0.);
    }
    CHECK(flow == doctest::Approx(g(0, n_sl - 1).q));
    const auto hub = yams::api::stream_tube(r.pb, r.res, 0);
    CHECK(hub.lower(hub.stations[2].m)[1] == doctest::Approx(g(2, 0).y));
    const auto shroud = yams::api::stream_tube(r.pb, r.res, n_sl - 1);
    CHECK(shroud.upper(shroud.stations[2].m)[1] == doctest::Approx(g(2, n_sl - 1).y));

    CHECK_THROWS_AS(yams::api::stream_tube(r.pb, r.res, n_sl), std::invalid_argument);
}

TEST_CASE("streamtube.invalid_flow")
{
    // a choked meridional run is not a flow to analyse along its stream lines
    auto pb = yams::tests::rotor37(180.);
    const auto res = yams::api::solve(pb);
    REQUIRE(res.status == yams::api::Status::Choked);
    CHECK_THROWS_AS(yams::api::stream_tube(pb, res, 5), std::invalid_argument);
}

TEST_CASE("streamtube.b2b_rotor37")
{
    const Rotor37Tube r;
    const auto &tube = r.tube;
    const auto sec = section(r.pb, tube, r.l);

    // Blended curvature (the API's default is the double differentiation)
    yams::api::B2BOptions blended;
    blended.curvature = yams::CurvatureMethod::Blended;
    const auto c = yams::api::make_b2b_case(tube, 0, sec, blended);
    CHECK(c.points.size() == 11 * c.n_planes);
    CHECK(c.j_le > 1);
    CHECK(c.j_te > c.j_le + 1);
    CHECK(c.j_te + 2 < c.n_planes);
    CHECK(c.m[c.j_le] == doctest::Approx(tube.rows[0].m_le));
    CHECK(c.m[c.j_te] == doctest::Approx(tube.rows[0].m_te));
    // the planes are the meridional stations (Novak and Hearsey 1977)
    REQUIRE(c.station.size() == c.n_planes);
    for (size_t k = 0; k < c.n_planes; k++)
    {
        if (k != c.j_le && k != c.j_te)
            CHECK(c.m[k] == tube.stations[c.station[k]].m);
        if (k > 0)
            CHECK(c.station[k] == c.station[k - 1] + 1);
    }
    CHECK(c.station[c.j_le] == r.res.rows[0].i_le);
    CHECK(c.station[c.j_te] == r.res.rows[0].i_te);
    // about 1.3 with 31 stations and 11 nodes across, above the solver's limit
    CHECK(c.spacing_ratio > 1.1);
    // relative inlet flow of a rotor without inlet swirl
    const auto in = tube.at(c.m.front());
    CHECK(c.beta_in == doctest::Approx(std::atan2(-c.omega * in.r, in.Vm)));
    for (auto dr : c.dr)
        CHECK(dr > 0.);

    auto solver = yams::api::make_b2b_solver(c);
    solver.computeW();
    REQUIRE(solver.converged());
    // 498 iterations with clang on arm64 (Blended curvature)
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
    // Outlet (Novak and Hearsey 1977): the far field angle is the meridional run's, the
    // metal angle in DIRECT mode, no deviation; the stagnation lines keep it. With no
    // deviation the Kutta condition fails at the trailing edge, and the outlet flow is
    // not uniform across the pitch: the inner stream lines turn away from the imposed
    // angle (kutta -0.25 with clang on arm64).
    const auto out = tube.at(c.m.back());
    CHECK(c.beta_out == doctest::Approx(std::atan2(out.Vu - c.omega * out.r, out.Vm)));
    CHECK(msh.BT[ni * (nj - 1)] == doctest::Approx(c.beta_out).epsilon(0.02));
    const auto o = yams::api::b2b_output(c, solver);
    CHECK(std::abs(o.kutta) > 0.1);
    REQUIRE(o.B.size() == nj);
    REQUIRE(o.beta_mean.size() == nj);
    CHECK(o.beta_out == msh.BT[ni / 2 + ni * (nj - 1)]);
    // B is 1 far upstream and dips inside the blade (Novak and Hearsey's Fig. 4)
    CHECK(std::abs(o.B.front() - 1.) < 0.05);
    CHECK(*std::ranges::min_element(o.B.begin() + c.j_le, o.B.begin() + c.j_te + 1) < 0.95);

    // The far field angle that closes the Kutta condition: 3.8 deg of deviation past
    // the metal angle (offset -3.80 deg in 3 solutions with clang on arm64), and an
    // outlet flow uniform across the pitch within 1.2 deg
    const auto k = yams::api::b2b_kutta(tube, 0, sec, blended);
    REQUIRE(k.converged);
    CHECK(std::abs(k.output.kutta) <= 1e-3);
    CHECK(k.runs <= 8);
    // stopped before closing it: not converged, the best of its solutions returned
    {
        const auto k2 = yams::api::b2b_kutta(tube, 0, sec, blended, 1e-3, 2);
        CHECK_FALSE(k2.converged);
        CHECK(k2.runs == 2);
        CHECK(k2.solver.converged());
        CHECK(std::abs(k2.output.kutta) < std::abs(o.kutta));
    }
    CHECK(k.outlet_offset < -1. * deg);
    CHECK(k.outlet_offset > -8. * deg);
    CHECK(k.case_.beta_out == doctest::Approx(c.beta_out + k.outlet_offset));
    {
        const auto [nik, njk] = k.solver.dimensions();
        const auto &mk = k.solver.mesh();
        const auto [lo, hi] = std::ranges::minmax(std::vector<double>(mk.BT.end() - nik, mk.BT.end()));
        CHECK(hi - lo < 1.5 * deg);
        CHECK(njk == nj);
    }
}

TEST_CASE("streamtube.b2b_errors")
{
    const Rotor37Tube r;
    const auto &tube = r.tube;
    const auto sec = section(r.pb, tube, r.l);
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 1, sec), std::invalid_argument);
    yams::api::B2BOptions opt;
    // a section whose edges are off the row's edge stations by half their spacing
    {
        auto shifted = tube;
        const auto i_le = r.res.rows[0].i_le;
        const auto dm = 0.5 * (tube.stations[i_le + 1].m - tube.stations[i_le].m);
        shifted.rows[0].m_le += dm;
        shifted.rows[0].m_te += dm;
        CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, section(r.pb, shifted, r.l)), std::invalid_argument);
    }
    opt.n_pitch = 2;
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, sec, opt), std::invalid_argument);
    // planes closer than the nodes across: 5 nodes give a spacing ratio about 0.5
    opt.n_pitch = 5;
    CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, sec, opt), std::invalid_argument);
    opt = {};
    // a neighbouring row one station downstream of the trailing edge: no two
    // stations before halfway to it, the domain does not reach into its blade
    {
        auto two_rows = tube;
        const auto i_te = r.res.rows[0].i_te;
        two_rows.rows.push_back({"Stator", tube.stations[i_te + 1].m, tube.stations[i_te + 5].m, 0., 46});
        CHECK_THROWS_AS(yams::api::make_b2b_case(two_rows, 0, sec), std::invalid_argument);
        // three stations further, it leaves room
        two_rows.rows.back().m_le = tube.stations[i_te + 5].m;
        two_rows.rows.back().m_te = tube.stations[i_te + 8].m;
        const auto c2 = yams::api::make_b2b_case(two_rows, 0, sec);
        CHECK(c2.m.back() <= 0.5 * (tube.rows[0].m_te + two_rows.rows.back().m_le));
    }
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
    // a blunt trailing edge: the grid takes pointed edges
    {
        const auto &s2 = sec.side2;
        auto [u1, u2] = s2.bounds();
        gbs::points_vector<double, 2> pts;
        std::vector<double> u;
        for (size_t k = 0; k < 21; k++)
        {
            const auto x = double(k) / 20;
            auto p = s2(u1 + (u2 - u1) * x);
            p[1] -= 0.02 * x; // opens up to 0.02 rad at the trailing edge
            pts.push_back(p);
            u.push_back(p[0]);
        }
        yams::api::BladeSection blunt{sec.side1, gbs::interpolate(pts, u, 3)};
        CHECK_THROWS_AS(yams::api::make_b2b_case(tube, 0, blunt), std::invalid_argument);
    }
}

TEST_CASE("streamtube.b2b_curvature")
{
    // The blade-to-blade stream lines' curvature (CurvatureMethod, as the meridional
    // solver's): the double differentiation leaves a zigzag of the stream lines from
    // one plane to the next, which Blended divides by 5 and ThreePoint by 15 (0.154,
    // 0.032, 0.010 with clang on arm64, theta's alternating part against its variation)
    const Rotor37Tube r;
    const auto sec = section(r.pb, r.tube, r.l);
    auto zigzag = [](const yams::api::B2BOptions &opt, const yams::api::StreamTube &tube, const yams::api::BladeSection &s)
    {
        const auto c = yams::api::make_b2b_case(tube, 0, s, opt);
        auto solver = yams::api::make_b2b_solver(c);
        CHECK(solver.curvatureMethod() == opt.curvature);
        solver.computeW();
        CHECK(solver.converged());
        const auto &th = solver.mesh().TH;
        const auto [ni, nj] = solver.dimensions();
        double alt = 0., var = 0.;
        for (size_t i = 1; i + 1 < ni; i++)
            for (size_t j = 2; j + 2 < nj; j++)
            {
                auto e = [&](size_t k) { return th[i + ni * k] - 0.5 * (th[i + ni * (k - 1)] + th[i + ni * (k + 1)]); };
                alt += std::abs(e(j) - 0.5 * (e(j - 1) + e(j + 1)));
                var += std::abs(th[i + ni * (j + 1)] - th[i + ni * (j - 1)]);
            }
        return alt / var;
    };
    yams::api::B2BOptions opt;
    opt.curvature = yams::CurvatureMethod::DoubleDifferentiation;
    const auto dd = zigzag(opt, r.tube, sec);
    opt.curvature = yams::CurvatureMethod::Blended;
    const auto blended = zigzag(opt, r.tube, sec);
    opt.curvature = yams::CurvatureMethod::ThreePoint;
    const auto three = zigzag(opt, r.tube, sec);
    CHECK(blended < dd / 3.);
    CHECK(three < blended);

    opt.curvature = yams::CurvatureMethod::Blended;
    opt.curvature_blend = 1.5;
    CHECK_THROWS_AS(yams::api::make_b2b_case(r.tube, 0, sec, opt), std::invalid_argument);
}
