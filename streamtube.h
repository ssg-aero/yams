#pragma once
/**
 * @file streamtube.h
 * @brief A stream tube of a meridional result, and the blade-to-blade case of a
 * blade row along it.
 *
 * The stream tube around a stream line of a meridional run carries what a
 * blade-to-blade analysis of that stream line needs, whatever its solver:
 * - the mean stream line and the two stream surfaces that bound the tube, (z, r)
 *   curves (a thin slab between the bounding surfaces is the domain of a
 *   quasi-3D Navier-Stokes run, SU2 on a QuadSmith mesh for instance);
 * - its thickness normal to the mean stream line, its mass flow;
 * - the meridional flow along the mean stream line (absolute frame), the gas,
 *   the blade rows' edges and rotation speeds.
 *
 * make_b2b_case builds from it, and from the two sides of a blade section on the
 * stream surface, the grid and the conditions of yams' blade-to-blade solver
 * (BladeToBladeCurvatureSolver); make_b2b_solver makes the solver.
 *
 * Coordinates on the stream surface: m is the meridional arc length along the
 * tube's mean stream line from the meridional inlet (GridPoint::m), theta the
 * angle about the axis, in rad, positive in the direction of positive Vu (yams'
 * convention, beta = atan2(Wu, Vm) = atan(r dtheta / dm) along a stream line).
 */
#include <meridionalapi.h>
#include <baldeToBlade/bladeToBladeCurvatureSolver.h>
#include <gbs/bscinterp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace yams::api
{
    /// The meridional flow at a station of a stream tube, on its mean stream line
    struct StreamTubeStation
    {
        double z{}, r{};  ///< position, m
        double m{};       ///< meridional arc length from the meridional inlet, m
        double b{};       ///< tube thickness normal to the mean stream line, m
        double phi{};     ///< mean stream line slope, atan2(dr, dz), rad
        double Vm{}, Vu{}; ///< absolute velocity components, m/s
        double Pt{}, Tt{}, Ps{}, Ts{}, rho{};
    };

    /// A blade row as the stream tube crosses it
    struct StreamTubeRow
    {
        std::string name;
        double m_le{}, m_te{}; ///< edges on the mean stream line (computing stations)
        double omega{};        ///< rotation speed, rad/s (Wu = Vu - omega r)
        size_t blade_count{};
    };

    struct StreamTube
    {
        size_t stream_line{};  ///< index of the meridional stream line, 0 at the hub
        Curve2d mean;          ///< (z, r) mean stream line, parameter ~ m
        Curve2d lower, upper;  ///< (z, r) bounding stream surfaces, hub side and shroud side
        std::vector<StreamTubeStation> stations; ///< at the meridional computing stations
        double mass_flow{};    ///< through the tube, all the annulus, kg/s
        Gas gas;
        std::vector<StreamTubeRow> rows;

        /// Station values at m, linear between the computing stations (clamped to the ends)
        StreamTubeStation at(double m) const
        {
            if (stations.empty())
                throw std::logic_error("Empty stream tube.");
            if (m <= stations.front().m)
                return stations.front();
            if (m >= stations.back().m)
                return stations.back();
            auto it = std::ranges::upper_bound(stations, m, {}, &StreamTubeStation::m);
            const auto &a = *std::prev(it), &c = *it;
            const auto t = (m - a.m) / (c.m - a.m);
            // every field by name: a new field of StreamTubeStation is added here
            using S = StreamTubeStation;
            S s;
            for (auto f : {&S::z, &S::r, &S::m, &S::b, &S::phi, &S::Vm, &S::Vu, &S::Pt, &S::Tt, &S::Ps, &S::Ts, &S::rho})
                s.*f = a.*f + t * (c.*f - a.*f);
            return s;
        }
    };

    namespace detail
    {
        inline Curve2d curve_through(const std::vector<std::array<double, 2>> &pts, const std::vector<double> &u)
        {
            return gbs::interpolate(gbs::points_vector<double, 2>(pts.begin(), pts.end()), u,
                                    std::min<size_t>(3, pts.size() - 1));
        }
    }

    /**
     * @brief The stream tube around the stream line j of a meridional result.
     *
     * Its bounding stream surfaces lie halfway in mass flow between the stream line
     * and its neighbours (the hub or the shroud for the first and the last), so the
     * tubes of all the stream lines share the channel and its mass flow. The
     * positions are linear in the cumulated mass flow along each station.
     *
     * The curves' parameter is the arc length of the mean curve (StreamTubeStation::m),
     * the abscissa BladeToBladeCurvatureSolver computes from it; it differs from the
     * meridional grid's m, the length of the polyline through the nodes.
     *
     * @throw std::invalid_argument for a result that is not a valid flow (any status
     * but Status::Converged: choked, flow reversal, not converged...)
     */
    inline StreamTube stream_tube(const Problem &pb, const Result &res, size_t j)
    {
        if (res.status != Status::Converged)
            throw std::invalid_argument(std::string("The meridional result is not a valid flow (") +
                                        to_string(res.status) + ").");
        if (!res.solver_case.gi || !res.solver_case.gi->g)
            throw std::invalid_argument("The result holds no meridional grid.");
        const auto &g = *res.solver_case.gi->g;
        const size_t n_st = g.nRows(), n_sl = g.nCols();
        if (j >= n_sl)
            throw std::invalid_argument("No such stream line: " + std::to_string(j) + ".");
        if (n_st < 2)
            throw std::invalid_argument("The meridional grid has less than 2 stations.");

        // point at the cumulated mass flow q along station i
        auto at_q = [&](size_t i, double q) -> std::array<double, 2>
        {
            size_t k = 1;
            while (k < n_sl - 1 && g(i, k).q < q)
                k++;
            const auto &a = g(i, k - 1), &c = g(i, k);
            const auto t = c.q != a.q ? std::clamp((q - a.q) / (c.q - a.q), 0., 1.) : 0.;
            return {a.x + t * (c.x - a.x), a.y + t * (c.y - a.y)};
        };

        StreamTube tube;
        tube.stream_line = j;
        tube.gas = pb.gas;
        std::vector<std::array<double, 2>> p_mean, p_low, p_up;
        std::vector<double> u;
        for (size_t i = 0; i < n_st; i++)
        {
            const auto &gp = g(i, j);
            const auto q_low = j > 0 ? 0.5 * (g(i, j - 1).q + gp.q) : g(i, 0).q;
            const auto q_up = j + 1 < n_sl ? 0.5 * (gp.q + g(i, j + 1).q) : g(i, n_sl - 1).q;
            p_mean.push_back({gp.x, gp.y});
            p_low.push_back(at_q(i, q_low));
            p_up.push_back(at_q(i, q_up));
            u.push_back(gp.m);
            if (i == 0)
                tube.mass_flow = q_up - q_low;

            StreamTubeStation st;
            st.z = gp.x;
            st.r = gp.y;
            st.m = gp.m;
            st.Vm = gp.Vm;
            st.Vu = gp.Vu;
            st.Pt = gp.Pt;
            st.Tt = gp.Tt;
            st.Ps = gp.Ps;
            st.Ts = gp.Ts;
            st.rho = gp.rho;
            tube.stations.push_back(st);
        }
        // parameter: the mean curve's own arc length at the nodes (a fixed point,
        // the polyline length to start with)
        for (int pass = 0; pass < 3; pass++)
        {
            const auto crv = detail::curve_through(p_mean, u);
            std::vector<double> s{u.front()};
            for (size_t i = 1; i < n_st; i++)
                s.push_back(s.back() + gbs::length(crv, u[i - 1], u[i]));
            u = std::move(s);
        }
        for (size_t i = 0; i < n_st; i++)
            tube.stations[i].m = u[i];
        tube.mean = detail::curve_through(p_mean, u);
        tube.lower = detail::curve_through(p_low, u);
        tube.upper = detail::curve_through(p_up, u);

        // slope and thickness normal to the mean stream line
        for (size_t i = 0; i < n_st; i++)
        {
            auto &st = tube.stations[i];
            const auto d = tube.mean(st.m, 1);
            st.phi = std::atan2(d[1], d[0]);
            const auto n = std::hypot(d[0], d[1]);
            const std::array<double, 2> w{p_up[i][0] - p_low[i][0], p_up[i][1] - p_low[i][1]};
            st.b = std::abs(w[0] * d[1] - w[1] * d[0]) / n;
        }

        for (size_t k = 0; k < pb.rows.size() && k < res.rows.size(); k++)
        {
            const auto &row = pb.rows[k];
            tube.rows.push_back({row.name, tube.stations[res.rows[k].i_le].m, tube.stations[res.rows[k].i_te].m,
                                 row.rpm * std::numbers::pi / 30., row.blade_count});
        }
        return tube;
    }

    /// A blade section on the stream surface of a tube: its two sides, (m, theta)
    /// curves from the leading to the trailing edge (m along the tube's mean stream
    /// line, see the file's description), m increasing along each.
    struct BladeSection
    {
        Curve2d side1, side2;
    };

    struct B2BOptions
    {
        size_t n_pitch{11}; ///< nodes across the passage, sides included (>= 3)
        /// Extent upstream and downstream of the blade, fractions of its meridional
        /// chord, cut halfway to a neighbouring row's edge and at the tube's ends; the
        /// planes are the meridional stations within it, two at least on each side.
        /// The stagnation lines need room to settle: on Rotor 37's real blade (cut in
        /// turbo_blade_forge, 31 stations), 0.3, 1 and 2 chords give outlet angles
        /// within 1.1 deg, 0.5 a wrong flow (-8 deg)
        double upstream{1.};
        double downstream{1.};
        /// Added to the meridional run's downstream far field angle, rad (b2b_kutta
        /// searches the one that satisfies the Kutta condition)
        double outlet_offset{0.};
        /// The smallest B2BCase::spacing_ratio accepted: below about 1.1 the solver
        /// diverges or settles on a wrong flow (B2BCase::spacing_ratio). More nodes
        /// across the passage, or a meridional run with fewer stations, raise it.
        double min_spacing_ratio{1.1};
        /// The blade-to-blade stream lines' curvature (BladeToBladeCurvatureSolver::
        /// setCurvature). DoubleDifferentiation by default, as the solver's: Blended
        /// sees the zigzag from one plane to the next that it leaves (NASA Rotor 37 at
        /// 55 % speed: the zigzag of theta divided by 5, in 498 iterations instead of
        /// 835), but does not converge near choke (solver.base_novak, relative Mach
        /// 0.98 with the double differentiation), whatever the relaxation
        CurvatureMethod curvature{CurvatureMethod::DoubleDifferentiation};
        double curvature_blend{0.2};
    };

    /// The grid and conditions of BladeToBladeCurvatureSolver for a blade row along a
    /// stream tube
    struct B2BCase
    {
        std::shared_ptr<Curve2d> stream_line; ///< (z, r) mean stream line of the tube
        /// Grid nodes (m, theta) (m is the stream line's parameter), across the passage first
        /// (i = 0 on side_a and the stagnation line, i = n_pitch - 1 one pitch
        /// further), plane after plane from upstream to downstream
        gbs::points_vector<double, 2> points;
        size_t n_planes{}, j_le{}, j_te{};
        std::vector<double> m;  ///< meridional abscissa of each plane
        /// The meridional computing station of each plane (StreamTube::stations): the
        /// planes are the stations, as in Novak and Hearsey's coordinated meridional
        /// and blade-to-blade programs (J. Fluids Eng. 99, 1977), so that each one's
        /// blade-to-blade results (flow angle, specific mass flow ratio) belong to a
        /// meridional station
        std::vector<size_t> station;
        /// The smallest spacing of the planes relative to the nodes' across the passage
        /// at the leading edge, r pitch / (n_pitch - 1). The solver's stream line
        /// relocation does not stand planes much closer than the nodes across (NASA
        /// Rotor 37 at mid span, 55 % speed, uniform planes: 0.7 diverges, 0.9 to 1
        /// may settle on a wrong flow, 1.1 to 2 agree): fewer stations, or more nodes
        /// across, raise it. 1.3 on Rotor 37 with 31 stations and 11 nodes.
        double spacing_ratio{};
        std::vector<double> dr; ///< stream sheet thickness of each plane, radial (BladeToBladeCurvatureSolver::setDr)
        double mass_flow{};     ///< through the tube, all the passages, kg/s
        double omega{};         ///< rad/s
        double Pt{}, Tt{};      ///< absolute, at the first plane
        /// relative flow angle of the meridional flow at the first plane, rad, which the
        /// grid's stagnation line follows upstream (its first segment averages it with
        /// the second plane's)
        double beta_in{};
        /// relative flow angle of the meridional flow at the last plane, rad: the
        /// downstream far field angle the grid imposes (its last segment's)
        double beta_out{};
        std::vector<double> b;  ///< tube thickness normal to the mean stream line at each plane
        CurvatureMethod curvature{CurvatureMethod::DoubleDifferentiation}; ///< B2BOptions::curvature
        double curvature_blend{0.2};
        size_t blade_count{};
        Gas gas;
    };

    namespace detail
    {
        /// theta(m) along a side, from a dense sampling of the curve
        struct SideLaw
        {
            std::vector<double> m, th;
            explicit SideLaw(const Curve2d &c, size_t n = 401)
            {
                auto [u1, u2] = c.bounds();
                for (size_t k = 0; k < n; k++)
                {
                    auto p = c(u1 + (u2 - u1) * k / (n - 1));
                    m.push_back(p[0]);
                    th.push_back(p[1]);
                }
                for (size_t k = 1; k < n; k++)
                    if (!(m[k] > m[k - 1]))
                        throw std::invalid_argument("A blade side's m does not increase from the leading to the trailing edge.");
            }
            double operator()(double x) const
            {
                x = std::clamp(x, m.front(), m.back());
                auto it = std::ranges::upper_bound(m, x);
                if (it == m.end())
                    return th.back();
                auto k = static_cast<size_t>(std::distance(m.begin(), it));
                auto t = (x - m[k - 1]) / (m[k] - m[k - 1]);
                return th[k - 1] + t * (th[k] - th[k - 1]);
            }
        };
    }

    /**
     * @brief The blade-to-blade case of row `row` along `tube`.
     *
     * The planes are the meridional computing stations crossing the tube (B2BCase::
     * station), the row's leading and trailing edges at its edge stations, as Novak
     * and Hearsey's blade-to-blade program takes the meridional program's stations.
     * The edge planes lie at the section's edges, which must be within a fifth of the
     * neighbouring stations' spacing of the edge stations. Upstream, the initial stagnation line
     * follows the meridional relative flow angle back from the leading edge; its first
     * segment sets the inlet flow angle of the blade-to-blade run (the solver shifts
     * the planes upstream together). Downstream, the far field angle is the meridional
     * run's (B2BCase::beta_out), as in Novak and Hearsey's system, where the meridional
     * program gives the blade-to-blade one its downstream air angle and the Kutta
     * condition at the trailing edge checks it (b2b_output): the stagnation line leaves
     * the trailing edge along the mean of the sides over the last 5 % of the chord and
     * blends into the meridional angle by the plane before the last. Across the passage
     * the nodes are uniform in theta.
     *
     * @throw std::invalid_argument on a wrong row, section or option, or where the
     * stream line is nearly radial (|cos phi| < 0.2): the solver's sheet thickness is
     * radial (dr = b / cos phi).
     */
    inline B2BCase make_b2b_case(const StreamTube &tube, size_t row, const BladeSection &section,
                                 const B2BOptions &opt = {})
    {
        if (row >= tube.rows.size())
            throw std::invalid_argument("No such blade row: " + std::to_string(row) + ".");
        if (opt.n_pitch < 3)
            throw std::invalid_argument("Too few blade-to-blade nodes across the passage.");
        if (!(opt.upstream > 0.) || !(opt.downstream > 0.))
            throw std::invalid_argument("The blade-to-blade domain must extend upstream and downstream of the blade.");
        const auto &rw = tube.rows[row];
        if (rw.blade_count < 1)
            throw std::invalid_argument("The blade row has no blades.");

        const detail::SideLaw s1{section.side1}, s2{section.side2};
        const auto m_le = std::max(s1.m.front(), s2.m.front());
        const auto m_te = std::min(s1.m.back(), s2.m.back());
        const auto chord = m_te - m_le;
        if (!(chord > 0.))
            throw std::invalid_argument("The blade section's sides do not overlap in m.");
        const auto m_mid = 0.5 * (m_le + m_te);
        // side_a, at the larger theta, starts the passage
        const bool one_is_a = s1(m_mid) >= s2(m_mid);
        const auto &sa = one_is_a ? s1 : s2;
        const auto &sb = one_is_a ? s2 : s1;
        const auto pitch = 2. * std::numbers::pi / rw.blade_count;
        if (!(sb(m_mid) + pitch > sa(m_mid)))
            throw std::invalid_argument("The blade section is thicker than the pitch.");
        // the grid takes pointed edges: the sides meet at both
        for (auto m_edge : {m_le, m_te})
            if (std::abs(sa(m_edge) - sb(m_edge)) > 1e-3 * pitch)
                throw std::invalid_argument("The blade section's sides do not meet at m = " + std::to_string(m_edge) +
                                            ": close the section to a point at both edges.");

        // domain ends: the tube's, halfway to the neighbouring rows' edges (hard
        // limits, never crossed), and the extent asked within them
        auto m_in_max = tube.stations.front().m, m_out_max = tube.stations.back().m;
        for (size_t k = 0; k < tube.rows.size(); k++)
        {
            if (k == row)
                continue;
            if (tube.rows[k].m_te <= m_le)
                m_in_max = std::max(m_in_max, 0.5 * (tube.rows[k].m_te + m_le));
            if (tube.rows[k].m_le >= m_te)
                m_out_max = std::min(m_out_max, 0.5 * (m_te + tube.rows[k].m_le));
        }
        const auto m_in = std::max(m_le - opt.upstream * chord, m_in_max);
        const auto m_out = std::min(m_te + opt.downstream * chord, m_out_max);
        if (!(m_le > tube.stations.front().m) || !(m_te < tube.stations.back().m))
            throw std::invalid_argument("The blade section is not inside the stream tube.");
        if (!(m_in < m_le) || !(m_out > m_te))
            throw std::invalid_argument("No room for the blade-to-blade domain between the blade and a "
                                        "neighbouring row: the rows' edges are too close.");

        three_point_share(opt.curvature, opt.curvature_blend); // throws for a wrong blend
        B2BCase c;
        c.curvature = opt.curvature;
        c.curvature_blend = opt.curvature_blend;
        c.stream_line = std::make_shared<Curve2d>(tube.mean);
        c.omega = rw.omega;
        c.blade_count = rw.blade_count;
        c.mass_flow = tube.mass_flow;
        c.gas = tube.gas;

        // the planes: the meridional stations, the row's edges at its edge stations
        const auto &st = tube.stations;
        auto nearest = [&st](double m)
        {
            size_t k = 0;
            for (size_t i = 1; i < st.size(); i++)
                if (std::abs(st[i].m - m) < std::abs(st[k].m - m))
                    k = i;
            return k;
        };
        const size_t i_le = nearest(rw.m_le), i_te = nearest(rw.m_te);
        if (!(i_te > i_le + 1))
            throw std::invalid_argument("The blade row spans less than two meridional station intervals.");
        for (auto [m_sec, i_edge] : {std::pair{m_le, i_le}, std::pair{m_te, i_te}})
        {
            const auto gap = std::min(i_edge > 0 ? st[i_edge].m - st[i_edge - 1].m : 1e300,
                                      i_edge + 1 < st.size() ? st[i_edge + 1].m - st[i_edge].m : 1e300);
            if (std::abs(m_sec - st[i_edge].m) > 0.2 * gap)
                throw std::invalid_argument("The blade section's edge at m = " + std::to_string(m_sec) +
                                            " is off the row's edge station at m = " + std::to_string(st[i_edge].m) + ".");
        }
        size_t i_in = i_le, i_out = i_te;
        // the stations within the extent asked, two at least on each side, but never
        // past the tube's ends or halfway to a neighbouring row
        while (i_in > 0 && st[i_in - 1].m >= m_in_max && (st[i_in - 1].m >= m_in || i_le - i_in < 2))
            i_in--;
        while (i_out + 1 < st.size() && st[i_out + 1].m <= m_out_max && (st[i_out + 1].m <= m_out || i_out - i_te < 2))
            i_out++;
        if (i_le - i_in < 2 || i_out - i_te < 2)
            throw std::invalid_argument("Less than two meridional stations upstream or downstream of the blade row "
                                        "before the tube's end or halfway to a neighbouring row.");
        for (size_t i = i_in; i <= i_out; i++)
        {
            // the edge planes at the section's edges, where its sides meet: a plane
            // through a round edge, off the point, starts the stagnation line off
            // the blade and the solver settles on a wrong flow (Rotor 37's real
            // blade cut in turbo_blade_forge: outlet angle -7 deg instead of -40)
            c.m.push_back(i == i_le ? m_le : i == i_te ? m_te : st[i].m);
            c.station.push_back(i);
        }
        c.n_planes = c.m.size();
        c.j_le = i_le - i_in;
        c.j_te = i_te - i_in;
        {
            double dm_min = 1e300;
            for (size_t j = 1; j < c.n_planes; j++)
                dm_min = std::min(dm_min, c.m[j] - c.m[j - 1]);
            c.spacing_ratio = dm_min / (tube.at(c.m[c.j_le]).r * pitch / (opt.n_pitch - 1));
            if (c.spacing_ratio < opt.min_spacing_ratio)
                throw std::invalid_argument("The blade-to-blade planes are too close for the nodes across the passage "
                                            "(spacing ratio " + std::to_string(c.spacing_ratio) + ", under " +
                                            std::to_string(opt.min_spacing_ratio) + "): more nodes across, or fewer "
                                            "meridional stations.");
        }
        const auto inlet = tube.at(c.m.front());
        c.Pt = inlet.Pt;
        c.Tt = inlet.Tt;

        // relative flow angle of the meridional flow, dtheta/dm = tan(beta) / r
        auto dth_dm = [&](double m)
        {
            const auto st = tube.at(m);
            return (st.Vu - c.omega * st.r) / st.Vm / st.r;
        };
        c.beta_in = std::atan2(inlet.Vu - c.omega * inlet.r, inlet.Vm);
        // stagnation line: trapezoidal integration from the edges
        std::vector<double> th_stag(c.n_planes);
        th_stag[c.j_le] = sa(m_le);
        for (size_t j = c.j_le; j-- > 0;)
            th_stag[j] = th_stag[j + 1] - 0.5 * (dth_dm(c.m[j]) + dth_dm(c.m[j + 1])) * (c.m[j + 1] - c.m[j]);
        // downstream: the far field angle is the meridional run's (Novak and Hearsey
        // 1977: the blade-to-blade program takes the downstream air angle from the
        // meridional one, and the Kutta condition at the trailing edge tells whether that
        // angle, hence the meridional deviation, is right). The solver moves the last
        // plane with the one before: the last segment's angle is kept, it is the
        // meridional one. From the trailing edge, the stagnation line leaves along the
        // mean of the sides over the last 5 % of the chord (a round edge's tangents are
        // anything) and blends linearly into the meridional angle by the plane before the
        // last: a kink at the trailing edge makes the solver diverge.
        const auto dm_te = 0.05 * chord;
        const auto slope_te = 0.5 * ((sa(m_te) - sa(m_te - dm_te)) + (sb(m_te) - sb(m_te - dm_te))) / dm_te;
        const auto r_te = tube.at(m_te).r;
        const auto m_blend = c.m[c.n_planes - 2];
        // the far field: the meridional angle, offset
        auto dth_dm_far = [&](double m)
        {
            const auto st = tube.at(m);
            return std::tan(std::atan2(st.Vu - c.omega * st.r, st.Vm) + opt.outlet_offset) / st.r;
        };
        auto dth_dm_down = [&](double m)
        {
            const auto r = tube.at(m).r;
            const auto w = std::clamp((m - m_te) / (m_blend - m_te), 0., 1.);
            return (1. - w) * slope_te * r_te * r_te / (r * r) + w * dth_dm_far(m);
        };
        th_stag[c.j_te] = sa(m_te);
        for (size_t j = c.j_te + 1; j < c.n_planes; j++)
            th_stag[j] = th_stag[j - 1] + 0.5 * (dth_dm_down(c.m[j - 1]) + dth_dm_down(c.m[j])) * (c.m[j] - c.m[j - 1]);
        {
            const auto out = tube.at(c.m.back());
            c.beta_out = std::atan2(out.Vu - c.omega * out.r, out.Vm) + opt.outlet_offset;
        }

        for (size_t j = 0; j < c.n_planes; j++)
        {
            const auto m = c.m[j];
            const bool blade = j >= c.j_le && j <= c.j_te;
            const auto th1 = blade ? sa(m) : th_stag[j];
            // at the edges both sides meet: the passage is one pitch wide there
            const auto th2 = blade ? sb(m) + pitch : th_stag[j] + pitch;
            for (size_t i = 0; i < opt.n_pitch; i++)
                c.points.push_back({m, th1 + (th2 - th1) * i / (opt.n_pitch - 1)});

            const auto st = tube.at(m);
            const auto cphi = std::cos(st.phi);
            if (std::abs(cphi) < 0.2)
                throw std::invalid_argument("The stream line is nearly radial at m = " + std::to_string(m) +
                                            ": the blade-to-blade solver takes a radial sheet thickness.");
            c.dr.push_back(st.b / std::abs(cphi));
            c.b.push_back(st.b);
        }
        return c;
    }

    /// What a blade-to-blade solution gives back to the meridional run (Novak and
    /// Hearsey 1977), at each plane of its case
    struct B2BOutput
    {
        /// Kutta condition at the trailing edge: (W_a - W_b) / mean, the relative
        /// speeds of the trailing edge's nodes on side a (i = 0) and side b (i =
        /// n_pitch - 1). Zero when the far field angle (the meridional deviation) is
        /// right; else that angle, hence the meridional deviation, is to be changed
        /// until it closes (Novak and Hearsey's pump: two more runs).
        double kutta{};
        /// relative flow angle on the mid stream line (i = (n_pitch - 1) / 2: the
        /// solver's nodes stay on stream lines), rad
        std::vector<double> beta_mean;
        /// specific mass flow ratio B = (rho Vm)_a / (rho Vm)_m (Novak and Hearsey's eq.
        /// 7): the passage's mass flow over the whole pitch (blade blockage counted)
        /// against the mid stream line's. 1 far from the blade
        std::vector<double> B;
        /// the planes' outlet angle against the imposed one: the solver kept it when
        /// they match (rad)
        double beta_out{};
    };

    /// The blade-to-blade solver of a case, ready for computeW()
    inline BladeToBladeCurvatureSolver<double> make_b2b_solver(const B2BCase &c)
    {
        // the stream line's parameter is m at its nodes (stream_tube)
        BladeToBladeCurvatureSolver<double> solver{c.points, c.n_planes, c.stream_line, c.blade_count, c.j_le, c.j_te};
        solver.setGas(c.gas);
        solver.setRotationSpeed(c.omega);
        solver.setFullBladePassageMassFlow(c.mass_flow);
        solver.setPtIn(c.Pt);
        solver.setTtIn(c.Tt);
        solver.setDr(c.dr);
        solver.setCurvature(c.curvature, c.curvature_blend);
        return solver;
    }

    /// The solution's Kutta condition, mid stream line angles and specific mass flow
    /// ratios (B2BOutput), at the planes of case `c`
    inline B2BOutput b2b_output(const B2BCase &c, const BladeToBladeCurvatureSolver<double> &solver)
    {
        const auto &msh = solver.mesh();
        const auto &dat = solver.data();
        const auto [ni, nj] = solver.dimensions();
        if (nj != c.n_planes)
            throw std::invalid_argument("The solver is not the case's.");
        B2BOutput o;
        const auto wa = dat.W[ni * c.j_te], wb = dat.W[ni - 1 + ni * c.j_te];
        o.kutta = (wa - wb) / (0.5 * (wa + wb));
        const auto i_mid = (ni - 1) / 2;
        // the solver's mass flow: it lowers it when a run gives NaN
        const auto mf_passage = solver.fullBladePassageMassFlow() / c.blade_count;
        const auto pitch = 2. * std::numbers::pi / c.blade_count;
        for (size_t j = 0; j < nj; j++)
        {
            const auto id = i_mid + ni * j;
            o.beta_mean.push_back(msh.BT[id]);
            const auto rho_vm_a = mf_passage / (c.b[j] * msh.R[id] * pitch);
            o.B.push_back(rho_vm_a / (dat.R[id] * dat.Vm[id]));
        }
        o.beta_out = msh.BT[i_mid + ni * (nj - 1)];
        return o;
    }

    /// A blade-to-blade solution satisfying the Kutta condition (b2b_kutta)
    struct B2BKutta
    {
        B2BCase case_;
        BladeToBladeCurvatureSolver<double> solver;
        B2BOutput output;
        /// the far field angle's offset from the meridional run's that closes the
        /// Kutta condition, rad: what the meridional deviation lacks (negative: more
        /// deviation, for a compressor rotor of negative angles)
        double outlet_offset{};
        size_t runs{};       ///< blade-to-blade solutions made
        /// the solution converged at the case's mass flow and |kutta| <= tol; else the
        /// solution is the converged one of smallest |kutta| (or the first, if none)
        bool converged{};
    };

    /**
     * @brief The downstream far field angle that satisfies the Kutta condition, and its
     * solution: secant iterations on B2BOptions::outlet_offset from opt's, until
     * |B2BOutput::kutta| <= tol, as Novak and Hearsey iterated on the discharge air
     * angles of their blade-to-blade runs (J. Fluids Eng. 99, 1977, p. 158).
     */
    inline B2BKutta b2b_kutta(const StreamTube &tube, size_t row, const BladeSection &section, B2BOptions opt = {},
                              double tol = 1e-3, size_t max_runs = 10)
    {
        constexpr auto deg = std::numbers::pi / 180.;
        size_t runs = 0;
        auto run = [&](double offset)
        {
            opt.outlet_offset = offset;
            B2BKutta k{.case_ = make_b2b_case(tube, row, section, opt), .outlet_offset = offset};
            k.solver = make_b2b_solver(k.case_);
            k.solver.computeW();
            k.output = b2b_output(k.case_, k.solver);
            k.converged = k.solver.converged() &&
                          std::abs(k.solver.fullBladePassageMassFlow() - k.case_.mass_flow) <= 1e-9 * k.case_.mass_flow;
            runs++;
            return k;
        };
        auto a = run(opt.outlet_offset);
        // the converged solution of smallest residual: what is returned
        auto best = a;
        // a second point a degree away, where the Kutta residual should get smaller:
        // it decreases as the far field angle increases (Rotor 37 at mid span, 55 %
        // speed: -0.23, -0.35, -0.49 at -1, 0, +1 deg); the secant corrects a wrong guess
        double step = (a.output.kutta > 0. ? 1. : -1.) * deg;
        while (a.converged && std::abs(a.output.kutta) > tol && runs < max_runs)
        {
            auto b = run(a.outlet_offset + step);
            if (!b.converged)
            {
                // a step too far: half of it from the last converged solution
                step *= 0.5;
                continue;
            }
            if (std::abs(b.output.kutta) < std::abs(best.output.kutta))
                best = b;
            const auto dk = b.output.kutta - a.output.kutta;
            if (std::abs(b.output.kutta) <= tol || dk == 0.)
                break;
            step = std::clamp(-b.output.kutta * (b.outlet_offset - a.outlet_offset) / dk, -5. * deg, 5. * deg);
            a = std::move(b);
        }
        best.runs = runs;
        best.converged = best.converged && std::abs(best.output.kutta) <= tol;
        return best;
    }
}
