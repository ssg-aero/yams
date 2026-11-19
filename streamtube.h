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
        /// Spacing of the planes (constant m), relative to the nodes' across the
        /// passage at the leading edge, r pitch / (n_pitch - 1). The solver's stream
        /// line relocation does not stand planes much closer than the nodes across:
        /// on NASA Rotor 37 at mid span (tests_streamtube, 55 % speed) 0.7 diverges,
        /// 0.9 to 1 may settle on a wrong flow (outlet Tt 12 % low with inlet swirl),
        /// 1.1 to 2 agree within 1.5 deg of outlet angle and 1 % of outlet Tt.
        double plane_spacing{1.25};
        /// Extent upstream and downstream of the blade, fractions of its meridional
        /// chord, cut halfway to a neighbouring row's edge and at the tube's ends
        double upstream{0.5};
        double downstream{0.5};
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
        std::vector<double> dr; ///< stream sheet thickness of each plane, radial (BladeToBladeCurvatureSolver::setDr)
        double mass_flow{};     ///< through the tube, all the passages, kg/s
        double omega{};         ///< rad/s
        double Pt{}, Tt{};      ///< absolute, at the first plane
        /// relative flow angle of the meridional flow at the first plane, rad, which the
        /// grid's stagnation line follows upstream (its first segment averages it with
        /// the second plane's)
        double beta_in{};
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
     * The planes are at constant m, uniform upstream, along the blade and downstream,
     * spaced as B2BOptions::plane_spacing asks. Upstream, the initial stagnation line
     * follows the meridional relative flow angle back from the leading edge; its first
     * segment sets the inlet flow angle of the blade-to-blade run (the solver shifts
     * the planes upstream together). Downstream, it leaves the trailing edge along
     * the bisector of the sides, r tan(beta) kept, and the solver's periodicity
     * correction bends it: the meridional angle there comes from the meridional
     * model of the row (its metal angle in DIRECT mode), and a kink at the trailing
     * edge makes the solver diverge. Across the passage the nodes are uniform in theta.
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
        if (!(opt.plane_spacing > 0.))
            throw std::invalid_argument("The blade-to-blade plane spacing must be positive.");
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

        // domain ends: the tube's, halfway to the neighbouring rows' edges
        auto m_in = std::max(m_le - opt.upstream * chord, tube.stations.front().m);
        auto m_out = std::min(m_te + opt.downstream * chord, tube.stations.back().m);
        for (size_t k = 0; k < tube.rows.size(); k++)
        {
            if (k == row)
                continue;
            if (tube.rows[k].m_te <= m_le)
                m_in = std::max(m_in, 0.5 * (tube.rows[k].m_te + m_le));
            if (tube.rows[k].m_le >= m_te)
                m_out = std::min(m_out, 0.5 * (m_te + tube.rows[k].m_le));
        }
        if (!(m_le > tube.stations.front().m) || !(m_te < tube.stations.back().m))
            throw std::invalid_argument("The blade section is not inside the stream tube.");
        if (!(m_in < m_le) || !(m_out > m_te))
            throw std::invalid_argument("No room for the blade-to-blade domain between the blade and a "
                                        "neighbouring row: the rows' edges are too close.");

        B2BCase c;
        c.stream_line = std::make_shared<Curve2d>(tube.mean);
        c.omega = rw.omega;
        c.blade_count = rw.blade_count;
        c.mass_flow = tube.mass_flow;
        c.gas = tube.gas;
        const auto inlet = tube.at(m_in);
        c.Pt = inlet.Pt;
        c.Tt = inlet.Tt;

        const auto dm = opt.plane_spacing * tube.at(m_le).r * pitch / (opt.n_pitch - 1);
        auto count = [dm](double length, size_t n_min) { return std::max(n_min, size_t(std::lround(length / dm))); };
        const auto n_up = count(m_le - m_in, 2), n_blade = count(chord, 2) + 1, n_down = count(m_out - m_te, 2);
        for (size_t k = 0; k < n_up; k++)
            c.m.push_back(m_in + (m_le - m_in) * k / n_up);
        for (size_t k = 0; k < n_blade; k++)
            c.m.push_back(m_le + chord * k / (n_blade - 1));
        for (size_t k = 1; k <= n_down; k++)
            c.m.push_back(m_te + (m_out - m_te) * k / n_down);
        c.n_planes = c.m.size();
        c.j_le = n_up;
        c.j_te = n_up + n_blade - 1;

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
        // downstream: along the bisector of the sides at the trailing edge, r tan(beta) kept
        const auto dm_te = 1e-3 * chord;
        const auto slope_te = 0.5 * ((sa(m_te) - sa(m_te - dm_te)) + (sb(m_te) - sb(m_te - dm_te))) / dm_te;
        const auto rtanb_te = slope_te * tube.at(m_te).r * tube.at(m_te).r;
        th_stag[c.j_te] = sa(m_te);
        for (size_t j = c.j_te + 1; j < c.n_planes; j++)
        {
            auto f = [&](double m) { const auto r = tube.at(m).r; return rtanb_te / (r * r); };
            th_stag[j] = th_stag[j - 1] + 0.5 * (f(c.m[j - 1]) + f(c.m[j])) * (c.m[j] - c.m[j - 1]);
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
        }
        return c;
    }

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
        return solver;
    }
}
