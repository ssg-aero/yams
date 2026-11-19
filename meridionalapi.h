#pragma once
/**
 * @file meridionalapi.h
 * @brief High level entry point of the meridional curvature solver, without VTK.
 *
 * Builds the meridional grid of a channel from its hub and shroud curves and the
 * axial positions of the blade rows leading / trailing edges, runs the solver and
 * returns a status. Meant for applications (e.g. Turbo Blade Forge) that own the
 * geometry and run the solver in a worker thread.
 *
 * The computing stations of the blade rows edges follow the edges' meridional
 * traces when they are given (BladeRow::le_trace / te_trace), else they are lines of
 * constant z; the inlet and outlet stations are lines of constant z, and the stations
 * in between blend their two neighbouring edges' (or the inlet's, the outlet's).
 * The stations keep their shapes while the solver moves the stream lines along them
 * (SolverCase::station_curves).
 *
 * Limits:
 * - the hub and shroud z must increase along the curves, and the stations' ends along
 *   them (axial or mixed flow channel, not radial);
 * - blade rows are in DIRECT mode (metal angle imposed, no deviation model).
 */
#include <curvaturesolver.h>
#include <gas/gasModel.h>
#include <gbs/bscurve.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace yams::api
{
    using Curve2d = gbs::BSCurve<double, 2>;

    /// Blade row: its edges at constant z, or along their meridional traces
    struct BladeRow
    {
        std::string name;
        double z_le{};          ///< leading edge axial position (without le_trace)
        double z_te{};          ///< trailing edge axial position (without te_trace)
        double rpm{};           ///< rotation speed, 0 for a stator
        size_t blade_count{2};
        /// Metal angle (rad, yams convention: beta = atan2(Wu, Vm), Wu = Vu - omega r),
        /// m in [0, 1] from leading to trailing edge, l in [0, 1] from hub to shroud.
        std::function<double(double m, double l)> metal_angle;
        /// Meridional traces of the edges, (z, r) curves from hub to shroud (either
        /// way). When given, the edge's computing station follows the part of the
        /// trace between the hub and the shroud, and z_le / z_te are not used. A trace
        /// going beyond a wall is cut there; one stopping short of it (a tip gap) is
        /// prolonged along its end's tangent.
        std::optional<Curve2d> le_trace;
        std::optional<Curve2d> te_trace;
    };

    /// Inlet conditions as functions of the relative span l in [0, 1] (hub to shroud).
    /// The inlet mean meridional velocity sets the mass flow; Tt and Pt are the
    /// inlet total temperature and pressure. The profiles are applied to the
    /// initial stream lines (evenly spread from hub to shroud), then each stream
    /// line keeps its values when the solver moves it.
    struct Inlet
    {
        /// The mass flow from a mean meridional velocity: vm_mean times the inlet
        /// station's area and its total density (the inlet is set from its totals
        /// before any velocity), not its static one: the flow's mean Vm is higher.
        double vm_mean{30.};
        /// The mass flow itself, kg/s, all the annulus: when set, vm_mean is not used
        /// (a design point is given by its mass flow)
        std::optional<double> mass_flow;
        std::function<double(double)> Tt = [](double) { return 288.15; };
        std::function<double(double)> Pt = [](double) { return 1.01325e5; };
        std::function<double(double)> Vu = [](double) { return 0.; };
    };

    /// Perfect gas: R, cp (gamma follows) and the entropy origin, shared with the
    /// blade-to-blade solver.
    using Gas = PerfectGas<double>;

    struct Options
    {
        size_t n_streams{11};       ///< stream lines from hub to shroud (>= 3)
        /// Indicative number of computing stations from inlet to outlet. The closer
        /// the stations compared with the span, the smaller the relaxation they need
        /// and the more iterations: see local_relaxation, Result::station_spacing.
        size_t n_stations{61};
        size_t max_geom{2000};      ///< maximal number of stream lines relocations
        /// Stream lines relocation relaxation at each station from its aspect ratio
        /// (span over spacing) and meridional Mach number, the optimum of
        /// Wilkinson's stability analysis (GridInfo::RF_local, station_relaxation):
        /// it keeps the solver stable as the stations get closer, its iterations
        /// growing as the aspect ratio squared. When false, relax_factor at every station, which
        /// diverges once the stations are too close for it.
        bool local_relaxation{true};
        /// Its constant (GridInfo::RF_B), > 0; by default the curvature method's
        /// (default_relax_constant(curvature, curvature_blend): 9.6 with the default
        /// curvature). At most max_relax_constant(curvature, curvature_blend), the
        /// method's stability limit (Wilkinson's 96 / (5 |k_min|): 19.2 for
        /// DoubleDifferentiation, 15.36 for the default Blended, 4.8 for ThreePoint):
        /// a constant set for one method is refused with a less stable one.
        std::optional<double> relax_constant;
        /// Its meridional Mach number squared limit (GridInfo::RF_Mm2_max), in [0, 1[.
        /// Not used when not compressible.
        double relax_mach2_max{0.6};
        double relax_factor{0.05};  ///< relaxation at every station, when not local_relaxation
        double tol_rel_pos{1e-5};   ///< convergence on stream lines positions, relative to the span
        bool compressible{true};
        /// Stream lines' curvature (GridInfo::curvature). Blended by default:
        /// DoubleDifferentiation leaves waves on the stream lines downstream of a
        /// leading edge (NASA Rotor 37), which Blended damps; ThreePoint damps them
        /// most, for about 4 times the iterations.
        CurvatureMethod curvature{CurvatureMethod::Blended};
        /// The share of ThreePoint in Blended (GridInfo::curvature_blend), in [0, 1].
        /// 0.2 divides the waves of NASA Rotor 37 by 3 (31 stations) to 8 (81) for
        /// about a tenth more iterations than DoubleDifferentiation; 0.5 by 5 to 11,
        /// for nearly twice the iterations.
        double curvature_blend{0.2};
    };

    using yams::default_relax_constant;
    using yams::max_relax_constant;

    struct Problem
    {
        Curve2d hub;      ///< (z, r) curve
        Curve2d shroud;   ///< (z, r) curve
        std::vector<BladeRow> rows;
        Inlet inlet;
        Gas gas;
        Options options;
    };

    enum class Status
    {
        Converged,     ///< stream lines positions within tol_rel_pos
        NotConverged,  ///< max_geom reached
        Cancelled,     ///< stopped by the progress callback
        Diverged,      ///< non finite values in the solution, or stream lines that could no longer be relocated (Result::divergence)
        /// stream lines positions within tol_rel_pos, but Vm held at the sonic
        /// limit at some points (Result::sonic_points): the mass flow is beyond
        /// the choke of the channel, the solution is not a valid flow
        Choked,
        /// stream lines positions within tol_rel_pos, but the radial equilibrium
        /// gives no positive Vm at some points (Result::floor_points): stream
        /// tubes stagnate or would flow back (blades far off their inlet flow,
        /// too much work for the mass flow), the solution is not a valid flow
        FlowReversal,
    };

    inline const char *to_string(Status status)
    {
        switch (status)
        {
        case Status::Converged: return "converged";
        case Status::NotConverged: return "not converged";
        case Status::Cancelled: return "cancelled";
        case Status::Diverged: return "diverged";
        case Status::Choked: return "choked";
        case Status::FlowReversal: return "flow reversal";
        }
        return "unknown";
    }

    /// Computing stations of a blade row edges
    struct RowStations
    {
        size_t i_le{};
        size_t i_te{};
    };

    struct Result
    {
        Status status{Status::NotConverged};
        size_t iterations{};
        /// Smallest distance between two neighbouring stations along any stream
        /// line of the initial grid, relative to the span there. With a fixed
        /// relaxation (relax_factor 0.05), below about 0.04 the stream lines
        /// oscillate and the solver does not converge; with local_relaxation it
        /// converges, in a number of iterations growing as its inverse squared.
        double station_spacing{};
        std::string divergence;             ///< why, when the stream lines could no longer be relocated
        /// Grid points where the solver holds Vm at the sonic limit (meridional
        /// Mach number 1): Status::Choked when the run converges with any.
        size_t sonic_points{};
        /// Grid points where the radial equilibrium gives no positive Vm, the
        /// integration holding it at its floor in a step (GridInfo::vm_floor; the
        /// final Vm there is small but positive): Status::FlowReversal when the
        /// run converges with any.
        size_t floor_points{};
        std::vector<double> residual_max;   ///< per iteration, relative to the span
        std::vector<double> residual_mean;
        std::vector<RowStations> rows;
        SolverCase<double> solver_case;     ///< solved case: solver_case.gi->g(i, j), i station, j stream line
    };

    /// Progress callback (iteration, residual max, residual mean); return false to cancel.
    using Progress = std::function<bool(size_t, double, double)>;

    namespace detail
    {
        using Vec2 = std::array<double, 2>;

        /// A wall of the channel (the hub or the shroud): a (z, r) curve, z increasing along it.
        class Wall
        {
            const Curve2d &crv_;
            const char *name_;
            double u1_{}, u2_{}, z1_{}, z2_{};

        public:
            /// @throw std::invalid_argument if z does not increase along the curve
            Wall(const Curve2d &crv, const char *name) : crv_{crv}, name_{name}
            {
                auto [u1, u2] = crv.bounds();
                u1_ = u1;
                u2_ = u2;
                const size_t n_check = 64;
                double z_prev = crv(u1)[0];
                for (size_t k = 1; k <= n_check; k++)
                {
                    double z_k = crv(u1 + (u2 - u1) * k / n_check)[0];
                    if (!(z_k > z_prev))
                        throw std::invalid_argument(std::string("yams: the ") + name + " z must increase along the curve (axial or mixed flow channel)");
                    z_prev = z_k;
                }
                z1_ = crv(u1)[0];
                z2_ = crv(u2)[0];
            }
            double z_begin() const { return z1_; }
            double z_end() const { return z2_; }
            /// The point of the wall at axial position z.
            /// @throw std::invalid_argument if z is outside the wall
            Vec2 at(double z) const
            {
                if (!(z >= z1_ && z <= z2_))
                    throw std::invalid_argument(std::string("yams: z = ") + std::to_string(z) + " is outside the " + name_);
                double a = u1_, b = u2_;
                for (int it = 0; it < 200 && (b - a) > 1e-14 * (u2_ - u1_); it++)
                {
                    double c = 0.5 * (a + b);
                    (crv_(c)[0] < z ? a : b) = c;
                }
                auto p = crv_(0.5 * (a + b));
                return {p[0], p[1]};
            }
            /// The radius of the wall at z, at its nearest end's for a z outside it.
            double r(double z) const { return at(std::clamp(z, z1_, z2_))[1]; }
        };

        /// A computing station: its ends on the hub and the shroud at z_hub and
        /// z_shroud, its points (span fractions t_k = k / (n - 1)) the chord's between
        /// them plus dev[k], the station's departure from its chord (none for a
        /// station of constant z).
        struct StationShape
        {
            double z_hub{};
            double z_shroud{};
            std::vector<Vec2> dev;
        };

        inline StationShape constant_z(double z, size_t n) { return {z, z, std::vector<Vec2>(n, Vec2{0., 0.})}; }

        /// The point where the segment a -> b leaves a region: inside(a), not inside(b).
        template <typename Inside>
        inline Vec2 crossing(Vec2 a, Vec2 b, const Inside &inside)
        {
            for (int it = 0; it < 60; it++)
            {
                const Vec2 c{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
                (inside(c) ? a : b) = c;
            }
            return {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
        }

        /// The station of an edge trace: the part of the trace between the hub and the
        /// shroud, sampled at n points evenly spread along its length from the hub.
        /// A trace going beyond a wall is cut there; one stopping short of a wall (a
        /// tip gap) is prolonged along its end's tangent, by its own length at most.
        /// @throw std::invalid_argument if the trace does not cross the channel
        inline StationShape trace_station(const Curve2d &trace, const Wall &hub, const Wall &shroud, size_t n)
        {
            auto [u1, u2] = trace.bounds();
            const size_t n_trace = 400;
            std::vector<Vec2> pts(n_trace + 1);
            double length = 0.;
            for (size_t k = 0; k <= n_trace; k++)
            {
                auto p = trace(u1 + (u2 - u1) * k / n_trace);
                pts[k] = {p[0], p[1]};
                if (k > 0)
                    length += std::hypot(pts[k][0] - pts[k - 1][0], pts[k][1] - pts[k - 1][1]);
            }
            if (!(length > 0.))
                throw std::invalid_argument("yams: an edge trace has no length");
            // From the hub to the shroud: by its ends' positions across the channel,
            // not by their radii (the walls' radii differ from one end's z to the other's)
            auto across = [&](const Vec2 &p)
            {
                const double r_hub = hub.r(p[0]);
                return (p[1] - r_hub) / (shroud.r(p[0]) - r_hub);
            };
            if (across(pts.front()) > across(pts.back()))
                std::reverse(pts.begin(), pts.end());

            // Prolonged at both ends, then cut at the walls
            auto prolonged = [length](const Vec2 &end, const Vec2 &before) -> Vec2
            {
                const double d = std::hypot(end[0] - before[0], end[1] - before[1]);
                if (!(d > 0.))
                    return end;
                return {end[0] + (end[0] - before[0]) * length / d, end[1] + (end[1] - before[1]) * length / d};
            };
            std::vector<Vec2> q;
            q.reserve(pts.size() + 2);
            q.push_back(prolonged(pts.front(), pts[1]));
            q.insert(q.end(), pts.begin(), pts.end());
            q.push_back(prolonged(pts.back(), pts[n_trace - 1]));
            auto above_hub = [&](const Vec2 &p) { return p[1] > hub.r(p[0]); };
            auto below_shroud = [&](const Vec2 &p) { return p[1] < shroud.r(p[0]); };
            const size_t mid = q.size() / 2;
            if (!(above_hub(q[mid]) && below_shroud(q[mid])))
                throw std::invalid_argument("yams: an edge trace must cross the channel from the hub to the shroud");
            size_t k_hub = mid, k_shroud = mid;
            while (k_hub > 0 && above_hub(q[k_hub - 1]))
                k_hub--;
            while (k_shroud + 1 < q.size() && below_shroud(q[k_shroud + 1]))
                k_shroud++;
            if (k_hub == 0 || k_shroud + 1 == q.size())
                throw std::invalid_argument("yams: an edge trace does not reach the hub or the shroud");
            std::vector<Vec2> cut;
            cut.reserve(k_shroud - k_hub + 3);
            cut.push_back(crossing(q[k_hub], q[k_hub - 1], above_hub));
            cut.insert(cut.end(), q.begin() + k_hub, q.begin() + k_shroud + 1);
            cut.push_back(crossing(q[k_shroud], q[k_shroud + 1], below_shroud));

            // Evenly spread along the length
            const size_t n_cut = cut.size() - 1;
            std::vector<double> len(n_cut + 1, 0.);
            for (size_t k = 1; k <= n_cut; k++)
                len[k] = len[k - 1] + std::hypot(cut[k][0] - cut[k - 1][0], cut[k][1] - cut[k - 1][1]);
            if (!(len[n_cut] > 0.))
                throw std::invalid_argument("yams: an edge trace has no length between the hub and the shroud");
            std::vector<Vec2> e(n);
            size_t k = 0;
            for (size_t j = 0; j < n; j++)
            {
                const double s = len[n_cut] * j / (n - 1);
                while (k + 1 < n_cut && len[k + 1] < s)
                    k++;
                const double dl = len[k + 1] - len[k];
                const double w = dl > 0. ? std::clamp((s - len[k]) / dl, 0., 1.) : 0.;
                e[j] = {cut[k][0] + w * (cut[k + 1][0] - cut[k][0]), cut[k][1] + w * (cut[k + 1][1] - cut[k][1])};
            }
            StationShape st{e.front()[0], e.back()[0], std::vector<Vec2>(n)};
            for (size_t j = 0; j < n; j++)
            {
                const double t = double(j) / (n - 1);
                for (int c = 0; c < 2; c++)
                    st.dev[j][c] = e[j][c] - (e.front()[c] + t * (e.back()[c] - e.front()[c]));
            }
            return st;
        }

        /// The station a fraction s of the way from a to b.
        inline StationShape blend(const StationShape &a, const StationShape &b, double s)
        {
            StationShape st{(1 - s) * a.z_hub + s * b.z_hub, (1 - s) * a.z_shroud + s * b.z_shroud, a.dev};
            for (size_t j = 0; j < st.dev.size(); j++)
                for (int c = 0; c < 2; c++)
                    st.dev[j][c] = (1 - s) * a.dev[j][c] + s * b.dev[j][c];
            return st;
        }

        /// The points of a station on the channel.
        inline std::vector<Vec2> station_points(const Wall &hub, const Wall &shroud, const StationShape &st)
        {
            // z from the station's own: a station of constant z exactly at its z
            const Vec2 h{st.z_hub, hub.at(st.z_hub)[1]}, s{st.z_shroud, shroud.at(st.z_shroud)[1]};
            if (!(s[1] > h[1]))
                throw std::invalid_argument("yams: the shroud must be above the hub");
            const size_t n = st.dev.size();
            std::vector<Vec2> pts(n);
            for (size_t j = 0; j < n; j++)
            {
                const double t = double(j) / (n - 1);
                for (int c = 0; c < 2; c++)
                    pts[j][c] = h[c] + t * (s[c] - h[c]) + st.dev[j][c];
            }
            return pts;
        }

        /// Station counts per block between the hard z positions, at least 2 intervals each.
        inline std::vector<size_t> block_intervals(const std::vector<double> &z_hard, size_t n_stations)
        {
            const size_t n_blocks = z_hard.size() - 1;
            const double length = z_hard.back() - z_hard.front();
            size_t n_intervals = std::max(n_stations, 2 * n_blocks + 1) - 1;
            std::vector<size_t> counts(n_blocks);
            for (size_t b = 0; b < n_blocks; b++)
                counts[b] = std::max<size_t>(2, std::lround(n_intervals * (z_hard[b + 1] - z_hard[b]) / length));
            return counts;
        }
    }

    /**
     * @brief Build the solver case of a problem.
     *
     * @param rows_stations if not null, receives the stations of each row edges
     * @throw std::invalid_argument for an inconsistent problem (rows order, z outside the channel...)
     */
    inline SolverCase<double> make_case(const Problem &pb, std::vector<RowStations> *rows_stations = nullptr)
    {
        const auto &opt = pb.options;
        if (opt.n_streams < 3)
            throw std::invalid_argument("yams: at least 3 stream lines are needed");
        if (!(opt.curvature_blend >= 0. && opt.curvature_blend <= 1.))
            throw std::invalid_argument("yams: curvature_blend must be between 0 and 1");
        const double relax_max = max_relax_constant(opt.curvature, opt.curvature_blend);
        const double relax_constant = opt.relax_constant.value_or(default_relax_constant(opt.curvature, opt.curvature_blend));
        if (opt.local_relaxation &&
            (!(relax_constant > 0.) || !(opt.relax_mach2_max >= 0. && opt.relax_mach2_max < 1.)))
            throw std::invalid_argument("yams: the local relaxation needs relax_constant > 0 and 0 <= relax_mach2_max < 1");
        if (opt.local_relaxation && relax_constant > relax_max)
            throw std::invalid_argument("yams: relax_constant " + std::to_string(relax_constant) +
                                        " is above the stability limit of the curvature method, " +
                                        std::to_string(relax_max) +
                                        " (leave it unset for the method's default)");

        const size_t nj = opt.n_streams;
        // The stations' shapes are held on more points than there are stream lines
        // (about a hundred intervals): the grid points are every n_sub-th of them.
        const size_t n_sub = (100 + nj - 2) / (nj - 1);
        const size_t n_shape = (nj - 1) * n_sub + 1;
        const detail::Wall hub{pb.hub, "hub"}, shroud{pb.shroud, "shroud"};

        // Hard stations: inlet, rows edges, outlet
        double z_in = std::max(hub.z_begin(), shroud.z_begin());
        double z_out = std::min(hub.z_end(), shroud.z_end());
        std::vector<detail::StationShape> hard{detail::constant_z(z_in, n_shape)};
        for (const auto &row : pb.rows)
        {
            if (!row.metal_angle)
                throw std::invalid_argument("yams: blade row " + row.name + " has no metal angle");
            hard.push_back(row.le_trace ? detail::trace_station(*row.le_trace, hub, shroud, n_shape) : detail::constant_z(row.z_le, n_shape));
            hard.push_back(row.te_trace ? detail::trace_station(*row.te_trace, hub, shroud, n_shape) : detail::constant_z(row.z_te, n_shape));
        }
        hard.push_back(detail::constant_z(z_out, n_shape));
        // In order along the hub, the shroud and in between
        std::vector<std::vector<detail::Vec2>> hard_pts;
        for (const auto &st : hard)
            hard_pts.push_back(detail::station_points(hub, shroud, st));
        for (size_t k = 1; k < hard.size(); k++)
            for (size_t d = 0; d < n_shape; d++)
                if (!(hard_pts[k][d][0] > hard_pts[k - 1][d][0]))
                    throw std::invalid_argument("yams: blade rows must be inside the channel, in axial order, with their leading edge before their trailing edge and no overlap");

        // Stations, counted per block along the mean of the hub and shroud z
        std::vector<double> z_mean;
        for (const auto &st : hard)
            z_mean.push_back(0.5 * (st.z_hub + st.z_shroud));
        auto counts = detail::block_intervals(z_mean, opt.n_stations);
        std::vector<std::vector<detail::Vec2>> pts{hard_pts.front()};
        std::vector<size_t> i_hard{0};
        for (size_t b = 0; b < counts.size(); b++)
        {
            // the last station of a block is the hard one itself: recomputing it can
            // land one ulp past it, i.e. outside the channel at the outlet
            for (size_t k = 1; k < counts[b]; k++)
                pts.push_back(detail::station_points(hub, shroud, detail::blend(hard[b], hard[b + 1], double(k) / counts[b])));
            pts.push_back(hard_pts[b + 1]);
            i_hard.push_back(pts.size() - 1);
        }
        const size_t ni = pts.size();

        // Grid: the stations, stream lines evenly spread along them from hub to shroud.
        // A leaning station's chord can cross a wall that is not straight: refused.
        MeridionalGrid<double> g(ni, nj);
        for (size_t i = 0; i < ni; i++)
        {
            for (size_t j = 0; j < nj; j++)
            {
                const auto &p = pts[i][j * n_sub];
                if (j > 0 && j + 1 < nj && !(p[1] > hub.r(p[0]) && p[1] < shroud.r(p[0])))
                    throw std::invalid_argument("yams: a computing station leaves the channel between the hub and the shroud (an edge leaning too much for the shape of the walls)");
                g(i, j).x = p[0];
                g(i, j).y = p[1];
            }
        }

        // Blade rows
        std::vector<BladeInfo<double>> blades;
        std::vector<RowStations> stations;
        for (size_t r = 0; r < pb.rows.size(); r++)
        {
            const auto &row = pb.rows[r];
            BladeInfo<double> b;
            b.name = row.name;
            b.i1 = static_cast<int>(i_hard[2 * r + 1]);
            b.i2 = static_cast<int>(i_hard[2 * r + 2]);
            b.is = (b.i1 + b.i2) / 2;
            b.omg = row.rpm * std::numbers::pi / 30.;
            b.z_ = row.blade_count;
            b.mode = MeridionalBladeMode::DIRECT;
            b.k = row.metal_angle;
            blades.push_back(b);
            stations.push_back({size_t(b.i1), size_t(b.i2)});
        }
        auto sc = make_solver_case<double>(std::move(g), blades);

        // The stations keep their shapes: the stream lines move along them
        std::vector<double> u_shape(n_shape), u_grid(nj);
        for (size_t d = 0; d < n_shape; d++)
            u_shape[d] = double(d) / (n_shape - 1);
        for (size_t j = 0; j < nj; j++)
            u_grid[j] = u_shape[j * n_sub];
        for (size_t i = 0; i < ni; i++)
            sc.station_curves.push_back(gbs::interpolate(pts[i], u_shape, 1));
        sc.station_params.assign(ni, u_grid);

        // Boundary conditions, gas and solver settings
        // INLET_VmMoy_Ts_Ps_Vu sets the inlet totals from its Ts / Ps before the
        // inlet velocity is computed (Vm = 0): they are the total values, which
        // the solver keeps. (A static inlet would need its own formulation.)
        if (pb.inlet.mass_flow)
        {
            // the same inlet, its mass flow given (apply_mf keeps it)
            if (!(*pb.inlet.mass_flow > 0.) || !std::isfinite(*pb.inlet.mass_flow))
                throw std::invalid_argument("The inlet mass flow must be positive.");
            sc.inlet.mode = MeridionalBC::INLET_Mf_Ts_Ps_Vu;
            sc.inlet.Mf = *pb.inlet.mass_flow;
            // the inlet totals come from its Ts / Ps before any velocity (Vm = 0, as
            // with vm_mean): no first guess of Vm here, it would raise them
        }
        else
        {
            sc.inlet.mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
            sc.inlet.Vm_moy = pb.inlet.vm_mean;
        }
        sc.inlet.Ts = pb.inlet.Tt;
        sc.inlet.Ps = pb.inlet.Pt;
        sc.inlet.Vu = pb.inlet.Vu;
        sc.gi->R = pb.gas.R;
        sc.gi->Tref = pb.gas.Tref;
        sc.gi->Pref = pb.gas.Pref;
        sc.gi->rho_cst = !opt.compressible;
        sc.gi->RF = opt.relax_factor;
        sc.gi->RF_local = opt.local_relaxation;
        sc.gi->RF_B = relax_constant;
        sc.gi->curvature = opt.curvature;
        sc.gi->curvature_blend = opt.curvature_blend;
        // the density follows the flow from the first iteration: held until the
        // stream lines settle (the solver's default), a run that stalls before
        // that never converges (a flow near choke downstream of a rotor)
        sc.gi->rho_update_residual = std::numeric_limits<double>::max();
        sc.gi->RF_Mm2_max = opt.relax_mach2_max;
        for (auto &gp : *sc.gi->g)
            gp.Cp = pb.gas.cp;
        sc.max_geom = opt.max_geom;
        sc.tol_rel_pos = opt.tol_rel_pos;
        sc.relocate = true;
        sc.mf_uniform = true;

        if (rows_stations)
            *rows_stations = stations;
        return sc;
    }

    /// Smallest distance between two neighbouring stations of a case along any
    /// stream line, relative to the span of the station (see
    /// Result::station_spacing): the inverse of the largest aspect ratio that
    /// station_relaxation sees.
    inline double station_spacing(const SolverCase<double> &sc)
    {
        const auto &g = *sc.gi->g;
        const size_t nj = g.nCols();
        double spacing = std::numeric_limits<double>::max();
        for (size_t i = 0; i < g.nRows(); i++)
            spacing = std::min(spacing, yams::station_min_spacing(g, i) / g(i, nj - 1).l);
        return spacing;
    }

    /**
     * @brief Solve a problem.
     *
     * Blocking: run it in a worker thread from a GUI. The progress callback is called
     * from the calling thread after each stream lines relocation.
     */
    inline Result solve(const Problem &pb, const Progress &progress = {})
    {
        Result res;
        res.solver_case = make_case(pb, &res.rows);
        auto &sc = res.solver_case;
        res.station_spacing = station_spacing(sc);
        if (progress)
            sc.on_iteration = [&progress](size_t it, double r_max, double r_moy) { return progress(it, r_max, r_moy); };

        curvature_solver(sc);
        sc.on_iteration = nullptr; // do not keep a reference to the callback

        res.iterations = sc.log.delta_pos_max.size();
        res.divergence = sc.divergence;
        res.residual_max = sc.log.delta_pos_max;
        res.residual_mean = sc.log.delta_pos_moy;

        const auto &g = *sc.gi->g;
        bool finite = std::all_of(g.begin(), g.end(), [](const auto &gp)
                                  { return std::isfinite(gp.Vm) && std::isfinite(gp.Ps) && std::isfinite(gp.Ts) && std::isfinite(gp.Vu); });
        // the solver clips Vm to sqrt(gamma R Ts): a point there is choked
        const double R = sc.gi->R;
        res.sonic_points = static_cast<size_t>(std::count_if(g.begin(), g.end(), [R](const auto &gp)
                                  { return gp.Vm >= (1. - 1e-3) * std::sqrt(gp.ga * R * gp.Ts); }));
        const auto &floor = sc.gi->vm_floor;
        res.floor_points = static_cast<size_t>(std::count(floor.begin(), floor.end(), 1));
        if (!finite || sc.diverged)
            res.status = Status::Diverged;
        else if (sc.cancelled)
            res.status = Status::Cancelled;
        else if (sc.converged)
            res.status = res.floor_points > 0 ? Status::FlowReversal
                       : res.sonic_points > 0 ? Status::Choked
                                              : Status::Converged;
        else
            res.status = Status::NotConverged;
        return res;
    }
}
