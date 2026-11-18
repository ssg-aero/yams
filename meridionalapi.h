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
 *
 * Limits:
 * - the hub and shroud z must increase along the curves, and the stations' ends along
 *   them (axial or mixed flow channel, not radial);
 * - blade rows are in DIRECT mode (metal angle imposed, no deviation model).
 */
#include <curvaturesolver.h>
#include <gbs/bscurve.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
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
        /// Meridional traces of the edges, (z, r) curves from hub to shroud (either
        /// way). When given, the edge's computing station follows it, its ends put
        /// on the hub and the shroud at the trace's ends' z, and z_le / z_te are not
        /// used.
        std::optional<Curve2d> le_trace;
        std::optional<Curve2d> te_trace;
        double rpm{};           ///< rotation speed, 0 for a stator
        size_t blade_count{2};
        /// Metal angle (rad, yams convention: beta = atan2(Wu, Vm), Wu = Vu - omega r),
        /// m in [0, 1] from leading to trailing edge, l in [0, 1] from hub to shroud.
        std::function<double(double m, double l)> metal_angle;
    };

    /// Inlet conditions as functions of the relative span l in [0, 1] (hub to shroud).
    /// The inlet mean meridional velocity sets the mass flow; Tt and Pt are the
    /// inlet total temperature and pressure. The profiles are applied to the
    /// initial stream lines (evenly spread from hub to shroud), then each stream
    /// line keeps its values when the solver moves it.
    struct Inlet
    {
        double vm_mean{30.};
        std::function<double(double)> Tt = [](double) { return 288.15; };
        std::function<double(double)> Pt = [](double) { return 1.01325e5; };
        std::function<double(double)> Vu = [](double) { return 0.; };
    };

    /// Perfect gas
    struct Gas
    {
        double R{287.04};
        double cp{1004.5};
    };

    struct Options
    {
        size_t n_streams{11};       ///< stream lines from hub to shroud (>= 3)
        size_t n_stations{61};      ///< indicative number of computing stations from inlet to outlet
        size_t max_geom{300};       ///< maximal number of stream lines relocations
        double relax_factor{0.05};  ///< stream lines relocation relaxation
        double tol_rel_pos{1e-5};   ///< convergence on stream lines positions, relative to the span
        bool compressible{true};
    };

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
        Diverged,      ///< non finite values in the solution
    };

    inline const char *to_string(Status status)
    {
        switch (status)
        {
        case Status::Converged: return "converged";
        case Status::NotConverged: return "not converged";
        case Status::Cancelled: return "cancelled";
        case Status::Diverged: return "diverged";
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
        std::vector<double> residual_max;   ///< per iteration, relative to the span
        std::vector<double> residual_mean;
        std::vector<RowStations> rows;
        SolverCase<double> solver_case;     ///< solved case: solver_case.gi->g(i, j), i station, j stream line
    };

    /// Progress callback (iteration, residual max, residual mean); return false to cancel.
    using Progress = std::function<bool(size_t, double, double)>;

    namespace detail
    {
        /// Parameter of the point of a curve at axial position z (z increasing along the curve).
        inline double param_at_z(const Curve2d &crv, double z, const char *name)
        {
            auto [u1, u2] = crv.bounds();
            const size_t n_check = 64;
            double z_prev = crv(u1)[0];
            for (size_t k = 1; k <= n_check; k++)
            {
                double z_k = crv(u1 + (u2 - u1) * k / n_check)[0];
                if (!(z_k > z_prev))
                    throw std::invalid_argument(std::string("yams: the ") + name + " z must increase along the curve (axial or mixed flow channel)");
                z_prev = z_k;
            }
            if (z < crv(u1)[0] || z > crv(u2)[0])
                throw std::invalid_argument(std::string("yams: z = ") + std::to_string(z) + " is outside the " + name);
            double a = u1, b = u2;
            for (int it = 0; it < 200 && (b - a) > 1e-14 * (u2 - u1); it++)
            {
                double c = 0.5 * (a + b);
                (crv(c)[0] < z ? a : b) = c;
            }
            return 0.5 * (a + b);
        }

        using Vec2 = std::array<double, 2>;

        /// A computing station: its ends on the hub and the shroud at z_hub and
        /// z_shroud, its points (span fractions t_j = j / (nj - 1)) the chord's between
        /// them plus dev[j], the station's departure from its chord (none for a
        /// station of constant z).
        struct StationShape
        {
            double z_hub{};
            double z_shroud{};
            std::vector<Vec2> dev;
        };

        inline StationShape constant_z(double z, size_t nj) { return {z, z, std::vector<Vec2>(nj, Vec2{0., 0.})}; }

        /// The station of an edge trace: sampled at nj points evenly spread along its
        /// length, from its hub end (its end of smaller r) to its shroud end.
        inline StationShape trace_station(const Curve2d &trace, size_t nj)
        {
            auto [u1, u2] = trace.bounds();
            const size_t n = 400;
            std::vector<Vec2> pts(n + 1);
            std::vector<double> len(n + 1, 0.);
            for (size_t k = 0; k <= n; k++)
            {
                auto p = trace(u1 + (u2 - u1) * k / n);
                pts[k] = {p[0], p[1]};
                if (k > 0)
                    len[k] = len[k - 1] + std::hypot(pts[k][0] - pts[k - 1][0], pts[k][1] - pts[k - 1][1]);
            }
            if (!(len[n] > 0.))
                throw std::invalid_argument("yams: an edge trace has no length");
            if (pts.front()[1] > pts.back()[1])
            {
                std::reverse(pts.begin(), pts.end());
                for (auto &l : len)
                    l = len[n] - l;
                std::reverse(len.begin(), len.end());
            }
            std::vector<Vec2> e(nj);
            size_t k = 0;
            for (size_t j = 0; j < nj; j++)
            {
                const double s = len[n] * j / (nj - 1);
                while (k + 1 < n && len[k + 1] < s)
                    k++;
                const double dl = len[k + 1] - len[k];
                const double w = dl > 0. ? std::clamp((s - len[k]) / dl, 0., 1.) : 0.;
                e[j] = {pts[k][0] + w * (pts[k + 1][0] - pts[k][0]), pts[k][1] + w * (pts[k + 1][1] - pts[k][1])};
            }
            StationShape st{e.front()[0], e.back()[0], std::vector<Vec2>(nj)};
            for (size_t j = 0; j < nj; j++)
            {
                const double t = double(j) / (nj - 1);
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
        inline std::vector<Vec2> station_points(const Problem &pb, const StationShape &st)
        {
            auto h = pb.hub(param_at_z(pb.hub, st.z_hub, "hub"));
            auto s = pb.shroud(param_at_z(pb.shroud, st.z_shroud, "shroud"));
            if (!(s[1] > h[1]))
                throw std::invalid_argument("yams: the shroud must be above the hub");
            const size_t nj = st.dev.size();
            std::vector<Vec2> pts(nj);
            // z from the station's own: a station of constant z exactly at its z
            const Vec2 h_ = {st.z_hub, h[1]}, s_ = {st.z_shroud, s[1]};
            for (size_t j = 0; j < nj; j++)
            {
                const double t = double(j) / (nj - 1);
                for (int c = 0; c < 2; c++)
                    pts[j][c] = h_[c] + t * (s_[c] - h_[c]) + st.dev[j][c];
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

        const size_t nj = opt.n_streams;

        // Hard stations: inlet, rows edges, outlet
        double z_in = std::max(pb.hub(pb.hub.bounds()[0])[0], pb.shroud(pb.shroud.bounds()[0])[0]);
        double z_out = std::min(pb.hub(pb.hub.bounds()[1])[0], pb.shroud(pb.shroud.bounds()[1])[0]);
        std::vector<detail::StationShape> hard{detail::constant_z(z_in, nj)};
        for (const auto &row : pb.rows)
        {
            if (!row.metal_angle)
                throw std::invalid_argument("yams: blade row " + row.name + " has no metal angle");
            hard.push_back(row.le_trace ? detail::trace_station(*row.le_trace, nj) : detail::constant_z(row.z_le, nj));
            hard.push_back(row.te_trace ? detail::trace_station(*row.te_trace, nj) : detail::constant_z(row.z_te, nj));
        }
        hard.push_back(detail::constant_z(z_out, nj));
        // In order along the hub, the shroud and every initial stream line
        std::vector<std::vector<detail::Vec2>> hard_pts;
        for (const auto &st : hard)
            hard_pts.push_back(detail::station_points(pb, st));
        for (size_t k = 1; k < hard.size(); k++)
            for (size_t j = 0; j < nj; j++)
                if (!(hard_pts[k][j][0] > hard_pts[k - 1][j][0]))
                    throw std::invalid_argument("yams: blade rows must be inside the channel, in axial order, with their leading edge before their trailing edge and no overlap");

        // Stations, counted per block along the mean of the hub and shroud z
        std::vector<double> z_mean;
        for (const auto &st : hard)
            z_mean.push_back(0.5 * (st.z_hub + st.z_shroud));
        auto counts = detail::block_intervals(z_mean, opt.n_stations);
        std::vector<detail::StationShape> shapes{hard.front()};
        std::vector<size_t> i_hard{0};
        for (size_t b = 0; b < counts.size(); b++)
        {
            // the last station of a block is the hard one itself: recomputing it can
            // land one ulp past it, i.e. outside the channel at the outlet
            for (size_t k = 1; k < counts[b]; k++)
                shapes.push_back(detail::blend(hard[b], hard[b + 1], double(k) / counts[b]));
            shapes.push_back(hard[b + 1]);
            i_hard.push_back(shapes.size() - 1);
        }
        const size_t ni = shapes.size();

        // Grid: the stations, stream lines evenly spread along them from hub to shroud
        MeridionalGrid<double> g(ni, nj);
        for (size_t i = 0; i < ni; i++)
        {
            const auto pts = detail::station_points(pb, shapes[i]);
            for (size_t j = 0; j < nj; j++)
            {
                g(i, j).x = pts[j][0];
                g(i, j).y = pts[j][1];
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

        // Boundary conditions, gas and solver settings
        // INLET_VmMoy_Ts_Ps_Vu sets the inlet totals from its Ts / Ps before the
        // inlet velocity is computed (Vm = 0): they are the total values, which
        // the solver keeps. (A static inlet would need its own formulation.)
        sc.inlet.mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
        sc.inlet.Vm_moy = pb.inlet.vm_mean;
        sc.inlet.Ts = pb.inlet.Tt;
        sc.inlet.Ps = pb.inlet.Pt;
        sc.inlet.Vu = pb.inlet.Vu;
        sc.gi->R = pb.gas.R;
        sc.gi->rho_cst = !opt.compressible;
        sc.gi->RF = opt.relax_factor;
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
        if (progress)
            sc.on_iteration = [&progress](size_t it, double r_max, double r_moy) { return progress(it, r_max, r_moy); };

        curvature_solver(sc);
        sc.on_iteration = nullptr; // do not keep a reference to the callback

        res.iterations = sc.log.delta_pos_max.size();
        res.residual_max = sc.log.delta_pos_max;
        res.residual_mean = sc.log.delta_pos_moy;

        const auto &g = *sc.gi->g;
        bool finite = std::all_of(g.begin(), g.end(), [](const auto &gp)
                                  { return std::isfinite(gp.Vm) && std::isfinite(gp.Ps) && std::isfinite(gp.Ts) && std::isfinite(gp.Vu); });
        if (!finite)
            res.status = Status::Diverged;
        else if (sc.cancelled)
            res.status = Status::Cancelled;
        else if (sc.converged)
            res.status = Status::Converged;
        else
            res.status = Status::NotConverged;
        return res;
    }
}
