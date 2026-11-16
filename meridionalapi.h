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
 * V1 limits:
 * - the computing stations are lines of constant z from hub to shroud: the hub and
 *   shroud z must increase along the curves (axial or mixed flow channel, not radial);
 * - blade rows are in DIRECT mode (metal angle imposed, no deviation model).
 */
#include <curvaturesolver.h>
#include <gbs/bscurve.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace yams::api
{
    using Curve2d = gbs::BSCurve<double, 2>;

    /// Blade row, edges at constant z
    struct BladeRow
    {
        std::string name;
        double z_le{};          ///< leading edge axial position
        double z_te{};          ///< trailing edge axial position
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

        // Hard axial positions: inlet, rows edges, outlet
        double z_in = std::max(pb.hub(pb.hub.bounds()[0])[0], pb.shroud(pb.shroud.bounds()[0])[0]);
        double z_out = std::min(pb.hub(pb.hub.bounds()[1])[0], pb.shroud(pb.shroud.bounds()[1])[0]);
        std::vector<double> z_hard{z_in};
        for (const auto &row : pb.rows)
        {
            if (!row.metal_angle)
                throw std::invalid_argument("yams: blade row " + row.name + " has no metal angle");
            z_hard.push_back(row.z_le);
            z_hard.push_back(row.z_te);
        }
        z_hard.push_back(z_out);
        for (size_t k = 1; k < z_hard.size(); k++)
            if (!(z_hard[k] > z_hard[k - 1]))
                throw std::invalid_argument("yams: blade rows must be inside the channel, in axial order, with z_le < z_te and no overlap");

        // Stations
        auto counts = detail::block_intervals(z_hard, opt.n_stations);
        std::vector<double> z_stations{z_in};
        std::vector<size_t> i_hard{0};
        for (size_t b = 0; b < counts.size(); b++)
        {
            // the last station of a block is the hard position itself: recomputing it
            // can land one ulp past it, i.e. outside the channel at the outlet
            for (size_t k = 1; k < counts[b]; k++)
                z_stations.push_back(std::min(z_hard[b + 1], z_hard[b] + (z_hard[b + 1] - z_hard[b]) * k / counts[b]));
            z_stations.push_back(z_hard[b + 1]);
            i_hard.push_back(z_stations.size() - 1);
        }
        const size_t ni = z_stations.size();
        const size_t nj = opt.n_streams;

        // Grid: straight stations at constant z, stream lines evenly spread from hub to shroud
        MeridionalGrid<double> g(ni, nj);
        for (size_t i = 0; i < ni; i++)
        {
            auto h = pb.hub(detail::param_at_z(pb.hub, z_stations[i], "hub"));
            auto s = pb.shroud(detail::param_at_z(pb.shroud, z_stations[i], "shroud"));
            if (!(s[1] > h[1]))
                throw std::invalid_argument("yams: the shroud must be above the hub");
            for (size_t j = 0; j < nj; j++)
            {
                double t = double(j) / (nj - 1);
                g(i, j).x = z_stations[i];
                g(i, j).y = h[1] + t * (s[1] - h[1]);
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
