#pragma once
#include <string>
#include <datastorage.h>
#include <curvaturemethod.h>
#include <optional>
#include <functional>
#include <gbs/curves>
#include <algorithm>
#include <stdexcept>
#include <utility>
namespace yams
{

    template <typename T>
    struct GridInfo{
        std::shared_ptr< MeridionalGrid<T> > g;
        std::shared_ptr< Array2d<Grid2dMetricsPoint<T>> >   g_metrics;
        T d_ksi;
        T d_eth;
        size_t ni;
        size_t nj;
        T j_0 {0.5};
        const char* gas_name;
        bool rho_cst = true;
        T R = 287.04; // perfect gas constant
        T Pref = 1.01325e5;
        T Tref = 288.;
        /// Stream lines relocation relaxation: the fraction of each computed
        /// correction applied per iteration, at every station unless RF_local.
        T RF = 0.01;
        /// Relaxation at each station from its grid aspect ratio and Mach number
        /// (station_relaxation), in place of RF: the stability analysis of
        /// Wilkinson (1969) as applied to axisymmetric flow by Hearsey (1975).
        bool RF_local = false;
        /// Its constant: 8 from the analysis, 6 (Hearsey's default) for a margin.
        /// It is not set from curvature: use default_relax_constant(curvature,
        /// curvature_blend), and stay under max_relax_constant(curvature,
        /// curvature_blend) (6 is too large for ThreePoint).
        T RF_B = 6.;
        /// Its meridional Mach number squared is limited to this, keeping some
        /// relaxation in transonic flow (Hearsey's XMMAX).
        T RF_Mm2_max = 0.6;
        /// Stream lines' curvature (compute_curvature). The less stable the
        /// method, the smaller RF_B it needs: default_relax_constant.
        CurvatureMethod curvature = CurvatureMethod::DoubleDifferentiation;
        /// The share of ThreePoint in the Blended curvature, in [0, 1]: 0 is
        /// DoubleDifferentiation, 1 is ThreePoint. A small share already removes
        /// most of the zigzag, for few more iterations (NASA Rotor 37, 31
        /// stations: zigzag 0.61, 0.19, 0.12 and 126, 143, 225 iterations for 0,
        /// 0.2 and 0.5, each with its default_relax_constant).
        T curvature_blend = 0.5;
        /// The density follows the flow (Ps / (R Ts)) only once the mean
        /// displacement of the stream lines is under this, relative to the span.
        /// A run whose displacement stalls above it never updates its density
        /// and never converges (a flow near choke downstream of a rotor): the
        /// largest value (std::numeric_limits<T>::max()) updates it at every
        /// iteration, which also saves about a tenth of the iterations.
        T rho_update_residual = 1e-2;
        T tol_newtow_mf_f = 1e-5;
        T tol_newtow_mf_u = 1e-6;
        size_t vm_distribution_max_count =200;
        /// Points (i, j at i * nj + j) where, in the last integration of the
        /// meridional velocity across their station, the radial equilibrium
        /// stepped from the previous point gave no positive Vm (Vm^2 / 2 <= 0), the
        /// step holding it at its floor: a stream tube that stagnates or would flow
        /// back. The point's final Vm is positive but small, an artefact of that
        /// floor (5 to 19 m/s against 105 to 154 m/s on Turbo Blade Forge's demo at
        /// 50 and 35 m/s). Sized by curvature_solver; left empty, not filled.
        std::vector<unsigned char> vm_floor;
        std::function<T(T)> KD = [](T m_rel){return T{};};
    };

    enum class MeridionalBladeMode
    {
        DESIGN_BETA_OUT,
        DESIGN_ALPHA_OUT,
        DESIGN_PSI,
        DIRECT
    };

    template <typename T>
    struct BladeInfo{
        std::string name;
        int i1 = -1;
        int is = -1; 
        int i2 = -1;
        T omg  = 0.;
        T Ksh  = 1.;
        size_t z_{2};
        std::function<T(T)> omg_;
        std::function<T(T)> dev;
        std::function<T(T)> mu;
        MeridionalBladeMode mode = MeridionalBladeMode::DIRECT;
        std::function<T(T)> beta_out;
        std::function<T(T)> alpha_out;
        std::function<T(T)> psi;
        std::function<T(T,T)> k;
        std::function<T(T,T)> tb;
        std::function<T(T,T)> eps;
        std::function<T(T)> gauge;
        std::function<T(T)> max_cam_pos;
        std::function<T(T)> max_thickness;
        bool compute_dev{false};
    };

    enum class MeridionalBC
    {
        INLET_Mf_Tt_Pt_Vu,
        INLET_Mf_Ts_Ps_Vu,
        INLET_VmMoy_Ts_Ps_Vu,
        INLET_Vm_Ts_Ps_Vu,
        CON,
    };

    template <typename T>
    struct InletBC{
        MeridionalBC mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
        std::function<T(T)> Ps = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Ts = [](auto l_rel){return 300.;};
        std::function<T(T)> Vm = [](auto l_rel){return 30.;};
        std::function<T(T)> Vu = [](auto l_rel){return 0.;};
        std::function<T(T)> Pt = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Tt = [](auto l_rel){return 300.;};
        T Mf =1.;
        T Vm_moy = 30.;
    };

    template <typename T>
    struct SolverLog
    {
        std::vector<T> delta_pos_max;
        std::vector<T> delta_pos_moy;
        std::vector<std::vector<T>> delta_pos;
        void clear(){delta_pos.clear(),delta_pos_max.clear();delta_pos_moy.clear();}
    };

    // enum class MassFlowDistrib
    template <typename T>
    struct 
    SolverCase
    {
        std::shared_ptr< GridInfo<T> > gi;
        std::vector<BladeInfo<T>> bld_info_lst;
        InletBC<T> inlet;
        std::vector<T> mf;
        SolverLog<T> log;
        size_t max_geom = 200;
        T eps = T{1}/std::pow<T>(2,10);
        T tol_rel_mf = 1e-4;
        T tol_rel_pos = 1e-5;
        bool relocate = true;
        size_t mf_ref_span = 0;
        bool mf_uniform = false;
        bool use_meridional_grad=false;
        /// Called after each geometric iteration with (iteration, delta_pos_max, delta_pos_moy);
        /// returning false stops the computation (cancelled = true).
        std::function<bool(size_t, T, T)> on_iteration;
        /// Set by curvature_solver: stream lines positions converged within tol_rel_pos
        bool converged = false;
        /// Set by curvature_solver: stopped by on_iteration
        bool cancelled = false;
        /// Set by curvature_solver: the stream lines could no longer be relocated
        /// (they cross or coincide on a station, or their positions or residuals are
        /// not finite). The run ends there, on the grid of its last complete iteration.
        bool diverged = false;
        /// Why, when diverged
        std::string divergence;
        /// Shapes of the computing stations, one curve per station from hub to shroud.
        /// When given, the stream lines move along these curves, which the stations
        /// then keep exactly; station_params[i][j] is the parameter of the grid point
        /// (i, j) on its curve, updated by the solver. When empty, they move along the
        /// broken line through the current points of each station (a curved station
        /// then loses its shape as the stream lines move).
        std::vector<gbs::BSCurve<T, 2>> station_curves;
        std::vector<std::vector<T>> station_params;
    };

}