#pragma once
#include <datastorage.h>
#include <optional>
#include <functional>
#include <gbs/curves>
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
        T RF = 0.01;
        T tol_newtow_mf_f = 1e-5;
        T tol_newtow_mf_u = 1e-6;
        size_t vm_distribution_max_count =100;
    };

    enum class MeridionalBladeMode
    {
        DESIGN_BETA_OUT,
        DESIGN_PSI,
        DIRECT
    };

    template <typename T>
    struct BladeInfo{
        std::string name;
        int i1 = -1;
        int is = -1; 
        int i2 = -1;
        T omg     = 0.;
        T omg_    = 0.;
        MeridionalBladeMode mode = MeridionalBladeMode::DIRECT;
        std::function<T(T)> beta_out;
        std::function<T(T)> psi;
    };

    enum class MeridionalBC
    {
        INLET_Mf_Ts_Ps_Vu,
        INLET_VmMoy_Ts_Ps_Vu,
    };

    template <typename T>
    struct InletBC{
        MeridionalBC mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
        std::function<T(T)> Ps = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Ts = [](auto l_rel){return 300.;};
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

    template <typename T>
    struct SolverCase
    {
        std::shared_ptr< GridInfo<T> > gi;
        std::vector<BladeInfo<T>> bld_info_lst;
        InletBC<T> inlet;
        std::vector<T> mf;
        SolverLog<T> log;
        size_t max_geom = 200;
        T eps = 0.00001;
        T tol_rel_mf = 1e-4;
        T tol_rel_pos = 1e-5;
        bool relocate = true;
    };

}