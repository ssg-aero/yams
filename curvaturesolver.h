#pragma once
#include <eqcurvaturesolver.h>
#include <meridionalsolvercase.h>
#include <gridmetrics.h>
#include <gbs/bscinterp.h>
#include <gbs/bscapprox.h>
#include <meshtools.h>
#include <parallel.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// const bool use_meridional_grad = false;
// const bool use_meridional_grad = true;
const bool verbose = false;
// const bool verbose = true;

namespace yams
{
    /// The stream lines of a station can no longer be relocated: they cross or
    /// coincide, or their positions are not finite. curvature_solver ends the run
    /// on it, with SolverCase::diverged; any other exception is an error and goes
    /// up to the caller.
    struct SolverDivergence : std::runtime_error
    {
        using std::runtime_error::runtime_error;
    };

    using gbs::operator*;
    using gbs::operator+;
    using gbs::operator-;

    const bool cap_vm_rk2 = false;

    inline auto f_sqV = [](const auto &gp)
    {
        return gp.Vm * gp.Vm + gp.Vu * gp.Vu;
    };

    inline auto f_sqW = [](const auto &gp)
    {
        auto Wu = gp.Vu - gp.omg * gp.y;
        return gp.Vm * gp.Vm + Wu * Wu;
    };

    inline auto f_Tt = [](const auto &gp)
    { 
        return gp.Ts + f_sqV(gp)  / 2. / gp.Cp; 
    };

    template <typename T>
    inline auto newton_solve = [](const auto &crv, auto p, T u0, T u1, T u2, T tol_f = 1.e-3, T tol_u = 1.e-4, size_t it_max = 100, T factor = 1.0)
    {
        auto delta = tol_u * 10., d0 = tol_f * 10.;
        auto u = u0;
        auto count = 0;
        while (fabs(delta) > tol_u && fabs(d0) > tol_f && count < it_max)
        {
            d0 = crv(u) - p;
            auto d1 = crv(u, 1);
            auto d2 = crv(u, 2);
            delta = d1 * d0 / (d2 * d0 + d1 * d1);
            u -= delta * factor;
            u = fmax(u1, fmin(u2, u));
            count++;
        }
        return std::make_tuple(u, delta, count);
    };

    template <typename T, typename _Func>
    void integrate_RK2_vm_sheet(size_t i,  GridInfo<T> &gi, _Func F, int j_beg, int j_end, int j_stp)
    {
        auto &g = *gi.g;
        auto &g_metrics= *gi.g_metrics;
        auto Vmi    = g(i, j_beg - j_stp).Vm;
        // auto sqVm_max_q2 = 3 * Vmi * Vmi / 2.;
        auto sqVm_max_q2 = 9. * Vmi * Vmi / 2.;
        // sqVm_max_q2 = sqVm_max_q2 * sqVm_max_q2 /2;
        auto nj = gi.nj;
        auto l = g(i,nj-1);
        auto Fjcap = 0.3*Vmi*Vmi/2;
        for (int j = j_beg; j != j_end; j+=j_stp)
        {

            auto Vm_prev =  g(i, j - j_stp).Vm;
            const auto &gp = g(i, j);
            const auto &gp_prev = g(i, j - j_stp);

            auto sqVmq2 = f_sqVmq2(gp_prev);
            auto dl = gp.l - gp_prev.l;

            auto Fjm= F(g,g_metrics, i, j - j_stp,gi.d_ksi,gi.d_eth);// assert(Fjm==Fjm);
            // Fjm = std::max<T>(-Fjcap/dl,std::min<T>(Fjcap/dl,Fjm));
            auto sqVmq2_1 = std::fmin(
                sqVm_max_q2, 
                std::fmax(1.0,sqVmq2 + Fjm * dl)
            );
            // the radial equilibrium gives no positive Vm from the previous point:
            // the stream tube stagnates (the step holds Vm at its floor)
            bool floored = sqVmq2 + Fjm * dl <= 0.;
            // auto sqVmq2_1 = std::fmax(0.1,sqVmq2 + Fjm * dl); 
            // auto sqVmq2_1 = sqVmq2 + Fjm * dl;
            g(i, j).Vm = sqrt(2. * sqVmq2_1);
            
            auto Fj = F(g,g_metrics, i, j,gi.d_ksi,gi.d_eth);// assert(Fj==Fj);
            // Fj = std::max<T>(-Fjcap/dl,std::min<T>(Fjcap/dl,Fj));
            auto sqVmq2_2 = std::fmin( 
                sqVm_max_q2, 
                std::fmax(1.0,sqVmq2 + Fj * dl)
            ); 
            floored = floored || sqVmq2 + Fj * dl <= 0.;
            // auto sqVmq2_2 = std::fmax(0.1,sqVmq2 + Fj * dl); 
            // auto sqVmq2_2 = sqVmq2 + Fj * dl; 
            // auto Vm_new = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
            // auto dVm = (Vm_new-Vm_prev)/Vm_prev;
            // g(i, j).Vm = Vm_prev + Vm_prev * cap(dVm,0.3);
            // auto dVm_dl = (Vm_new-Vm_prev)/dl;
            // g(i, j).Vm = Vm_prev + Vm_prev * cap(dVm,0.3)
            // g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));

            // g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
            // auto dVm = ( g(i, j).Vm - Vm_prev ) / Vmi;
            // dVm = cap(dVm,0.3);
            // g(i, j).Vm = dVm * Vmi + Vm_prev;

            g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));

            g(i, j).Vm = std::min(g(i, j).Vm, std::sqrt(g(i, j).ga * gi.R * g(i, j).Ts));
            if (gi.vm_floor.size() == g.size())
                gi.vm_floor[i * nj + j] = floored;

        }
    }

    template <typename T, typename _Func>
    void integrate_RK2_vm_sheet(T vmi, size_t i, GridInfo<T> &gi, _Func F, bool integrate)
    {
        auto &g = *gi.g;
        // this integration replaces the station's floor flags
        if (gi.vm_floor.size() == g.size())
            std::fill_n(gi.vm_floor.begin() + i * g.nCols(), g.nCols(), 0);

        if(!integrate)
        {
            std::for_each(g.begin(i),g.end(i),[vmi](auto &gp){gp.Vm=vmi;});
            return;
        }

        size_t nj = g.nCols();
        int j_0 = std::round((nj - 1 ) * gi.j_0);
        g(i, j_0).Vm = vmi;
        integrate_RK2_vm_sheet(i,gi,F,j_0+1,nj, 1);
        integrate_RK2_vm_sheet(i,gi,F,j_0-1,-1,-1);
    }

    /**
     * @brief Integrate mass-flow function over span (computation plane)
     * 
     * @tparam Iterator 
     * @param begin 
     * @param end 
     */
    template <typename Iterator>
    void compute_massflow_distribution(Iterator begin, Iterator end)
    {
        std::transform(
            begin,
            std::next(end, -1),
            std::next(begin),
            std::next(begin),
            [](const auto &gp_prev, auto &gp)
            {
                gp.q = (f_mf(gp_prev) + f_mf(gp)) * (gp.l - gp_prev.l) * 0.5 + gp_prev.q;
                return gp;
            });
    }

    template <typename T>
    auto compute_massflow(const GridInfo<T> &gi, int i)
    {
        const auto &g = *gi.g;
        auto nj = g.nCols();
        auto ni = g.nRows();
        auto mf = 0.;
        for (auto j = 1; j < nj; j++)
        {
            const auto &gp = g(i, j);
            const auto &gpp= g(i, j -1); 
            mf += (f_mf(gp) + f_mf(gpp)) * (gp.l - gpp.l) * 0.5;
        }
        size_t j_mean_line = std::round(nj / 2.);
        auto m_rel = g(i, j_mean_line).m / g(ni-1, j_mean_line).m;
        return mf * (1. - gi.KD(m_rel));
    }

    template <typename T>
    auto compute_gas_properties(SolverCase<T> &solver_case, int i)
    {
        auto &gi = *solver_case.gi;
        auto &g = *gi.g;
        auto nj = g.nCols();

        if (i != 0) // except inlet
        {
            // Compute totals values
            for (auto j = 0; j < nj; j++)
            {
                const auto &g1 = g(i - 1, j);
                auto &g2 = g(i, j);

                // Compute enthalpy variation
                g2.H = g1.H + g2.omg * (g2.y * g2.Vu - g1.y * g1.Vu); // So it also works if g1 is not a blade
                g2.I = g2.H - g2.omg * g2.y * g2.Vu;

                g2.Tt = g1.Tt + (g2.H - g1.H) / ( 0.5 * ( g1.Cp + g2.Cp) );
                auto ga = 0.5 * (g1.ga + g2.ga);
                if(g2.iB<0)
                {
                    auto P2is = g1.Pt * std::pow(g2.Tt / g1.Tt, ga / (ga - 1));
                    g2.Pt = P2is - 0.5 * g1.rho * g1.omg_ * f_sqW(g1);
                }
                else
                {
                    auto i1    = solver_case.bld_info_lst[g2.iB].i1;
                    auto i2    = solver_case.bld_info_lst[g2.iB].i2;
                    if(i>i1) // apply losses
                    {
                        T omg_{};
                        const auto &g_le = g(i1, j);
                        const auto &g_te = g(i2, j);
                        if (solver_case.bld_info_lst[g2.iB].omg_) // losses defined
                        {
                            omg_ = solver_case.bld_info_lst[g2.iB].omg_(g_le.l / g(i1, nj - 1).l);
                            omg_ *= (g2.m - g_le.m) / (g_te.m - g_le.m);
                        }
                        auto P2is = g_le.Pt * std::pow(g2.Tt / g_le.Tt, ga / (ga - 1));
                        g2.Pt = P2is - 0.5 * g_le.rho * omg_ * f_sqW(g_le);
                    }
                    else
                    {
                        auto P2is = g1.Pt * std::pow(g2.Tt / g1.Tt, ga / (ga - 1));
                        g2.Pt = P2is - 0.5 * g1.rho * g1.omg_ * f_sqW(g1);
                    }
                }
            }
        }

        // Compute static values
        for (auto j = 0; j < nj; j++)
        {
            auto &gp = g(i, j);
            // TODO use Coolprop
            gp.Ts = gp.Tt - f_sqV(gp) / 2. / gp.Cp;
            gp.Ps = gp.Pt * std::pow(gp.Ts / gp.Tt, gp.ga / (gp.ga - 1));
            // gp.rho = gp.rho + 0.01*(gp.Ps / gi.R / gp.Ts-gp.rho);
            // TODO update cp
            gp.ga = 1. / (1. - gi.R / gp.Cp);
            // Compute entropy rise
            // TODO compute elementary entropy rises from different losses and then compute Ps then Pt
            gp.s = gp.Cp * std::log(gp.Ts / gi.Tref) - gi.R * std::log(gp.Ps/gi.Pref);
        }
    }

    template<typename T>
    void eval_span_grad(SolverCase<T> &solver_case, int i)
    {
        auto &gi=*(solver_case.gi);
        auto &g_metrics = *gi.g_metrics;
        auto &g = *gi.g;
        auto nj = gi.nj;
        int j_0 = std::round((nj - 1 ) * gi.j_0);

        // for (size_t j{}; j < nj; j++)
        // {
        //     // g(i, j).dH_dl  = D1_O2_dx2(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_H); 
        //     // g(i, j).dI_dl  = D1_O2_dx2(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_I); 
        //     // g(i, j).dS_dl  = D1_O2_dx2(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_S_);
        //     // g(i, j).drVu_dl= D1_O2_dx2(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_rVu);
        //     // g(i, j).drTb_dl= D1_O2_dx2(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_rTanBeta);
        // }
        for (int j{1}; j < nj; j++)
        {
            g(i, j).dH_dl= ( g(i, j).H - g(i, j-1).H ) / ( g(i, j).l - g(i, j-1).l );
            g(i, j).dI_dl= ( g(i, j).I - g(i, j-1).I ) / ( g(i, j).l - g(i, j-1).l );
            g(i, j).dS_dl= ( g(i, j).s - g(i, j-1).s ) / ( g(i, j).l - g(i, j-1).l );
            g(i, j).drVu_dl= ( g(i, j).y * g(i, j).Vu - g(i, j-1).y * g(i, j-1).Vu ) / ( g(i, j).l - g(i, j-1).l );
            g(i, j).drTb_dl= ( g(i, j).y * tan(g(i, j).bet) - g(i, j-1).y * tan(g(i, j-1).bet) ) / ( g(i, j).l - g(i, j-1).l );
        }
        // for (int j{j_0-1}; j >= 0; j--)
        // {
        size_t j = 0 ; 
            g(i, j).dH_dl= ( g(i, j+1).H - g(i, j).H ) / ( g(i, j+1).l - g(i, j).l );
            g(i, j).dI_dl= ( g(i, j+1).I - g(i, j).I ) / ( g(i, j+1).l - g(i, j).l );
            g(i, j).dS_dl= ( g(i, j+1).s - g(i, j).s ) / ( g(i, j+1).l - g(i, j).l );
            g(i, j).drVu_dl= ( g(i, j+1).y * g(i, j+1).Vu - g(i, j).y * g(i, j).Vu ) / ( g(i, j+1).l - g(i, j).l );
            g(i, j).drTb_dl= ( g(i, j+1).y * tan(g(i, j+1).bet) - g(i, j).y * tan(g(i, j).bet) ) / ( g(i, j+1).l - g(i, j).l );
        // }
        // if(j_0==0)
        // {

        // }
        // else if
        // g(i, j0).drVu_dl = 0.5 * (  g(i, fmin(j0+1,nj-1)).drVu_dl + g(i, fmax(j0-1,0)).drVu_dl );
        // g(i, j0).drTb_dl = 0.5 * (  g(i, fmin(j0+1,nj-1)).drTb_dl + g(i, fmax(j0-1,0)).drTb_dl );
    }

    template <typename T>
    auto eq_massflow_no_blade(T vmi, SolverCase<T> &solver_case, int i, bool integrate) -> void
    {
        auto &gi=*(solver_case.gi);
        auto &g = *gi.g;
        auto nj = g.nCols();
        for (auto j = 0; j < nj; j++)
        {
            if (i > 0)
            {
                g(i, j).Vu = g(i, j).y > 0. ? g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y : 0.;
            }
            g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg , g(i, j).Vm); // <- lag from previous
        }

        compute_gas_properties(solver_case,i); // needed for span grad
        eval_span_grad(solver_case,i);
        integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
    }

    template <typename T>
    auto eq_massflow_blade_direct(T vmi, SolverCase<T> &solver_case, int i, int i1, int i2, bool integrate) -> void
    {
        auto &gi=*(solver_case.gi);
        auto &g = *gi.g;
        auto nj = g.nCols();
        if(i == i1)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).y > 0. ? g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y : 0.;
                g(i, j).bet = cap_angle( atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm)); // <- lag from previous
            }
            compute_gas_properties(solver_case,i); // needed for span grad
            eval_span_grad(solver_case,i);
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
        else{
            for (auto j = 0; j < nj; j++)
            {
                // apply deviation
                auto iB = g(i, j).iB;
                T dev{};
                if(iB>=0)
                {
                    auto bld_info = solver_case.bld_info_lst[iB];
                    if(bld_info.dev)
                    {
                        dev = bld_info.dev(g(i, j).l / g(i, nj - 1).l);
                        dev *= (g(i, j).m - g(i1, j).m) / (g(i2, j).m - g(i1, j).m);
                    }
                }
                if(g(i, j).omg>0.)
                {
                    dev = -dev;
                }
                g(i, j).bet = g(i, j).k + dev; 
                g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet) + g(i, j).y * g(i, j).omg;
            }
            compute_gas_properties(solver_case,i); // needed for span grad
            eval_span_grad(solver_case,i);
            integrate_RK2_vm_sheet(vmi, i, gi, eq_bet, integrate);
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet) + g(i, j).y * g(i, j).omg;
            }
        }
    }

    template <typename T>
    auto eq_massflow_blade_beta_out(T vmi, SolverCase<T> &solver_case, int i, int i1, int i2, const auto &bet_out, bool integrate) -> void
    {
        auto &gi=*(solver_case.gi);
        auto &g = *gi.g;
        auto nj = g.nCols();
        if(i == i1)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).y > 0. ? g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y : 0.;
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm); // <- lag from previous
            }
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
        else
        {
            for (auto j = 0; j < nj; j++)
            {
                auto m_rel_loc = (g(i, j).m - g(i1, j).m) / (g(i2, j   ).m - g(i1, j).m);
                auto l_rel     = (g(i, j).l - g(i , 0).l) / (g(i , nj-1).l - g(i , 0).l);
                auto bet_in  = g(i1, j).bet;
                g(i, j).bet = bet_in *(1.-m_rel_loc) + bet_out(l_rel) * m_rel_loc;
            }
            integrate_RK2_vm_sheet(vmi, i, gi, eq_bet, integrate);
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet) + g(i, j).y * g(i, j).omg; // <- lag from previous
            }
        }
    }

// Buggy
    template <typename T>
    auto eq_massflow_blade_alpha_out(T vmi, SolverCase<T> &solver_case, int i, int i1, int i2, const auto &alf_out, bool integrate) -> void
    {
        auto &gi=*(solver_case.gi);
        auto &g = *gi.g;
        auto nj = g.nCols();
        if(i == i1)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).y > 0. ? g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y : 0.;
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm); // <- lag from previous
            }
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
        else
        {
            for (auto j = 0; j < nj; j++)
            {
                auto m_rel_loc = (g(i, j).m - g(i1, j).m) / (g(i2, j   ).m - g(i1, j).m);
                auto l_rel     = (g(i, j).l - g(i , 0).l) / (g(i , nj-1).l - g(i , 0).l);
                auto alpha_in = atan2(g(i1, j).Vu, g(i1, j).Vm);
                auto alpha = alpha_in *(1.-m_rel_loc) + alf_out(l_rel) * m_rel_loc;
                g(i, j).Vu = g(i, j).Vm * tan(alpha);                                  // <- lag from previous
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm); // <- lag from previous
            }
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
    }

    template <typename T>
    auto eq_massflow_blade_design_psi(T vmi, SolverCase<T> &solver_case, int i, int i1, int i2, const auto &f_psi, bool integrate) -> void
    {
        auto &gi=*(solver_case.gi);
        auto &g = *gi.g;
        auto nj = g.nCols();

        if(i == i1)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).y > 0. ? g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y : 0.;
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm); // <- lag from previous
            }
            compute_gas_properties(solver_case,i); // needed for span grad
            eval_span_grad(solver_case,i);
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
        else
        {
            for (auto j = 0; j < nj; j++)
            {
                auto m_rel_loc = (g(i, j).m - g(i1, j).m) / (g(i2, j   ).m - g(i1, j).m);
                auto l_rel     = (g(i, j).l - g(i , 0).l) / (g(i , nj-1).l - g(i , 0).l);
                auto phi_out =f_psi(l_rel);
                g(i, j).Vu = g(i1, j).Vu + m_rel_loc * phi_out * g(i, j).y * g(i, j).omg;
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg , g(i, j).Vm); // <- lag from previous
            }
            compute_gas_properties(solver_case,i);// needed for span grad
            eval_span_grad(solver_case,i);
            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg , g(i, j).Vm); // <- lag from previous
            }
        }
    }
    /**
     * @brief 
     * 
     * @tparam T 
     * @param vmi : initial Vm value for integration, if no integration performed this value is the constant value on comptation plane
     * @param solver_case 
     * @param i : computation plane index
     * @param integrate : activate use contant value for Vm or not
     * @return T : Current massflow crossing computation plane
     */
    template <typename T>
    auto eq_massflow(T vmi, SolverCase<T> &solver_case, int i, bool integrate) -> T
    {
        auto &gi= *solver_case.gi;
        auto &g = *gi.g;
        auto nj = g.nCols();
        
        if (g(i, 0).iB == -1)
        {
            eq_massflow_no_blade(vmi, solver_case, i, integrate);
        }
        else
        {
            auto i1 = solver_case.bld_info_lst[g(i, 0).iB].i1;
            auto i2 = solver_case.bld_info_lst[g(i, 0).iB].i2;

            if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DIRECT)
            {
                eq_massflow_blade_direct(vmi,solver_case, i, i1, i2, integrate);
            }
            else if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DESIGN_BETA_OUT)
            {
                eq_massflow_blade_beta_out(vmi, solver_case, i, i1, i2, solver_case.bld_info_lst[g(i, 0).iB].beta_out, integrate);
            }
            else if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DESIGN_ALPHA_OUT)
            {
                eq_massflow_blade_alpha_out(vmi, solver_case, i, i1, i2, solver_case.bld_info_lst[g(i, 0).iB].alpha_out, integrate);
            }
            else if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DESIGN_PSI)
            {
                eq_massflow_blade_design_psi(vmi, solver_case, i, i1, i2, solver_case.bld_info_lst[g(i, 0).iB].psi, integrate);
            }
        }
        compute_gas_properties(solver_case,i);
        return compute_massflow(gi, i);
    }

    template <typename T>
    auto compute_span_curve(const MeridionalGrid<T> &g, int i)
    {
        auto nj = g.nCols();
        std::vector<T> u(nj);
        gbs::points_vector<T, 2> X(nj);
        auto l_tot = g(i, nj - 1).l;
        for (auto j = 0; j < nj; j++)
        {
            const auto &gp = g(i, j);
            X[j][0] = gp.x;
            X[j][1] = gp.y;
            u[j]    = gp.l / l_tot;
            // stream lines crossed, coincident or not finite: said here, before the
            // curve's construction refuses its parameters with its own exception
            if (j > 0 && !(u[j] > u[j - 1]))
                throw SolverDivergence("the stream lines cross or coincide on station " + std::to_string(i));
        }
        size_t p = fmax(fmin(3, nj), 1);
        return std::make_tuple( gbs::interpolate(X, u, 1), u );
        // return std::make_tuple( gbs::approx(X, u, p, true), u );
    }

    // template <typename T, typename _Func>
    // inline auto streamsheet_value_vector(const MeridionalGrid<T> &g, size_t i, _Func f)
    // {
    //     std::vector<T> vec(g.nCols());
    //     std::transform(
    //         GBS_PAR_EXEC
    //         g.begin(i),
    //         g.end(i),
    //         vec.begin(),
    //         [&f](const auto &gp) { return f(gp); });
    //         return vec;
    // }

    // template <typename T, typename _Func>
    // inline auto find_streamsheet_max(const MeridionalGrid<T> &g, size_t i, _Func f)
    // {
    //     auto vec = streamsheet_value_vector(g,i,f);
    //     return *std::max_element(
    //         GBS_PAR_EXEC
    //         vec.begin(),
    //         vec.end());
    // }

    // template <typename T, typename _Func>
    // inline auto find_streamsheet_max(const MeridionalGrid<T> &g, size_t i, size_t stride, _Func f)
    // {
    //     auto vec = streamsheet_value_vector(g,i,stride,f);
    //     return *std::max_element(
    //         GBS_PAR_EXEC
    //         vec.begin(),
    //         vec.end());
    // }

    /// Relaxation of the stream lines' relocation at station i: Wilkinson's
    /// optimum factor (Proc. IMechE 184 Pt 3G, 1969), in Hearsey's axisymmetric
    /// form (AD-A009 273, 1975),
    ///     1 / (1 + (1 - Mm^2) A^2 / B),
    /// A the aspect ratio of the grid at the station, its span over its spacing
    /// (station_min_spacing: the smallest to a neighbouring station along any
    /// stream line, so one place where two stations come close slows the whole
    /// station), Mm^2 the mean over the station of its meridional Mach number
    /// squared, limited to RF_Mm2_max and 0 for a constant density, B = RF_B.
    /// The correction of a stream line overshoots the more the stations are close
    /// compared with the span: a fixed factor that suits a coarse grid makes a
    /// fine one diverge. gi.RF when not RF_local.
    /// @throw std::invalid_argument unless RF_B > 0 and 0 <= RF_Mm2_max < 1
    /// @throw SolverDivergence if the station touches a neighbour (no spacing)
    template <typename T>
    T station_relaxation(const GridInfo<T> &gi, size_t i)
    {
        if (!gi.RF_local)
            return gi.RF;
        if (!(gi.RF_B > 0.) || !(gi.RF_Mm2_max >= 0. && gi.RF_Mm2_max < 1.))
            throw std::invalid_argument("yams: the local relaxation needs RF_B > 0 and 0 <= RF_Mm2_max < 1");
        const auto &g = *gi.g;
        const size_t nj = g.nCols();
        const T span = g(i, nj - 1).l;
        // a grid not yet measured (no abscissas): the fixed factor
        if (!(span > 0.))
            return gi.RF;
        T Mm2 = 0.;
        if (!gi.rho_cst) // a constant density has no Mach number effect
        {
            for (size_t j = 0; j < nj; j++)
            {
                const auto &gp = g(i, j);
                Mm2 += gp.Vm * gp.Vm / (gp.ga * gi.R * gp.Ts);
            }
            // not finite (no flow yet): none, the smallest factor
            Mm2 = std::isfinite(Mm2) ? std::clamp(Mm2 / nj, T{0}, gi.RF_Mm2_max) : T{0};
        }
        const T A = span / station_min_spacing(g, i);
        const T RF = 1. / (1. + (1. - Mm2) * A * A / gi.RF_B);
        // no spacing: the factor would be 0, and a fixed one the largest where the
        // smallest is needed
        if (!(RF > 0.))
            throw SolverDivergence("station " + std::to_string(i) + " touches a neighbouring station on a stream line, or its spacing is not finite");
        return std::min(RF, T{1});
    }

    template <typename T>
    auto balance_massflow(SolverCase<T> &solver_case, int i, T tol_mf)
    {
        auto i_ref = solver_case.mf_ref_span;
        if(!solver_case.mf_uniform && i==i_ref)
        {
            return T{};
        }
        auto &gi= *solver_case.gi;
        auto &g = *gi.g;
        auto nj = g.nCols();
        std::vector<T> q(nj);
        auto l_tot = g(i, nj - 1).l;
        for (auto j = 0; j < nj; j++)
        {
            // Use last cumulative mass flow rather than specified to perfectly match and then solve better
            q[j] = g(i, j).q * g(0, nj - 1).q / g(i, nj - 1).q; 
        }

        auto [f_X, u] = compute_span_curve(g,i);
        // the station's own curve, if the case gives it: the points stay on it
        const bool fixed_shape = !solver_case.station_curves.empty();
        if (fixed_shape)
            u = solver_case.station_params[i];
        // The points' parameters must increase from hub to shroud (they do not when
        // stream lines have crossed, coincide, or are not finite): said here, rather
        // than by the curve's construction refusing its knots.
        for (size_t j = 1; j < nj; j++)
            if (!(u[j] > u[j - 1]))
                throw SolverDivergence("the stream lines cross or coincide on station " + std::to_string(i));
        auto f_Q = gbs::interpolate(q, u, f_X.degree());

        auto span_geom_residual = 0.;
        auto RF = station_relaxation(gi, static_cast<size_t>(i));
        auto tol_f = gi.tol_newtow_mf_f;
        auto tol_u = gi.tol_newtow_mf_u;

        auto count = nj-1;
        auto [u1, u2] = f_Q.bounds();
        for (auto j = 1; j < count; j++)
        {
            auto q  = solver_case.mf_uniform ? (solver_case.inlet.Mf * j) / count : g(i_ref, j).q;
            auto [l, delta, count] = newton_solve<T>(f_Q, q, u[j], u1, u2, tol_f, tol_u);
            
            auto X = fixed_shape ? solver_case.station_curves[i].value(l) : f_X.value(l);
            // fmax below drops a NaN: a point that is not finite would count as not moved
            if (!std::isfinite(X[0]) || !std::isfinite(X[1]))
                throw SolverDivergence("a stream line's position is not finite on station " + std::to_string(i));
            auto dx = g(i, j).x - X[0];
            auto dy = g(i, j).y - X[1];
            if (fixed_shape)
            {
                // relax along the curve, not along the chord to the target
                auto &u_j = solver_case.station_params[i][j];
                u_j += RF * (l - u_j);
                auto P = solver_case.station_curves[i].value(u_j);
                g(i, j).x = P[0];
                g(i, j).y = P[1];
            }
            else
            {
                g(i, j).x += RF * (X[0] - g(i, j).x);
                g(i, j).y += RF * (X[1] - g(i, j).y);
            }
            
            span_geom_residual = fmax(fmax(fabs(dx), fabs(dy)), span_geom_residual);
        }
        if (fixed_shape)
        {
            // crossed by this relocation: refused now, the solver gives back the grid before it
            const auto &u_new = solver_case.station_params[i];
            for (size_t j = 1; j < nj; j++)
                if (!(u_new[j] > u_new[j - 1]))
                    throw SolverDivergence("the stream lines cross or coincide on station " + std::to_string(i));
        }
        return span_geom_residual;
    }

    template <typename T>
    auto compute_vm_distribution(SolverCase<T> &solver_case, T vmi, size_t i, T tol_rel_mf, T eps, bool integrate, T Vmi_max=360.)
    {
        auto mf     = solver_case.mf[i];
        auto err_mf = tol_rel_mf * 10.;
        auto mf_    = 0., mf_pre = 0.; // mf shall allways be strictly positive
        int count   = 0;
        auto max_count = solver_case.gi->vm_distribution_max_count;
        while (err_mf > tol_rel_mf && count < max_count)
        {
            if(i==0 && solver_case.inlet.mode == MeridionalBC::INLET_Mf_Tt_Pt_Vu)
            {
                apply_bc(solver_case);
            }
            mf_pre = eq_massflow(vmi - eps, solver_case, i, integrate);
            mf_ = eq_massflow(vmi, solver_case, i, integrate);
            vmi = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
            vmi = fmin(fmax(1e-3,vmi),Vmi_max); // TODO improve test
            err_mf = fabs(mf_ - mf) / mf;
            count++;
        }
        if (verbose)
        {
            if (count == max_count && err_mf > tol_rel_mf)
            {
                std::cout << "Warning span: " << i << " did not converged after " << count << " err_mf: " << err_mf * 100 << "%" << "vmi:" << vmi <<"m/s" <<  std::endl;
            }
            if (integrate)
            {
                std::cout << " i: " << i << " count: " << count << " err_mf: " << err_mf << std::endl;
            }
        }
    }

    template <typename T>
    auto compute_vm_distribution(SolverCase<T> &solver_case, T tol_rel_mf, T eps, bool integrate, size_t i0)
    {
        auto &gi = *solver_case.gi;
        auto &g = *gi.g;
        // auto &g_metrics = *gi.g_metrics;
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        T vmi{};
        auto j0 = std::round((nj - 1 ) * gi.j_0);
        for (auto i = i0; i < ni; i++)
        {
            if( ( solver_case.inlet.mode != MeridionalBC::INLET_Vm_Ts_Ps_Vu && solver_case.inlet.mode != MeridionalBC::CON )
                    || i != 0 )
            {
                vmi = g(i, j0).Vm;
                auto aj  = std::sqrt( g(i, j0).ga * gi.R * g(i, j0).Ts );
                compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps,integrate, aj);
            }
            // compute_gas_properties(solver_case,i);
        }
    }

    template <typename T>
    auto apply_bc(SolverCase<T> &solver_case)
    {
        auto &gi = *solver_case.gi;
        auto &g  = *gi.g;
        const auto &inlet = solver_case.inlet;
        size_t nj = g.nCols();
        auto l_tot = g(0, nj - 1).l;
        auto Pref = gi.Pref;
        auto Tref = gi.Tref;
        if (solver_case.inlet.mode == MeridionalBC::INLET_VmMoy_Ts_Ps_Vu ||
            solver_case.inlet.mode == MeridionalBC::INLET_Mf_Ts_Ps_Vu)
        {
            std::for_each(g.begin(0), g.end(0), [l_tot, Tref, Pref, &inlet, &gi](auto &gp)
                          {
                              auto l_rel = gp.l / l_tot;
                              gp.Vu = inlet.Vu(l_rel);
                              gp.Ts = inlet.Ts(l_rel);
                              gp.Ps = inlet.Ps(l_rel);
                              if(!gi.rho_cst)
                                    gp.rho = gp.Ps / (gi.R) / gp.Ts;
                              gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / gp.Cp;
                            //   gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
                                gp.Pt = gp.Ps / std::pow(gp.Ts / gp.Tt, gp.ga / (gp.ga - 1));
                              gp.H = gp.Tt * gp.Cp;
                              gp.s = std::log(pow(gp.Ts / Tref, gp.Cp) / std::pow(gp.Ps / Pref, gi.R));
                                                      });
        }
        if(solver_case.inlet.mode == MeridionalBC::INLET_Vm_Ts_Ps_Vu)
        {
           std::for_each(g.begin(0), g.end(0), [l_tot, Tref, Pref, &inlet, &gi](auto &gp)
                          {
                            auto l_rel = gp.l / l_tot;
                            gp.Vm = inlet.Vm(l_rel);
                            gp.Vu = inlet.Vu(l_rel);
                            gp.Ts = inlet.Ts(l_rel);
                            gp.Ps = inlet.Ps(l_rel);
                            if(!gi.rho_cst)
                                gp.rho = gp.Ps / (gi.R) / gp.Ts;
                            gp.Tt = gp.Ts + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / gp.Cp;
                            //   gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
                            gp.Pt = gp.Ps / std::pow(gp.Ts / gp.Tt, gp.ga / (gp.ga - 1));
                            gp.H = gp.Tt * gp.Cp;
                            gp.s = std::log(pow(gp.Ts / Tref, gp.Cp) / std::pow(gp.Ps / Pref, gi.R));
                          });
        }
        if(solver_case.inlet.mode == MeridionalBC::INLET_Mf_Tt_Pt_Vu)
        {
           std::for_each(g.begin(0), g.end(0), [l_tot, Tref, Pref, &inlet, &gi](auto &gp)
                          {
                            auto l_rel = gp.l / l_tot;
                            gp.Vm = inlet.Vm(l_rel);
                            gp.Vu = inlet.Vu(l_rel);
                            gp.Tt = inlet.Tt(l_rel);
                            gp.Pt = inlet.Pt(l_rel);
                            gp.Ts = gp.Tt - (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. / gp.Cp;
                            gp.Ps = gp.Pt / std::pow(gp.Tt / gp.Ts, gp.ga / (gp.ga - 1));
                            if(!gi.rho_cst)
                                gp.rho = gp.Ps / (gi.R) / gp.Ts;
                            gp.H = gp.Tt * gp.Cp;
                            gp.s = std::log(pow(gp.Ts / Tref, gp.Cp) / std::pow(gp.Ps / Pref, gi.R));
                          });
        }
        compute_gas_properties(solver_case,0);
    }

    template <typename T>
    auto apply_mf(SolverCase<T> &solver_case)
    {
        auto &gi = *solver_case.gi;
        auto &g  = *gi.g;
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        if(solver_case.inlet.mode == MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
        {
            std::for_each(g.begin(0),g.end(0),
                [Vm = solver_case.inlet.Vm_moy](auto &gp)
                {
                    gp.Vm=Vm;
                }
            );
            solver_case.inlet.Mf = compute_massflow(gi, 0);
            // std::cout << "Mass flow set to: " << solver_case.inlet.Mf <<std::endl;
        }
        if(solver_case.inlet.mode == MeridionalBC::INLET_Vm_Ts_Ps_Vu)
        {
            std::for_each(g.begin(0),g.end(0),
                [l_tot = g(0,nj-1).l,&solver_case](auto &gp)
                {
                    gp.Vm=solver_case.inlet.Vm(gp.l / l_tot );
                }
            );
            solver_case.inlet.Mf = compute_massflow(gi, 0);
            // std::cout << "Mass flow set to: " << solver_case.inlet.Mf <<std::endl;
        }
        if(solver_case.inlet.mode == MeridionalBC::CON)
        {
            solver_case.inlet.Mf = compute_massflow(gi, 0);
        }
        solver_case.mf.resize(ni);
        std::fill(solver_case.mf.begin(),solver_case.mf.end(),solver_case.inlet.Mf); // Todo add leakage and reintroduction
    }  

    template <typename T>
    auto init_values(SolverCase<T> &solver_case, T tol_rel_mf, T eps)
    {
        compute_vm_distribution(solver_case, tol_rel_mf, eps, false, 0);
        compute_vm_distribution(solver_case, tol_rel_mf, eps, false, 0);
        auto ni = solver_case.gi->ni;
        auto nj = solver_case.gi->nj;
        auto &gi   = *solver_case.gi;
        auto &g    = *gi.g;
        
        // auto span_range = gbs::make_range<size_t>(0,ni-1);
        // // update density
        // if( !gi.rho_cst  )
        // {
        //     std::for_each(
        //             // GBS_PAR_EXEC
        //             span_range.begin(), span_range.end(),
        //             [&](const auto &i){
        //                 for (size_t j{}; j < nj; j++)
        //                 {
        //                     auto rho = g(i, j).Ps / gi.R / g(i,j).Ts;
        //                     g(i, j).rho = rho;
        //                     // std::cout<<"i: " << i << " rho: " << rho << std::endl;
        //                 }
        //             }
        //     );
        // }
    }

    template <typename T>
    auto apply_rotation_speeds(SolverCase<T> &solver_case)
    {
        auto &gi   = *solver_case.gi;
        auto &g    = *gi.g;
        size_t nj = g.nCols();
        for (const auto &bld_info : solver_case.bld_info_lst)
        {
            auto i1 = bld_info.i1;
            auto i2 = bld_info.i2;
            auto omg= bld_info.omg;
            for (auto i = i1; i <= i2; i++)
            {
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).omg = omg;
                }
            }
        }
    }

    template <typename T>
    auto curvature_solver(SolverCase<T> &solver_case)
    {
        // init values
        auto &gi = *solver_case.gi;
        auto &g = *gi.g;
        auto &g_metrics = *gi.g_metrics;
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        if (ni < 3 || nj < 3)
        {
            throw std::length_error("Grid must have dimensions >= 3");
        }
        if (!solver_case.station_curves.empty())
        {
            const auto &params = solver_case.station_params;
            if (solver_case.station_curves.size() != ni || params.size() != ni ||
                std::any_of(params.begin(), params.end(), [nj](const auto &p) { return p.size() != nj; }))
                throw std::invalid_argument("yams: station_curves and station_params must give one curve and nj parameters per station");
        }
        size_t max_geom=solver_case.max_geom;
        auto eps = solver_case.eps;
        auto tol_rel_mf =solver_case.tol_rel_mf;
        auto tol_pos = solver_case.tol_rel_pos; // delta_pos_array is already relative to the span length
        T vmi;
        int count_geom = 0;
        auto converged = false;
        solver_case.converged = false;
        solver_case.cancelled = false;
        solver_case.diverged = false;
        solver_case.divergence.clear();
        // auto i_0 = 0;
        // solver_case.log.clear();
        T delta_pos_max {};
        T delta_pos {};
        T delta_pos_moy {};
        std::vector<T> delta_pos_array(ni);
        // The stream lines' positions before their relocation, to give back the
        // grid of the last complete iteration when one fails
        std::vector<std::array<T, 2>> saved_xy(ni * nj);
        std::vector<std::vector<T>> saved_params;
        bool relocating = false;
        // the case's curvature method from the first iteration: make_grid_info
        // computes the double differentiation, which is left as the caller gives it
        if (gi.curvature != CurvatureMethod::DoubleDifferentiation)
            compute_curvature(g, true, gi.curvature, gi.curvature_blend);
        // apply boundary conditions
        apply_bc(solver_case);
        // compute spans mass flow
        apply_mf(solver_case);
        // apply blade info
        apply_blade_info(solver_case);
        // innit values
        gi.vm_floor.assign(g.size(), 0);
        init_values(solver_case,tol_rel_mf, eps);
        // apply rotation speeds
        apply_rotation_speeds(solver_case);
        // run computation
        while (!converged && (count_geom < max_geom))
        {
          // An iteration whose stream lines can no longer be relocated (they cross,
          // or their positions or residuals are not finite) ends the run, diverged,
          // on the grid of the last complete iteration. Any other exception is an
          // error and goes up to the caller.
          relocating = false;
          try
          {
            // integrate radial eq equation and update gas properties
            size_t grad_count{};
            size_t max_grd_count = solver_case.use_meridional_grad ? 10 : 1;
            while(grad_count < max_grd_count) // TODO add criteria
            {
                compute_vm_distribution(solver_case, tol_rel_mf, eps, true, 0 );
                // compute_vm_distribution2(solver_case, tol_rel_mf, eps, true);
                if(solver_case.use_meridional_grad)
                {
                    detail::par_for(ni, // update meridional gradients
                        [&](const auto &i){
                            for (size_t j{}; j < nj; j++)
                            {
                                // if(i==0 || i == ni-1)
                                // {
                                //     g(i, j).dsqVm_dm_2 = 0.;
                                //     g(i, j).ds_dm      = 0.;
                                // }
                                // else
                                // {
                                //     g(i, j).dsqVm_dm_2 = D1_O2_dx1(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_sqVmq2);
                                //     g(i, j).ds_dm      = D1_O2_dx1(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_S_);
                                //     g(i,j).drtb_dm     = D1_O2_dx1(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_rTanBeta);
                                //     g(i,j).drVu_dm     = D1_O2_dx1(g, g_metrics, i, j, gi.d_ksi, gi.d_eth, f_rVu);
                                // }
                                // if(i==0)
                                // {
                                //     // g(i, j).dsqVm_dm_2 = 0.;
                                //     g(i, j).dsqVm_dm_2 = ( f_sqVmq2(g(i+1 , j)) - f_sqVmq2(g(i, j)) / (g(i+1, j).m -  g(i, j).m) );
                                //     g(i, j).ds_dm      = 0.;
                                // }

                                // else if( i == ni-1)
                                // {
                                //     g(i, j).dsqVm_dm_2 = ( f_sqVmq2(g(i , j)) - f_sqVmq2(g(i-1, j)) / (g(i, j).m -  g(i-1, j).m) );
                                // }
                                if(i==0 || i == ni-1)
                                {
                                    g(i, j).dsqVm_dm_2 = 0.;
                                    g(i, j).ds_dm      = 0.;
                                    g(i,j).drtb_dm     = 0.;
                                    g(i,j).drVu_dm     = 0.;
                                }
                                else
                                {
                                    auto dm1 =  g(i+1, j).m -  g(i, j).m;
                                    auto dm2 =  g(i, j).m -  g(i-1, j).m;
                                    const auto &gp1 = g(i+1, j);
                                    const auto &gp2 = g(i , j);
                                    const auto &gp3 = g(i-1, j);
                                    g(i, j).dsqVm_dm_2 = ( ( f_sqVmq2(gp1) - f_sqVmq2(gp2) ) / dm1 + ( f_sqVmq2(gp2) - f_sqVmq2(gp3) )/dm2 )/2;
                                    g(i, j).ds_dm      = ( ( f_S_(gp1) - f_S_(gp2) ) / dm1 + ( f_S_(gp2) - f_S_(gp3) )/dm2 )/2;
                                    g(i,j).drtb_dm     = ( ( f_rTanBeta(gp1) - f_rTanBeta(gp2) ) / dm1 + ( f_rTanBeta(gp2) - f_rTanBeta(gp3) )/dm2 )/2;
                                    g(i,j).drVu_dm     = ( ( f_rVu(gp1) - f_rVu(gp2) ) / dm1 + ( f_rVu(gp2) - f_rVu(gp3) )/dm2 )/2;
                                }
                            }
                        }
                    );
                }
                grad_count++;
            }
            // compute_vm_distribution2(solver_case, tol_rel_mf, eps, true);
            if( !solver_case.relocate )
            {
                break;
            }
            // update density
            if(!gi.rho_cst && delta_pos_moy < gi.rho_update_residual )
            {
                detail::par_for(ni,
                        [&](const auto &i){
                            for (size_t j{}; j < nj; j++)
                            {
                                g(i, j).rho = g(i, j).Ps / gi.R / g(i,j).Ts;
                                // g(i, j).rho = g(i, j).rho + 0.01 * ( g(i, j).Ps / gi.R / g(i,j).Ts - g(i, j).rho );
                            }
                        }
                );
            }
            // Compute mass flow distribution
           detail::par_for(ni,
                [&](const auto &i){
                    compute_massflow_distribution(g.begin(i), g.end(i));
                    // for (size_t j{}; j < nj; j++)
                    // {
                    //     g(i, j).Vm_pre = g(i, j).Vm;
                    // }
                }
           );
           // relocate streams to balance mass flow
           for (size_t i{}; i < ni; i++)
               for (size_t j{}; j < nj; j++)
                   saved_xy[i * nj + j] = {g(i, j).x, g(i, j).y};
           saved_params = solver_case.station_params;
           relocating = true;
           detail::par_for(ni,
                [&](const auto &i)
                {
                    delta_pos_array[i] = balance_massflow(solver_case, i, tol_rel_mf * solver_case.mf[i]) / g(i,nj-1).l;
                }
           );
            // // apply boundary conditions if static conditions
            // if(solver_case.inlet.mode != MeridionalBC::INLET_Mf_Tt_Pt_Vu)
            //     apply_bc(solver_case);
            // compute residuals
            delta_pos_moy = std::reduce(
                    GBS_PAR_EXEC
                    delta_pos_array.begin(),delta_pos_array.end()
            ) / delta_pos_array.size();

            delta_pos_max = *std::max_element(
                delta_pos_array.begin(),delta_pos_array.end()
            );
            if (!std::isfinite(delta_pos_moy) || !std::isfinite(delta_pos_max))
                throw SolverDivergence("the stream lines' displacement is not finite");
            relocating = false;
           // update blades info
           apply_blade_info(solver_case);
            // add residuals to logger
            solver_case.log.delta_pos.push_back(delta_pos_array);
            solver_case.log.delta_pos_max.push_back(delta_pos_max);
            solver_case.log.delta_pos_moy.push_back(delta_pos_moy);
            // update metrics
            compute_grid_metrics(g,g_metrics,f_m,f_l,gi.curvature,gi.curvature_blend);// TODO run in // 
            // update convergence criteria
            converged = delta_pos_moy < tol_pos;
            count_geom++;
          }
          catch (const SolverDivergence &e)
          {
              if (relocating)
              {
                  // some stations have moved, others not: back to the last grid
                  for (size_t i{}; i < ni; i++)
                      for (size_t j{}; j < nj; j++)
                      {
                          g(i, j).x = saved_xy[i * nj + j][0];
                          g(i, j).y = saved_xy[i * nj + j][1];
                      }
                  solver_case.station_params = saved_params;
              }
              solver_case.diverged = true;
              solver_case.divergence = e.what();
              break;
          }
            if (solver_case.on_iteration && !solver_case.on_iteration(count_geom, delta_pos_max, delta_pos_moy))
            {
                solver_case.cancelled = true;
                break;
            }
        }
        solver_case.converged = converged;
    }

}
