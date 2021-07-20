#pragma once
#include <eqcurvaturesolver.h>
#include <datastorage.h>
#include <gridmetrics.h>
#include <gbs/bscinterp.h>
#include <optional>
#include <functional>

namespace quiss
{

    using gbs::operator*;
    using gbs::operator+;
    using gbs::operator-;

    const bool cap_vm_rk2 = false;

    auto f_sqV = [](const auto &gp)
    {
        return gp.Vm * gp.Vm + gp.Vu * gp.Vu;
    };

    auto f_sqW = [](const auto &gp)
    {
        auto Wu = gp.Vu - gp.omg * gp.y;
        return gp.Vm * gp.Vm + Wu * Wu;
    };

    auto f_Tt = [](const auto &gp)
    { 
        return gp.Ts + f_sqV(gp)  / 2. / gp.Cp; 
    };

    template <typename T>
    struct GridInfo{
        MeridionalGrid<T> &g;
        Array2d<Grid2dMetricsPoint<T>>   &g_metrics;
        T d_ksi;
        T d_eth;
        size_t ni;
        size_t nj;
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
        DESIGN_PHI,
        DIRECT
    };
    template <typename T>
    struct BladeInfo{
        std::string name;
        size_t i1 = -1;
        size_t i2 = -1;
        T omg     = 0.;
        T omg_    = 0.;
        MeridionalBladeMode mode = MeridionalBladeMode::DIRECT;
        std::optional< gbs::BSCfunction<T> > beta_out;
        std::optional< gbs::BSCfunction<T> > phi;
    };

    enum class MeridionalBC
    {
        INLET_Mf_Ts_Ps_Vu,
        INLET_VmMoy_Ts_Ps_Vu,
    };

    template <typename T>
    struct Inlet_BC{
        MeridionalBC mode = MeridionalBC::INLET_VmMoy_Ts_Ps_Vu;
        std::function<T(T)> Ps = [](auto l_rel){return 1.01325e5;};
        std::function<T(T)> Ts = [](auto l_rel){return 300.;};
        std::function<T(T)> Vu = [](auto l_rel){return 0.;};
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
        GridInfo<T> &gi;
        std::vector<BladeInfo<T>> bld_info_lst;
        Inlet_BC<T> inlet;
        std::vector<T> mf;
        SolverLog<T> log;
        size_t max_geom = 200;
        T eps = 0.00001;
        T tol_rel_mf = 1e-4;
        T tol_rel_pos = 1e-5;
    };

    template <typename T>
    auto newton_solve = [](const auto &crv, auto p, T u0, T u1, T u2, T tol_f = 1.e-3, T tol_u = 1.e-4, size_t it_max = 100, T factor = 1.0)
    {
        auto delta = tol_u * 10., d0 = tol_f * 10.;
        auto u = u0;
        auto count = 0;
        while (fabs(delta) > tol_u && fabs(d0) > tol_f && count < it_max)
        {
            auto d0 = crv(u) - p;
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
    void integrate_RK2_vm_sheet(T vmi, size_t i, GridInfo<T> &gi, _Func F)
    {
        auto &g = gi.g;
        auto &g_metrics= gi.g_metrics;

        g(i, 0).Vm = vmi;
        size_t nj = g.nCols();

        // eval_H_s(gi,i,0);

        // for (auto j = 1; j < nj; j++)
        // {
        //     const auto &gp = g(i, j);
        //     const auto &gp_prev = g(i, j - 1);

        //     auto sqVmq2 = f_sqVmq2(gp_prev);
        //     auto dl = gp.l - gp_prev.l;
        //     auto Fjm= F(g,g_metrics, i, j - 1,gi.d_ksi,gi.d_eth);
        //     assert(Fjm==Fjm);

        //     auto sqVmq2_1 = std::fmax(0.1,sqVmq2 + Fjm * dl);
        //     // g(i, j).Vm = std::fmin(sqrt(2. * sqVmq2_1),320.);// TODO add H and s to eq
        //     g(i, j).Vm = sqrt(2. * sqVmq2_1);
        //     // eval_H_s(gi,i,j); // equations are using H and s (enthalpy and entropy)
        //     auto Fj = F(g,g_metrics, i, j,gi.d_ksi,gi.d_eth);
        //     assert(Fj==Fj);

        //     auto sqVmq2_2 = std::fmax(0.1,sqVmq2 + Fj * dl); 
        //     // g(i, j).Vm = std::fmin(0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2)),320.);
        //     g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
        //     // eval_H_s(gi,i,j);
        // }

        // int j_beg = nj-2;
        // int j_end = -1;
        // int j_stp = -1;
        int j_beg = 1;
        int j_end = nj;
        int j_stp = 1;
        for (int j = j_beg; j != j_end; j+=j_stp)
        {
            const auto &gp = g(i, j);
            const auto &gp_prev = g(i, j - j_stp);

            auto sqVmq2 = f_sqVmq2(gp_prev);
            auto dl = gp.l - gp_prev.l;

            auto Fjm= F(g,g_metrics, i, j - j_stp,gi.d_ksi,gi.d_eth); assert(Fjm==Fjm);
            auto sqVmq2_1 = std::fmax(0.1,sqVmq2 + Fjm * dl);
            g(i, j).Vm = sqrt(2. * sqVmq2_1);
            
            auto Fj = F(g,g_metrics, i, j,gi.d_ksi,gi.d_eth); assert(Fj==Fj);
            auto sqVmq2_2 = std::fmax(0.1,sqVmq2 + Fj * dl); 
            g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
        }

    }


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
    auto compute_massflow(const MeridionalGrid<T> &g, int i)
    {
        auto nj = g.nCols();
        auto mf = 0.;
        for (auto j = 1; j < nj; j++)
        {
            mf += (f_mf(g(i, j)) + f_mf(g(i, j - 1))) * (g(i, j).l - g(i, j - 1).l) * 0.5;
        }
        return mf;
    }

    template <typename T>
    auto compute_gas_properties(GridInfo<T> &gi, int i)
    {
        auto &g = gi.g;
        auto nj = g.nCols();

        if (i != 0) // except inlet
        {
            for (auto j = 0; j < nj; j++)
            {
                const auto &g1 = g(i - 1, j);
                auto &g2 = g(i, j);

                // Compute enthalpy variation
                g2.H = g1.H + g2.omg * (g2.y * g2.Vu - g1.y * g1.Vu); // So it also works if g1 is not a blade
                g2.I = g2.H - g2.omg * g2.y * g2.Vu;

                g2.Tt = g1.Tt + (g2.H - g1.H) / ( 0.5 * ( g1.Cp + g2.Cp) );
                auto ga = 0.5 * (g1.ga + g2.ga);
                auto P2is = g1.Pt * std::pow(g2.Tt / g1.Tt, ga / (ga - 1));
                g2.Pt = P2is - 0.5 * g2.rho * g2.omg_ * f_sqW(g2); // using rho from previous
            }
        }

        for (auto j = 0; j < nj; j++)
        {
            auto &gp = g(i, j);
            gp.Ts = gp.Tt - f_sqV(gp) / 2. / gp.Cp;
            gp.Ps = gp.Pt * std::pow(gp.Ts / gp.Tt, gp.ga / (gp.ga - 1));
            // TODO use Coolprop
            if(!gi.rho_cst)
            {
                gp.rho = gp.Ps / gi.R / gp.Ts;
            }
            // TODO update cp
            gp.ga = 1. / (1. - gi.R / gp.Cp);
            // Compute entropy rise
            // TODO compute elemetary entropy rises from different losses and the compute Ps then Pt
            gp.s = std::log(pow(gp.Ts/gi.Tref,gp.Cp)/std::pow(gp.Ps/gi.Pref,gi.R)) ;
        }
    }

    template <typename T>
    auto eq_massflow(T vmi, GridInfo<T> &gi, int i)
    {
        auto &g = gi.g;
        auto nj = g.nCols();
        if (g(i, 0).iB == -1)
        {
            for (auto j = 0; j < nj; j++)
            {
                if (i > 0)
                {
                    g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
                }
                g(i, j).bet = atan2(g(i, j).Vu, g(i, j).Vm);
            }

            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu);
        }
        else
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).bet = g(i, j).k;
            }
            integrate_RK2_vm_sheet(vmi, i, gi, eq_bet);
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet + g(i, j).y * g(i, j).omg );
            }
        }
        // compute_gas_properties(gi,i);
        return compute_massflow(g, i);
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
        }
        size_t p = fmax(fmin(3, nj), 1);
        return std::make_tuple( gbs::interpolate(X, u, p), u );
    }

    template <typename T>
    auto balance_massflow(GridInfo<T> &gi, int i, T tol_mf)
    {
        auto &g = gi.g;
        auto nj = g.nCols();
        std::vector<T> q(nj);
        // gbs::points_vector<T,1> q(nj);
        gbs::points_vector<T, 2> X(nj);
        auto l_tot = g(i, nj - 1).l;
        for (auto j = 0; j < nj; j++)
        {
            q[j] = g(i, j).q * g(0, nj - 1).q / g(i, nj - 1).q; // To perfectly match and then solve better
        }
        // size_t p = fmax(fmin(3, nj), 1);
        auto [f_X, u] = compute_span_curve(g,i);
        auto f_Q = gbs::interpolate(q, u, f_X.degree());

        // auto f_Q = gbs::BSCfunction( gbs::approx(q,2,fmax(nj / 3,3), u,true) );
        // auto f_X =  gbs::approx(X,2,fmax(nj / 3,3), u,true);

        auto delta_pos = 0.;
        auto RF = gi.RF;
        auto tol_f = gi.tol_newtow_mf_f;
        auto tol_u = gi.tol_newtow_mf_u;

        // T B_ = fmin(8., 0.5 * 60. / g.nRows());
        // auto RF = eval_RF(g,i,B_);

        for (auto j = 1; j < nj - 1; j++)
        {
            // auto l = newton_solve<T>(f_Q, gbs::point<T,1>{g(0, j).q}, u[j]);
            auto [u1, u2] = f_Q.bounds();
            auto [l, delta, count] = newton_solve<T>(f_Q, g(0, j).q, u[j], u1, u2, tol_f, tol_u);
            assert(l <= u2 && l >= u1);
            auto X = f_X.value(l);
            auto dx = g(i, j).x - X[0];
            auto dy = g(i, j).y - X[1];
            // std::cout << delta << " " << dx << " " << dy  << std::endl;
            delta_pos = fmax(fmax(fabs(dx), fabs(dy)), delta_pos);
            g(i, j).x += RF * (X[0] - g(i, j).x);
            g(i, j).y += RF * (X[1] - g(i, j).y);
        }
        return delta_pos;
    }

    template <typename T>
    auto compute_vm_distribution(T mf, T vmi, size_t i,GridInfo<T> &gi, T tol_rel_mf, T eps)
    {
        auto err_mf = tol_rel_mf * 10.;
        auto mf_ = 0., mf_pre = 0.; // mf shall allways be strictly positive
        int count = 0;
        auto max_count = gi.vm_distribution_max_count;
        while (err_mf > tol_rel_mf && count < max_count)
        {
            mf_pre = eq_massflow(vmi - eps, gi, i);
            mf_ = eq_massflow(vmi, gi, i);
            vmi = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
            assert(vmi >= 0.);
            vmi = fmin(fmax(0.1,vmi),360.);
            err_mf = fabs(mf_ - mf) / mf;
            count++;
        }
        if(count==max_count && err_mf <= tol_rel_mf)
        {
            std::cout << "Warning span: " << i << " did not converged." << std::endl;
        }
    }

    template <typename T>
    auto apply_bc(quiss::SolverCase<T> &solver_case)
    {
        auto &gi = solver_case.gi;
        auto &g  = gi.g;
        const auto &inlet = solver_case.inlet;
        size_t nj = g.nCols();
        auto l_tot = g(0, nj - 1).l;
        auto Pref = solver_case.gi.Pref;
        auto Tref = solver_case.gi.Tref;
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
                              gp.Pt = gp.Ps + (gp.Vm * gp.Vm + gp.Vu * gp.Vu) / 2. * gp.rho;
                              gp.H = gp.Tt * gp.Cp;
                              gp.s = std::log(pow(gp.Ts / Tref, gp.Cp) / std::pow(gp.Ps / Pref, gi.R));
                          });
        }
    }

    template <typename T>
    auto apply_mf(quiss::SolverCase<T> &solver_case)
    {
        auto &gi = solver_case.gi;
        size_t ni = gi.g.nRows();
        if(solver_case.inlet.mode == MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
        {
            // std::fill(gi.g.begin(0),gi.g.end(0),solver_case.inlet.Vm_moy);
            std::for_each(gi.g.begin(0),gi.g.end(0),
                [Vm = solver_case.inlet.Vm_moy](auto &gp)
                {
                    gp.Vm=Vm;
                }
            );
            solver_case.inlet.Mf = compute_massflow(gi.g, 0);
        }
        solver_case.mf.resize(ni);
        std::fill(solver_case.mf.begin(),solver_case.mf.end(),solver_case.inlet.Mf); // Todo add leakage and reintroduction
    }  

    template <typename T>
    auto curvature_solver(quiss::SolverCase<T> &solver_case)
    {
        auto &gi = solver_case.gi;
        size_t ni = gi.g.nRows();
        size_t nj = gi.g.nCols();
        if (ni < 3 && nj < 3)
        {
            throw std::length_error("Grid must have dimensions >= 3");
        }
        size_t max_geom=solver_case.max_geom;
        auto eps = solver_case.eps;
        auto tol_rel_mf =solver_case.tol_rel_mf;
        auto tol_pos = solver_case.tol_rel_pos * gi.g(0, nj - 1).l;

        apply_mf(solver_case);

        auto vmi = gi.g(0, 0).Vm;
        int count_geom = 0;

        solver_case.log.clear();
        auto converged = false;

        auto i_0 = 0;
        for (auto i = 0; i < ni; i++) // ensure value coherence
        {
            compute_gas_properties(gi, i);
        }

        apply_bc(solver_case);

        while (!converged)
        {

            // apply_bc(solver_case);

            for (auto i = i_0; i < ni; i++)
            {
                vmi = gi.g(i, 0).Vm;
                compute_vm_distribution(solver_case.mf[i], vmi, i, gi, tol_rel_mf, eps);
            }

            count_geom++;
            T delta_pos_max {};
            T delta_pos {};
            T delta_pos_moy {};
            solver_case.log.delta_pos.push_back(std::vector<T>{});
            for (auto i = i_0; i < ni; i++) // TODO run in //
            {
                compute_gas_properties(gi,i);
                compute_massflow_distribution(gi.g.begin(i), gi.g.end(i));
                if (i > 0 && count_geom < max_geom)
                {
                    delta_pos = balance_massflow(gi, i, tol_rel_mf * solver_case.mf[i]);
                    delta_pos /= gi.g(i,nj-1).l; // To get relative length
                    delta_pos_moy += delta_pos / (ni - 2.);
                    delta_pos_max = fmax(delta_pos_max,delta_pos);
                    solver_case.log.delta_pos.back().push_back(delta_pos);
                }
            }
            solver_case.log.delta_pos_max.push_back(delta_pos_max);
            solver_case.log.delta_pos_moy.push_back(delta_pos_moy);

            compute_grid_metrics(gi.g,gi.g_metrics,f_m,f_l);// TODO run in // 

            // apply_bc(solver_case);

            converged = (delta_pos_moy < tol_pos) || (count_geom >= max_geom);
            std::cout << count_geom << " " << delta_pos_max << " " << delta_pos_moy << std::endl;
        }
    }
}
