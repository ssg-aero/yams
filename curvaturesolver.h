#pragma once
#include <eqcurvaturesolver.h>
#include <meridionalsolvercase.h>
#include <gridmetrics.h>
#include <gbs/bscinterp.h>


namespace yams
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
    void integrate_RK2_vm_sheet(size_t i,  GridInfo<T> &gi, _Func F, int j_beg, int j_end, int j_stp)
    {
        auto &g = *gi.g;
        auto &g_metrics= *gi.g_metrics;
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

    template <typename T, typename _Func>
    void integrate_RK2_vm_sheet(T vmi, size_t i, GridInfo<T> &gi, _Func F, bool integrate)
    {
        auto &g = *gi.g;

        if(!integrate)
        {
            std::for_each(g.begin(i),g.end(i),[vmi](auto &gp){gp.Vm=vmi;});
            return;
        }

        size_t nj = g.nCols();
        int j_0 = nj * gi.j_0 + 1;
        g(i, j_0).Vm = vmi;
        integrate_RK2_vm_sheet(i,gi,F,j_0+1,nj, 1);
        integrate_RK2_vm_sheet(i,gi,F,j_0-1,-1,-1);
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
        auto &g = *gi.g;
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
    auto eq_massflow(T vmi, SolverCase<T> &solver_case, int i, bool integrate)
    {
        auto &gi= *solver_case.gi;
        auto &g = *gi.g;
        auto nj = g.nCols();
        
        if (g(i, 0).iB == -1)
        {
            for (auto j = 0; j < nj; j++)
            {
                if (i > 0)
                {
                    g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
                }
                g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg , g(i, j).Vm); // <- lag from previous
            }

            integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
        }
        else
        {
            if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DIRECT)
            {
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).bet = g(i, j).k;
                }
                integrate_RK2_vm_sheet(vmi, i, gi, eq_bet, integrate);
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet) + g(i, j).y * g(i, j).omg; // <- lag from previous
                }
            }
            else if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DESIGN_BETA_OUT)
            {
                auto i1 = solver_case.bld_info_lst[g(i, 0).iB].i1;
                auto i2 = solver_case.bld_info_lst[g(i, 0).iB].i2;
                if(i == i1)
                {
                    for (auto j = 0; j < nj; j++)
                    {
                        g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
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
                        auto bet_out = solver_case.bld_info_lst[g(i, j).iB].beta_out(l_rel);
                        auto bet_in  = g(i1, j).bet;
                        g(i, j).bet = bet_in *(1.-m_rel_loc) + bet_out * m_rel_loc;
                    }
                    integrate_RK2_vm_sheet(vmi, i, gi, eq_bet, integrate);
                    for (auto j = 0; j < nj; j++)
                    {
                        g(i, j).Vu = g(i, j).Vm * tan(g(i, j).bet) + g(i, j).y * g(i, j).omg; // <- lag from previous
                    }
                }
            }
            else if(solver_case.bld_info_lst[g(i, 0).iB].mode == MeridionalBladeMode::DESIGN_PHI)
            {
                auto i1 = solver_case.bld_info_lst[g(i, 0).iB].i1;
                auto i2 = solver_case.bld_info_lst[g(i, 0).iB].i2;
                auto omg= solver_case.bld_info_lst[g(i, 0).iB].omg;
                auto f_phi=solver_case.bld_info_lst[g(i, 0).iB].phi;
                if(i == i1)
                {
                    for (auto j = 0; j < nj; j++)
                    {
                        g(i, j).omg = omg; // TODO rem omg from grid
                        g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
                        g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg, g(i, j).Vm); // <- lag from previous
                    }
                    integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
                }
                else
                {
                    for (auto j = 0; j < nj; j++)
                    {
                        g(i, j).omg = omg; // TODO rem omg from grid
                        auto m_rel_loc = (g(i, j).m - g(i1, j).m) / (g(i2, j   ).m - g(i1, j).m);
                        auto l_rel     = (g(i, j).l - g(i , 0).l) / (g(i , nj-1).l - g(i , 0).l);
                        auto phi_out =f_phi(l_rel);
                        g(i, j).Vu = g(i1, j).Vu + m_rel_loc * phi_out * g(i, j).y * omg;
                    }
                    integrate_RK2_vm_sheet(vmi, i, gi, eq_vu, integrate);
                    for (auto j = 0; j < nj; j++)
                    {
                        g(i, j).bet = atan2(g(i, j).Vu - g(i, j).y * g(i, j).omg , g(i, j).Vm); // <- lag from previous
                    }
                }
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
        auto &g = *gi.g;
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
    // auto compute_vm_distribution(T mf, T vmi, size_t i,GridInfo<T> &gi, T tol_rel_mf, T eps, bool integrate)
    auto compute_vm_distribution(SolverCase<T> &solver_case, T vmi, size_t i, T tol_rel_mf, T eps, bool integrate)
    {
        auto mf     = solver_case.mf[i];
        auto err_mf = tol_rel_mf * 10.;
        auto mf_    = 0., mf_pre = 0.; // mf shall allways be strictly positive
        int count   = 0;
        auto max_count = solver_case.gi->vm_distribution_max_count;
        while (err_mf > tol_rel_mf && count < max_count)
        {
            mf_pre = eq_massflow(vmi - eps, solver_case, i, integrate);
            mf_ = eq_massflow(vmi, solver_case, i, integrate);
            vmi = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
            // assert(vmi >= 0.);
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
    }

    template <typename T>
    auto apply_mf(SolverCase<T> &solver_case)
    {
        auto &gi = *solver_case.gi;
        auto &g  = *gi.g;
        size_t ni = g.nRows();
        if(solver_case.inlet.mode == MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
        {
            // std::fill(gi.g.begin(0),gi.g.end(0),solver_case.inlet.Vm_moy);
            std::for_each(g.begin(0),g.end(0),
                [Vm = solver_case.inlet.Vm_moy](auto &gp)
                {
                    gp.Vm=Vm;
                }
            );
            solver_case.inlet.Mf = compute_massflow(g, 0);
            std::cout << "Mass flow set to: " << solver_case.inlet.Mf <<std::endl;
        }
        solver_case.mf.resize(ni);
        std::fill(solver_case.mf.begin(),solver_case.mf.end(),solver_case.inlet.Mf); // Todo add leakage and reintroduction
    }  

    template <typename T>
    // auto init_values(GridInfo<T> &gi,const std::vector<T> &mf, T tol_rel_mf, T eps)
     auto init_values(SolverCase<T> &solver_case, T tol_rel_mf, T eps)
    {
        auto &gi   = *solver_case.gi;
        auto &g    = *gi.g;
        size_t ni = g.nRows();
        auto vmi  = g(0, 0).Vm;
        for (auto i = 0; i < ni; i++)
        {
            compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps, false);
            compute_gas_properties(gi, i);
        }
        for (auto i = 0; i < ni; i++)
        {
            compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps, false);
            compute_gas_properties(gi, i);
        }
    }

    template <typename T>
    auto curvature_solver(SolverCase<T> &solver_case)
    {
        auto &gi = *solver_case.gi;
        auto &g = *gi.g;
        auto &g_metrics = *gi.g_metrics;
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        if (ni < 3 && nj < 3)
        {
            throw std::length_error("Grid must have dimensions >= 3");
        }
        size_t max_geom=solver_case.max_geom;
        auto eps = solver_case.eps;
        auto tol_rel_mf =solver_case.tol_rel_mf;
        auto tol_pos = solver_case.tol_rel_pos * g(0, nj - 1).l;

        apply_mf(solver_case);

        T vmi;
        int count_geom = 0;
        auto converged = false;
        auto i_0 = 0;

        solver_case.log.clear();

        apply_bc(solver_case);

        init_values(solver_case,tol_rel_mf, eps);

        while (!converged && (count_geom < max_geom))
        {

            // apply_bc(solver_case);

            for (auto i = i_0; i < ni; i++)
            {
                vmi = g(i, nj * gi.j_0 + 1).Vm;
                // compute_vm_distribution(solver_case.mf[i], vmi, i, gi, tol_rel_mf, eps,true);
                compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps,true);
                compute_gas_properties(gi,i);
            }

            count_geom++;
            T delta_pos_max {};
            T delta_pos {};
            T delta_pos_moy {};
            solver_case.log.delta_pos.push_back(std::vector<T>{});
            for (auto i = i_0; i < ni; i++) // TODO run in //
            {
                compute_massflow_distribution(g.begin(i), g.end(i));
                if (i > 0 && count_geom < max_geom)
                {
                    delta_pos = balance_massflow(gi, i, tol_rel_mf * solver_case.mf[i]);
                    delta_pos /= g(i,nj-1).l; // To get relative length
                    delta_pos_moy += delta_pos / (ni - 2.);
                    delta_pos_max = fmax(delta_pos_max,delta_pos);
                    solver_case.log.delta_pos.back().push_back(delta_pos);
                }
            }
            solver_case.log.delta_pos_max.push_back(delta_pos_max);
            solver_case.log.delta_pos_moy.push_back(delta_pos_moy);

            compute_grid_metrics(g,g_metrics,f_m,f_l);// TODO run in // 

            // apply_bc(solver_case);

            converged = delta_pos_moy < tol_pos;

        }
    }
}
