#pragma once
#include <eqcurvaturesolver.h>
#include <datastorage.h>
#include <gridmetrics.h>
#include <gbs/bscinterp.h>
#include <optional>

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

    template <typename T>
    struct SolverCase
    {
        GridInfo<T> &gi;
        std::vector<BladeInfo<T>> bld_info_lst;
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

    auto eval_H_s = [](auto &gi, size_t i, size_t j) {
        if(i>0)
        {
            auto &g = gi.g;
            const auto &g1 = g(i - 1, j);
            auto &g2 = g(i, j);
            g2.H = g1.H + g2.omg * (g2.y * g2.Vu - g1.y * g1.Vu); // So it also works if g1 is not a blade
            g2.I = g2.H - g2.omg * g2.y * g2.Vu;
            // g2.s = g1.s; // TODO modify entropy
            // g2.s = g(0,j).s + g2.Cp/g2.ga * std::log((g2.Ps/g1.Ps)/std::pow(g2.rho/g1.rho,g2.ga)) ;
            // g2.s = g1.s + g2.Cp/g2.ga * std::log((g2.Ps/g1.Ps)/std::pow(g2.rho/g1.rho,g2.ga)) ;
            g2.s = std::log(pow(g2.Ts/288.,g2.Cp)/std::pow(g2.Ps/1.01325e5,gi.R)) ;
        }
    };

    template <typename T, typename _Func>
    void integrate_RK2_vm_sheet(T vmi, size_t i, GridInfo<T> &gi, _Func F)
    {
        auto &g = gi.g;
        auto &g_metrics= gi.g_metrics;

        g(i, 0).Vm = vmi;
        size_t nj = g.nCols();

        eval_H_s(gi,i,0);

        for (auto j = 1; j < nj; j++)
        {
            const auto &gp = g(i, j);
            const auto &gp_prev = g(i, j - 1);

            auto sqVmq2 = f_sqVmq2(gp_prev);
            auto dl = gp.l - gp_prev.l;
            auto Fjm= F(g,g_metrics, i, j - 1,gi.d_ksi,gi.d_eth);
            assert(Fjm==Fjm);

            auto sqVmq2_1 = std::fmax(0.1,sqVmq2 + Fjm * dl);
            // g(i, j).Vm = std::fmin(sqrt(2. * sqVmq2_1),320.);// TODO add H and s to eq
            g(i, j).Vm = sqrt(2. * sqVmq2_1);// TODO add H and s to eq
            eval_H_s(gi,i,j); // equations are using H and s (enthalpy and entropy)
            auto Fj = F(g,g_metrics, i, j,gi.d_ksi,gi.d_eth);
            assert(Fj==Fj);

            auto sqVmq2_2 = std::fmax(0.1,sqVmq2 + Fj * dl); 
            // g(i, j).Vm = std::fmin(0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2)),320.);
            g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
            
            eval_H_s(gi,i,j);
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
                eval_H_s(gi, i, j);
                const auto &g1 = g(i - 1, j);
                auto &g2 = g(i, j);
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
            if(i!=0) eval_H_s(gi, i, j);
            // TODO update cp ga
        }
    }

    template <typename T, typename _Func>
    auto eq_massflow(T vmi, GridInfo<T> &gi, int i, _Func F)
    {
        auto &g = gi.g;
        auto nj = g.nCols();
        if (i > 0)
        {
            if (g(i, 0).iB == -1)
            {
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
                    g(i, j).bet = atan2(g(i, j).Vu, g(i, j).Vm);
                    // g(i, j).H = g(i-1, j).H;
                }
            }
            else
            {
                // TODO ecart flux-profile
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).bet = g(i, j).k;
                }
            }
                
        }
        integrate_RK2_vm_sheet(vmi, i, gi, F);
        if (i > 0)
        {
            if (g(i, 0).iB != -1)
            {
                for (auto j = 0; j < nj; j++)
                {
                    g(i, j).Vu = g(i, j).Vm * tan( g(i, j).bet );
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

    template <typename T, typename _Func>
    auto compute_vm_distribution(T mf, T &vmi, size_t i,GridInfo<T> &gi, _Func F, T tol_rel_mf, T eps)
    {
        auto err_mf = tol_rel_mf * 10.;
        auto mf_ = 0., mf_pre = 0.; // mf shall allways be strictly positive
        int count = 0;
        auto max_count = gi.vm_distribution_max_count;
        while (err_mf > tol_rel_mf && count < max_count)
        {
            mf_pre = eq_massflow(vmi - eps, gi, i, F);
            mf_ = eq_massflow(vmi, gi, i, F);
            vmi = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
            assert(vmi >= 0.);
            // vmi = fmin(fmax(0.1,vmi),360.);
            err_mf = fabs(mf_ - mf) / mf;
            count++;
        }
        if(count==max_count && err_mf <= tol_rel_mf)
        {
            std::cout << "Warning span: " << i << " did not converged." << std::endl;
        }
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

        auto mf = compute_massflow(gi.g, 0);
        auto vmi = gi.g(0, 0).Vm;
        int count_geom = 0;
        auto delta_pos_max = 0.;
        auto delta_pos = 0.;
        auto delta_pos_moy = 0.;
        auto converged = false;

        auto i_0 = 0;
        for (auto i = 0; i < ni; i++) // ensure value coherence
        {
            compute_gas_properties(gi, i);
        }

        while (!converged)
        {
            for (auto i = i_0; i < ni; i++)
            {
                if (gi.g(i, 0).iB == -1)
                {
                    compute_vm_distribution(mf, vmi, i, gi, eq_vu, tol_rel_mf, eps);
                }
                else
                {
                    compute_vm_distribution(mf, vmi, i, gi, eq_bet, tol_rel_mf, eps);
                }
            }
            count_geom++;
            delta_pos_max = 0.;
            delta_pos_moy = 0.;
            for (auto i = i_0; i < ni; i++) // TODO run in //
            {
                compute_gas_properties(gi,i);
                compute_massflow_distribution(gi.g.begin(i), gi.g.end(i));
                if (i > 0 && count_geom < max_geom)
                {
                    delta_pos = balance_massflow(gi, i, tol_rel_mf * mf);
                    delta_pos_moy += delta_pos / (ni - 2.);
                    delta_pos_max = fmax(delta_pos_max, delta_pos);
                }
            }

            compute_grid_metrics(gi.g,gi.g_metrics,f_m,f_l);// TODO run in // 

            converged = (delta_pos_moy < tol_pos) || (count_geom >= max_geom);
            std::cout << count_geom << " " << delta_pos_max << " " << delta_pos_moy << std::endl;
        }
    }
}