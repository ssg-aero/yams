#pragma once

#include "curvaturesolver.h"
#include "meridionalsolvercase_experiemental.h"

namespace yams
{

    template <typename T>
    auto balance_massflow(SolverCaseBiPass<T> &bip_case, int i, T tol_mf)
    {
        T span_geom_residual{};
        if( i==bip_case.iBip())
        {
            return span_geom_residual;
        }

        if( i >= bip_case.iBip())
        {
            span_geom_residual = std::max(
                balance_massflow(bip_case.primary(), i, tol_mf),
                balance_massflow(bip_case.secondary(), i, tol_mf)
            );
        }
        else
        {
            auto nj_prim = bip_case.primary().gi->nj;
            auto nj_sec = bip_case.secondary().gi->nj;
            auto nj = nj_prim + nj_sec - 1;
            auto & g_prim = *(bip_case.primary().gi->g);
            auto & g_sec  = *(bip_case.secondary().gi->g);
            auto l_prim = std::next(g_prim.end(i),-1)->l;
            auto l_sec = std::next(g_sec.end(i),-1)->l;
            auto l_tot = l_prim + l_sec;
            auto mf_prim_ref = std::next(g_prim.end(bip_case.iBip()),-1)->q;
            auto mf_sec_ref  =std::next(g_sec.end(bip_case.iBip()),-1)->q;
            auto mf_prim = std::next(g_prim.end(i),-1)->q;
            auto mf_sec  = std::next(g_sec.end(i),-1)->q;

            // Compute full mass_flow
            std::vector<T> q(nj);
            for (auto j = 0; j < nj; j++)
            {
                // Use last cumulative mass flow rather than specified to perfectly match and then solve better
                if(j<nj_prim)
                {
                    q[j] = g_prim(i, j).q * mf_prim_ref / mf_prim;
                }
                else
                {
                    q[j] = g_sec(i, j - nj_prim + 1).q * mf_sec_ref / mf_sec + mf_prim_ref;
                }
            }
            // Compute span pos law
            gbs::points_vector<T, 2> X(nj);
            std::vector<T> u(nj);
            for (auto j = 0; j < nj; j++)
            {
                if(j<nj_prim)
                {
                    const auto &gp = g_prim(i, j);
                    X[j][0] = gp.x;
                    X[j][1] = gp.y;
                    u[j]    = gp.l / l_tot;
                }
                else
                {
                    const auto &gp = g_sec(i, j - nj_prim + 1);
                    X[j][0] = gp.x;
                    X[j][1] = gp.y;
                    u[j]    = (gp.l + l_prim) / l_tot;
                }
            }
            size_t p = fmax(fmin(3, nj), 1);
            auto f_X =  gbs::interpolate(X, u, p);
            auto f_Q = gbs::interpolate(q, u, f_X.degree());
            //

            // auto span_geom_residual = 0.;
            // TODO use specifc params
            auto RF = bip_case.primary().gi->RF;
            auto tol_f = bip_case.primary().gi->tol_newtow_mf_f;
            auto tol_u = bip_case.primary().gi->tol_newtow_mf_u;

            auto nj_ = nj-1;
            for (auto j = 1; j < nj_; j++)
            {
                auto [u1, u2] = f_Q.bounds();
                auto q  = j<nj_prim ? g_prim(bip_case.iBip(), j).q : g_sec(bip_case.iBip(), j - nj_prim +1).q + mf_prim;
                auto [l, delta, count] = newton_solve<T>(f_Q, q, u[j], u1, u2, tol_f, tol_u);
                T dx{},dy{};
                if(j<nj_prim)
                {
                    auto X = f_X.value(l);
                    dx = g_prim(i, j).x - X[0];
                    dy = g_prim(i, j).y - X[1];
                    g_prim(i, j).x += RF * (X[0] - g_prim(i, j).x);
                    g_prim(i, j).y += RF * (X[1] - g_prim(i, j).y);
                }
                else
                {
                    auto X = f_X.value(l);
                    auto j_ = j - nj_prim + 1;
                    dx = g_sec(i, j_).x - X[0];
                    dy = g_sec(i, j_).y - X[1];
                    g_sec(i, j_).x += RF * (X[0] - g_sec(i, j_).x);
                    g_sec(i, j_).y += RF * (X[1] - g_sec(i, j_).y);
                }
                span_geom_residual = fmax(fmax(fabs(dx), fabs(dy)), span_geom_residual);               
            }
            g_sec(i, 0).x = g_prim(i,nj_prim-1).x;
            g_sec(i, 0).y = g_prim(i,nj_prim-1).y;
        }

        return span_geom_residual;

    }


    template <typename T>
    auto apply_mf(SolverCaseBiPass<T> &bip)
    {
        bip.primary().inlet.mode = bip.inlet().inlet.mode;
        bip.secondary().inlet.mode = bip.inlet().inlet.mode;

        if(bip.inlet().inlet.mode == MeridionalBC::INLET_VmMoy_Ts_Ps_Vu)
        {
            std::for_each(bip.primary().gi->g->begin(0),bip.primary().gi->g->end(0),
                [Vm = bip.inlet().inlet.Vm_moy](auto &gp)
                {
                    gp.Vm=Vm;
                }
            );
            std::for_each(bip.secondary().gi->g->begin(0),bip.secondary().gi->g->end(0),
                [Vm = bip.inlet().inlet.Vm_moy](auto &gp)
                {
                    gp.Vm=Vm;
                }
            );
            bip.inlet().inlet.Mf = compute_massflow(*(bip.primary().gi->g), 0) + compute_massflow(*(bip.secondary().gi->g), 0);
            bip.primary().inlet.Mf = bip.inlet().inlet.Mf / ( 1 + bip.BPR);
            bip.secondary().inlet.Mf =bip.primary().inlet.Mf * bip.BPR;
            // std::cout << "Mass flow set to: " << solver_case.inlet.Mf <<std::endl;
        }
        else
        {
            throw std::invalid_argument("Unsupported yet");
        }
        bip.primary().mf.resize(bip.primary().gi->ni);
        bip.secondary().mf.resize(bip.secondary().gi->ni);
        std::fill(bip.primary().mf.begin(),bip.primary().mf.end(),bip.primary().inlet.Mf); // Todo add leakage and reintroduction
        std::fill(bip.secondary().mf.begin(),bip.secondary().mf.end(),bip.secondary().inlet.Mf); // Todo add leakage and reintroduction
    }  
    template <typename T, auto ExPo = std::execution::par>
    auto curvature_solver(SolverCaseBiPass<T> &bip_case)
    {
        auto &prim = bip_case.primary();
        auto &sec  = bip_case.secondary();

        size_t max_geom=prim.max_geom;
        auto eps = prim.eps;
        auto tol_rel_mf =prim.tol_rel_mf;
        // auto tol_pos = prim.tol_rel_pos * ( prim.gi->g(0, prim.gi->nj - 1).l + sec.gi->g(0, sec.gi->nj - 1).l );
        auto l_prim = prim.gi->g->end(0)->l;
        auto l_sec  = sec.gi->g->end(0)->l;
        auto tol_pos = prim.tol_rel_pos * ( l_prim + l_sec );
        T vmi;
        int count_geom = 0;
        auto converged = false;
        // auto i_0 = 0;
        bip_case.primary().log.clear();
        bip_case.secondary().log.clear();
        T delta_pos_max {};
        T delta_pos {};
        T delta_pos_moy {};

        // compute spans mass flow
        apply_mf(bip_case);
        // apply boundary conditions
        apply_bc(prim);
        apply_bc(sec);
        // innit values
        init_values(prim,tol_rel_mf, eps);
        init_values(sec,tol_rel_mf, eps);
        // apply rotation sppeds
        apply_rotation_speeds(prim);
        apply_rotation_speeds(sec);

        std::vector<T> delta_pos_array_prim(prim.gi->ni);
        std::vector<T> delta_pos_array_sec(sec.gi->ni);
        auto span_range_prim = gbs::make_range<size_t>(0,prim.gi->ni-1);
        auto span_range_sec  = gbs::make_range<size_t>(0,sec.gi->ni-1);
        while (!converged && (count_geom < 1000))
        {
            // integrate radial eq equation and update gas properties
            compute_vm_distribution(prim, tol_rel_mf, eps, true, 0 );
            compute_vm_distribution(sec, tol_rel_mf, eps, true, 0 );
            // Compute mass flow distribution
           std::for_each(
                ExPo,
                span_range_prim.begin(), span_range_prim.end(),
                [&](const auto &i){
                    compute_massflow_distribution(prim.gi->g->begin(i), prim.gi->g->end(i));
                }
           );
           std::for_each(
                ExPo,
                span_range_sec.begin(), span_range_sec.end(),
                [&](const auto &i){
                    compute_massflow_distribution(sec.gi->g->begin(i), sec.gi->g->end(i));
                }
           );
            // relocate streams to balance mass flow
           std::transform(
                // ExPo,
                span_range_prim.begin(),
                std::next(span_range_prim.begin(), bip_case.iBip()),
                delta_pos_array_prim.begin(),
                [&](const auto &i)
                {
                    return balance_massflow(bip_case, i, tol_rel_mf * bip_case.inlet().inlet.Mf);// / g(i,nj-1).l;
                }
           );
            std::transform(
                // ExPo,
                std::next(span_range_prim.begin(), bip_case.iBip()+1),
                span_range_prim.end(),
                std::next(delta_pos_array_prim.begin(), bip_case.iBip()+1),
                [&](const auto &i)
                {
                    return balance_massflow(prim, i, tol_rel_mf * bip_case.inlet().inlet.Mf);// / g(i,nj-1).l;
                }
           );
           std::transform(
                // ExPo,
                std::next(span_range_sec.begin(), bip_case.iBip()+1),
                span_range_sec.end(),
                std::next(delta_pos_array_sec.begin(), bip_case.iBip()+1),
                [&](const auto &i)
                {
                    return balance_massflow(sec, i, tol_rel_mf * bip_case.inlet().inlet.Mf);// / g(i,nj-1).l;
                }
           );


            count_geom++;
        }
    }




    template <typename T, auto ExPo = std::execution::par>
    auto compute_vm_distribution2(SolverCase<T> &solver_case, T tol_rel_mf, T eps, bool integrate)
    {
        auto &gi = *solver_case.gi;
        auto &g = *gi.g;
        size_t ni = g.nRows();
        size_t nj = g.nCols();

        gbs::VectorX<T> vmi(ni), vmi_1(ni), vmi_2(ni), F(ni);
        auto indexes = gbs::make_range<size_t>(0, ni - 1);

        std::transform(
            ExPo,
            indexes.begin(), indexes.end(),
            vmi.begin(),
            [nj, &g, &gi] (size_t i)
            {
                return g(i, std::round((nj - 1 ) * gi.j_0)).Vm;
            }
        );

        gbs::MatrixX<T> J(ni, ni);
        T err{};
        size_t count{};
        do{
            for(size_t i{}; i < ni; i++)
            {
                std::copy(
                    ExPo,
                    vmi.begin(), vmi.end(),vmi_1.begin()
                );
                std::copy(
                    ExPo,
                    vmi.begin(), vmi.end(),vmi_2.begin()
                );
                vmi_1(i) -= eps * 0.5;
                vmi_2(i) += eps * 0.5;
                // for(size_t j{}; j < ni; j++)
                std::for_each(
                    ExPo,
                    indexes.begin(), indexes.end(),
                    [&](size_t j)
                        {
                            auto F2 = eq_massflow(vmi_2(j), solver_case, j, true) - solver_case.mf[j];
                            auto F1 = eq_massflow(vmi_1(j), solver_case, j, true) - solver_case.mf[j];
                            J(i,j) = ( F2 - F1 ) / eps;
                        }
                );
            }

            // std::cout << J << std::endl;
            auto J_inv = J.partialPivLu();
            // auto J_inv = J.llt();
            

            for(size_t i{}; i < ni; i++) // run in par
            {
                F(i) = eq_massflow(vmi(i), solver_case, i, true) - solver_case.mf[i];
            }

            auto delta = J_inv.solve(F);
            err = std::sqrt(delta.squaredNorm());
            vmi -= delta;
            // std::cout << delta << std::endl << std::endl << err << std::endl;

            for(size_t i{}; i < ni; i++)
            {
                eq_massflow(vmi(i), solver_case, i, true);
                int j_0 = std::round((nj - 1 ) * gi.j_0);
                g(i, j_0).Vm = vmi(i);
                compute_gas_properties(solver_case,i);
            }
            count++;
        } while (err>1e-6  && count < 500);
        std::cout << "Count: " << count << " err: " << err << std::endl;
    }

   template <typename T>
    auto init_values(SolverCaseSet<T> &solver_case, T tol_rel_mf, T eps)
    {
        auto &gi   = *solver_case.gi;
        auto &g    = *gi.g;
        size_t ni = g.nRows();
        auto vmi  = g(0, 0).Vm;
        for (auto i = 0; i < ni; i++)
        {
            compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps, false);
            compute_gas_properties(solver_case, i);
        }
        for (auto i = 0; i < ni; i++)
        {
            compute_vm_distribution(solver_case, vmi, i, tol_rel_mf, eps, false);
            compute_gas_properties(solver_case, i);
        }
    }

    template <typename T, auto ExPo = std::execution::par>
    auto curvature_solver(SolverCaseSet<T> &set)
    {
        auto solver_case_in = set.cases().front();


        // apply boundary conditions
        apply_bc(*solver_case_in);

        solver_case_in->mf_ref_span = solver_case_in->gi->ni - 1;

        std::for_each(
            std::next(set.cases().begin()), set.cases().end(),
            [](auto &solver_case)
            {
                solver_case->inlet.mode = MeridionalBC::CON;
                solver_case->mf_ref_span = 0; // The starting span will be updated by connection
            }
        );

        for(int i = 0 ; i < 1000 ; i++)
        {
            std::for_each( // Todo solve according connection tree
                set.cases().begin(), set.cases().end(),
                [&set](auto &solver_case)
                {
                    auto &gi = *solver_case->gi;
                    auto &g = *gi.g;
                    auto &g_metrics = *gi.g_metrics;
                    size_t ni = g.nRows();
                    size_t nj = g.nCols();
                    // size_t max_geom=solver_case.max_geom;
                    auto eps = solver_case->eps;
                    auto tol_rel_mf =solver_case->tol_rel_mf;
                    // auto tol_pos = solver_case.tol_rel_pos * g(0, nj - 1).l;

                    // compute spans mass flow
                    apply_mf(*solver_case);
                    // // get connections curvature influence
                    // std::for_each(
                    //     set.connections().begin(),set.connections().end(),
                    //     [&solver_case](auto & connection)
                    //     {
                    //         if(connection.isLeft(solver_case))
                    //         {
                    //             connection.interpolateCurvature();
                    //         }
                    //     }
                    // );
                    // integrate radial eq equation and update gas properties
                    for (auto i = 0; i < ni; i++)
                    {
                        if( ( solver_case->inlet.mode != MeridionalBC::INLET_Vm_Ts_Ps_Vu && solver_case->inlet.mode != MeridionalBC::CON )
                             || i != 0)
                        {
                            auto vmi = g(i, std::round((nj - 1 ) * gi.j_0)).Vm;
                            compute_vm_distribution(*solver_case, vmi, i, tol_rel_mf, eps,true);
                        }
                        // compute_gas_properties<T>(solver_case,i);
                    }
                    if( !solver_case->relocate )
                    {
                        return;
                    }
                    auto span_range = gbs::make_range<size_t>(0,ni-1);
                    // Compute mass flow distribution
                    std::for_each(
                            ExPo,
                            span_range.begin(), span_range.end(),
                            [&](const auto &i){
                                compute_massflow_distribution(g.begin(i), g.end(i));
                            }
                    );
                    // relocate streams to balance mass flow
                    // std::transform(
                    std::for_each(
                            ExPo,
                            span_range.begin(),
                            span_range.end(),
                            // delta_pos_array.begin(),
                            [&](const auto &i)
                            {
                                return balance_massflow(*solver_case, i, tol_rel_mf * solver_case->mf[i]) / g(i,nj-1).l;
                            }
                    );

                    // copy to rights sides, this includes curvature
                    std::for_each(
                        set.connections().begin(),set.connections().end(),
                        [&solver_case](auto & connection)
                        {
                            if(connection.isLeft(solver_case))
                            {
                                connection.copyLeftToRight();
                            }
                        }
                    );
                    compute_grid_metrics(g,g_metrics,f_m,f_l);// TODO run in // 
                //  // set connections curvature influence
                //     std::for_each(
                //         set.connections().begin(),set.connections().end(),
                //         [&solver_case](auto & connection)
                //         {
                //             if(connection.isLeft(solver_case))
                //             {
                //                 connection.interpolateCurvature();
                //             }
                //         }
                //     );
                }
            );

            std::for_each(
                set.connections().begin(),set.connections().end(),
                [](auto & connection)
                {
                    connection.interpolateCurvature();
                }
            );
        }
    }

}