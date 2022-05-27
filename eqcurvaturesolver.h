#pragma once
#include <diffop.h>
#include <numbers>
#include <gbs/maths.h>
// Using Novak and Hearsey 1977
namespace yams
{
    using std::cos;
    using std::sin;
    using std::sqrt;
    using std::tan;

    template <typename T>
    T cap(T a, T a_max_abs)
    {
        return gbs::sgn(a) * std::min( std::abs(a), a_max_abs );
    }
    template <typename T>
    T cap_angle(T a)
    {
        return cap(a, 80 * (std::numbers::pi_v<T> / 180) );
    }

    const double c_r = 287.04;

    auto f_sqVmq2 = [](const auto &gp)
    { return 0.5 * gp.Vm * gp.Vm; };

    auto f_mf = [](const auto &gp)
    {
        // return gp.Vm * gp.rho * cos(gp.phi + gp.gam) * ( 2 * std::numbers::pi - gp.th_ ) * gp.y; //Schobeiri p. 273 adapted to Novak 1977 angle convention
        return gp.Vm * gp.rho * gp.cgp * ( 2 * std::numbers::pi - gp.th_ ) * gp.y; //Schobeiri p. 273 adapted to Novak 1977 angle convention
    };

    auto f_I = [](const auto &gp)
    { return gp.I; };
    auto f_H = [](const auto &gp)
    { return gp.H; };
    auto f_S_ = [](const auto &gp)
    { return gp.s; };
    auto f_rVu = [](const auto &gp)
    { return gp.y * gp.Vu; };
    auto f_rTanBeta = [](const auto &gp)
    { return gp.y * tan(gp.bet); };
    auto f_l = [](const auto &gp)
    { return gp.l; };
    auto f_m = [](const auto &gp)
    { return gp.m; };
    auto f_Mm = [](const auto &gp)
    { return gp.Vm / sqrt(gp.ga * c_r * gp.Ts); };
    auto f_Vm = [](const auto &gp)
    { return gp.Vm; };
    auto f_Wu = [](const auto &gp)
    { return gp.Vu - gp.omg * gp.y; };
    auto f_rWu = [](const auto &gp)
    { return f_Wu(gp) * gp.y; };


    auto D = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improves speed
        auto tb = tan(gp.bet);
        auto ce = cos(gp.eps);
        auto te = tan(gp.eps);
        // auto drtb_dl = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_rTanBeta);
        auto drtb_dl = g(i, j).drTb_dl;

        auto D1 = gp.cgp * gp.cur;
        auto D2 = gp.y > 0. ? -tb / gp.y * drtb_dl : 0.;
        auto D3 = gp.y > 0. ?  te / gp.y * gp.drtb_dm  : 0.;
        return cb * cb * (D1 + D2 + D3);
    };

    auto E = [](const auto &gp)
    {
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improves speed
        auto te = tan(gp.eps);
        auto sp = sin(gp.phi);
        auto cg = cos(gp.gam);
        auto tb = tan(gp.bet);
        return 2. * gp.omg * cb * cb * (te * sp - cg * tb);
    };

    auto F = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improves speed
        auto sb = sin(gp.bet);
        auto te = tan(gp.eps);
        auto ce = cos(gp.eps);
        // auto dI_dl = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_I);
        // auto dS_dl = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_S_);
        auto dI_dl = gp.dI_dl; 
        auto dS_dl = gp.dS_dl;

        auto F1 = cb * cb * (dI_dl - gp.Ts *dS_dl) / ce;
        auto F2 = cb * gp.sgp + te * sb;
        auto F3 = cb * ( gp.dsqVm_dm_2 +  cb * cb * gp.Ts * gp.ds_dm );
        return F1 + F2 * F3;
    };

    auto G = [](const auto &gp)
    {
        return gp.cgp * gp.cur;
    };

    auto J = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        auto te = tan(gp.eps);
        return gp.y > 0. ? te / gp.y * gp.drVu_dm : 0.;
    };

    auto K = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        auto beta = abs(gp.Vm) > 1e-1 ? atan2(gp.Vu, gp.Vm) : std::numbers::pi / 2.;
        // auto beta = abs(gp.Vm) > 1e-4 ? atan2(gp.Vu, gp.Vm) : 0.;
        // assert(abs(beta) < std::numbers::pi / 2.);
        auto cb = cos(beta);
        auto sb = sin(beta);
        auto te = tan(gp.eps);
        auto sgp= gp.sgp;
        
        // -> compute gradients in a separate subroutine for vectorization
        // auto dH_qdl = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_H); 
        // // auto dH_qdl = 0;
        // auto dS_qdl = D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_S_);
        // auto drVu_dl= D1_O2_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_rVu);
        auto dH_qdl = gp.dH_dl; 
        auto dS_qdl = gp.dS_dl;
        auto drVu_dl= gp.drVu_dl;

        auto K1 = ( dH_qdl - gp.Ts * dS_qdl );
        auto K2 = gp.y > 0. ?  -gp.Vu / gp.y * drVu_dl : 0.;
        auto K3 = sgp * gp.dsqVm_dm_2; 
        auto K4 = (cb * sgp + te * sb) * gp.Ts * cb * gp.ds_dm;
        return K1 + K2 + K3 + K4;
    };

    auto eq_bet = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        const auto Vm = gp.Vm;
        return D(g, g_metrics, i, j, d_ksi, d_eth) * Vm * Vm + E(gp) * Vm + F(g, g_metrics, i, j, d_ksi, d_eth);
    };

    auto eq_vu = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto &gp = g(i, j);
        auto Vm = gp.Vm;
        return G(gp) * Vm * Vm + J(g, g_metrics, i, j, d_ksi, d_eth) * Vm + K(g, g_metrics, i, j, d_ksi, d_eth);
    };

    // auto eq_vu = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    // {
    //     const auto &gp = g(i, j);
    //     auto cgp= gp.cgp;
    //     auto sgp= gp.sgp;
    //     auto cg = cos(gp.gam);
    //     auto Vm = gp.Vm;
    //     auto Vu = gp.Vu;
    //     auto Wu = Vu - gp.omg * gp.y;
    //     auto cu = gp.cur;

    //     auto K2 = gp.y > 0. ?  gp.Vu / gp.y * gp.drVu_dl : 0.;

    //     auto beta = abs(gp.Vm) > 1e-1 ? atan2(gp.Vu, gp.Vm) : std::numbers::pi / 2.;
    //     auto cb = cos(beta);

    //     return 
    //         sgp * gp.dsqVm_dm_2 
    //         + cgp * cu * Vm * Vm 
    //         - K2 
    //         + gp.dI_dl 
    //         - gp.Ts * gp.dS_dl 
    //         - 2 * gp.omg * Wu * cg 
    //         + sgp * cb * cb * gp.Ts * gp.ds_dm; 
    // };

}