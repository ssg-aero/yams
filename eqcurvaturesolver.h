#pragma once
#include <diffop.h>
#include <numbers>
// Using Novak and Hearsey 1977
namespace quiss
{
    using std::cos;
    using std::sin;
    using std::sqrt;
    using std::tan;

    auto f_sqVmq2 = [](const auto &gp)
    { return 0.5 * gp.Vm * gp.Vm; };

    auto f_mf = [](const auto &gp)
    {
        return gp.Vm * gp.rho * cos(gp.phi + gp.gam) * 2 * std::numbers::pi * gp.y; //Schobeiri p. 273 adapted to Novak 1977 angle convention
    };

    auto f_I = [](const auto &gp)
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

    auto D = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
        auto tb = tan(gp.bet);
        auto ce = cos(gp.eps);
        auto se = sin(gp.eps);
        auto D1 = gp.cgp * ce * (gp.cur + tb / gp.y * gp.Dphi_Dth);
        auto D2 = gp.y > 0. ? -tb / gp.y * ce * D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_rTanBeta) : 0.;
        auto D3 = gp.y > 0. ? se / gp.y / cb * cb * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_rTanBeta) : 0.;
        return cb * cb * (D1 + D2 + D3);
    };

    auto E = [](const auto &gp)
    {
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
        auto ce = cos(gp.eps);
        auto se = sin(gp.eps);
        auto sp = sin(gp.phi);
        auto cg = cos(gp.gam);
        auto tb = tan(gp.bet);
        return 2. * gp.omg * cb * cb * (se * sp - cg * ce * tb);
    };

    auto F = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
        auto sb = sin(gp.bet);
        auto ce = cos(gp.eps);
        auto se = sin(gp.eps);
        auto F1 = ce * cb * cb * (D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_I) - gp.Ts * D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_S_));
        auto F2 = ce * cb * gp.sgp + se * sb;
        auto F3 = cb * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_sqVmq2) + cb * cb * cb * gp.Ts * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_S_);
        return F1 + F2 * F3;
    };

    auto eq_bet = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        const auto Vm = gp.Vm;
        return D(g, g_metrics, i, j, d_ksi, d_eth) * Vm * Vm + E(gp) * Vm + F(g, g_metrics, i, j, d_ksi, d_eth);
    };

    auto G = [](const auto &gp)
    {
        auto tb = tan(gp.bet);
        auto ce = cos(gp.eps);
        auto G1 = ce * gp.cgp * gp.cur;
        auto G2 = tb / gp.y * gp.Dphi_Dth;
        return G1 + G2;
    };

    auto J = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        auto se = sin(gp.eps);                                                    //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
        return se / gp.y * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_rVu); // simplification of cos beta with dS -> dm
    };

    auto K = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        auto beta = atan2(gp.Vu, gp.Vm);
        assert(abs(beta) < std::numbers::pi / 2.);
        auto cb = cos(beta);
        auto ce = cos(gp.eps);
        auto sb = sin(beta);
        auto se = sin(gp.eps);
        auto K1 = ce * (D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_I) - gp.Ts * D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_S_));
        auto K2 = -gp.Vu / gp.y * cb * D1_O2_so_dx2(g, g_metrics, i, j, d_ksi, d_eth, f_rVu);
        auto K3 = ce * gp.sgp * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_sqVmq2); // simplification of cos beta with dS -> dm
        auto K4 = (ce * cb * gp.sgp + se * sb) * gp.Ts * cb * D1_O2_so_dx1(g, g_metrics, i, j, d_ksi, d_eth, f_S_);
        return K1 + K2 + K3 + K4;
        /*
        auto tg_part = 0.;
        if (gp.y > 0.)
        {
            tg_part = -gp.Vu / gp.y;
            tg_part *= D1_O2_so_dx2(g,g_metrics,i,j,ksi,eth,f_rVu);

        }
        assert(tg_part == tg_part);
        auto m_part = 0.;
        if (i > 0)
        {
            m_part = sin(gp.gam + gp.phi) / cos(beta);
                        T v_ksi, v_eth;
            m_part *= D1_O2_so_dx2(g,g_metrics,i,j,ksi,eth,f_rVu) / cos(beta);

        }
        assert(m_part == m_part);
        return tg_part + m_part;
        */
    };

    auto eq_vu = [](const auto &g, const auto &g_metrics, size_t i, size_t j, auto d_ksi, auto d_eth)
    {
        const auto gp = g(i, j);
        const auto Vm = gp.Vm;
        return G(gp) * Vm * Vm + J(g, g_metrics, i, j, d_ksi, d_eth) * Vm + K(g, g_metrics, i, j, d_ksi, d_eth);
    };

}