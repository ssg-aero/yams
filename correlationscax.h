#pragma once
#include <gbs/maths.h>
auto Kti(auto tb, auto c)
{    
    auto q = 0.28/(0.1+(tb/c)**0.3);
    return (10*tb/c)**q;
}

auto Ktd( auto tb, auto c )
{
    return 6.25*(tb/c)+37.5*std::pow(tb/c,2);
}

auto design_angle_of_attack_herrig(auto sig,auto th,auto a,auto c,auto Kti,auto Ksh=1.0)
{
    auto e = 0.65-0.002*th;
    return (3.6*Ksh*Kti+0.3532*th*std::pow((a/c),0.25))*sig**e;
}

auto i0_opt_10(auto b1, auto sig)
{
    auto p = 0.914+std::pow(sig,3)/160;
    return (b1**p) / (5+46*std::exp(-2.3*sig)) - (0.1 * std::pow(sig,3)) * std::exp((b1-70)/4);
}

auto i_opt_lieblein(auto th, auto b1, auto sig, auto Kti, auto Ksh=1.0)
{
    auto n = 0.025 * sig - 0.06 - ((b1/90)**(1+1.2*sig))/(1.5+0.43*sig);
    auto i = i0_opt_10(b1,sig);
    return Ksh * Kti * i + n * th;
}

auto m10_naca65 = [](auto b1)
{
    auto x = b1 / 100;
    return 0.17 * 0.0333 * x + 0.333 * x*x;
}


auto m(const auto m10, auto b1, auto sig)
{
    auto x = b1 / 100;
    auto b = 0.9625 - 0.17 * x - 0.85 * x*3;
    return m10(b1) * pow(sig, b);
}


auto d_opt_10(auto sig, auto b1)
{
    return 0.01 * sig * b1 + (0.74 * std::pow(sig, 1.9) + 3 * sig)*std::pow((b1/90),(1.67+1.09*sig));
}


auto d_opt(auto sig, auto b1, auto th, auto Ktd, auto Ksh=1.0, auto m10=m10_naca65)
{
    return Ksh * Ktd * d_opt_10(sig, b1) + m(m10, b1, sig) * th;
}


auto d_d_opt_q_d_i(auto sig, auto b1)
{
    return (1+(sig+0.25*std::pow(sig,4)) * std::pow((b1/53),2.5))/std::exp(3.1*sig);
}


auto  d(auto d_opt, auto sig, auto b1, auto i, auto i_opt, auto Vm1, auto Vm2)
{
    auto slope_ = d_d_opt_q_d_i(sig, b1);
    return d_opt+slope_*(i-i_opt)+10*(1-Vm2/Vm1);
}

auto dev_cmp_ax(auto b1,auto k1,auto k2,auto sig12,auto th12,auto g12,auto a12,auto c12,auto tb12,auto Vm1,auto Vm2,auto Ksh=1.0):
{
    auto a1_opt_herrig = design_angle_of_attack_herrig(sig12,th12,a12,c12,Kti(tb12,c12));
    auto b1_opt = a1_opt_herrig + g12;
    auto i1_opt_herrig = b1_opt - k1;

    auto eq_a1_opt_lieblein = [](auto x)
    {
        return i_opt_lieblein(th12,x,sig12,Kti(tb12,c12),Ksh) - (x-k1);
    }

    auto i1_opt_lieblein = gbs::newtown_solve(eq_a1_opt_lieblein,0,b1)[0] - k1;
    auto b1_opt = i1_opt_lieblein + k1;

    auto d2_opt = d_opt(sig12,b1_opt,th12,Ktd(tb12,c12));
    auto d2     = d(d2_opt,sig12,b1,b1-k1,i1_opt_lieblein,Vm1,Vm2);
    auto b2     = k2+d2;
    auto b2_opt = k2+d2_opt;
}