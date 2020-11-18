#include <gtest/gtest.h>
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>
#include <gridrender.h>

const double PI = acos(-1.);

using namespace quiss;

auto f_sqVmq2 = [](const auto &gp) { return 0.5 * gp.Vm * gp.Vm; };

auto f_mf = [](const auto &gp) {
    return gp.Vm * gp.rho * cos(gp.phi - gp.gam) * 2 * PI * gp.y; //Schobeiri p. 273
};

auto f_rVu = [](const auto &gp) { return gp.y * gp.Vu; };
auto f_rTanBeta = [](const auto &gp) { return gp.y * tan(gp.bet); };
auto f_l = [](const auto &gp) { return gp.l; };
auto f_m = [](const auto &gp) { return gp.m; };

auto G = [](const auto &gp) {
    return cos(gp.phi + gp.gam) * gp.cur;
};

auto J = [](const auto &gp) {
    return 0.;
};

auto K = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    auto beta = atan2(gp.Vu, gp.Vm);
    auto tg_part = 0.;
    if (gp.y > 0.)
    {
        tg_part = -gp.Vu / gp.y;
        if(j>1)
        {
            tg_part *= D1_O2_j_bw(g, i, j, f_rVu, f_l);
        }
        else if(j==0) // warning, using value from previous iteration or imposed by solver
        {
            tg_part *= D1_O1_j_fw(g, i, j, f_rVu, f_l);
        }
        else // j == 1
        {
            tg_part *= D1_O1_j_bw(g, i, j, f_rVu, f_l);
        }
    }
    auto m_part = 0.;
    if(i>0)
    {
        m_part = sin(gp.gam + gp.phi) / cos(beta);
        if(i>1)
        {
            m_part *=  D1_O2_i_bw(g, i, j, f_sqVmq2, f_m);
        }
        else
        {
            m_part *=  D1_O1_i_bw(g, i, j, f_sqVmq2, f_m);
        }
    } 
    return tg_part + m_part;
};

auto eq_vu = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    const auto Vm = gp.Vm;
    return G(gp) * Vm * Vm + J(gp) * Vm + K(g, i, j);
};

auto D = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
    auto result = cos(gp.phi + gp.gam) * gp.cur;
    auto tg_part = 0.;
    if (gp.y > 0.)
    {
        tg_part = -tan(gp.bet) / gp.y;
        if (j > 1)
        {
            tg_part *= D1_O2_j_bw(g, i, j, f_rTanBeta, f_l);
        }
        else if (j == 0) // warning, using value from previous iteration or imposed by solver
        {
            tg_part *= D1_O1_j_fw(g, i, j, f_rTanBeta, f_l);
        }
        else // j == 1
        {
            tg_part *= D1_O1_j_bw(g, i, j, f_rTanBeta, f_l);
        }
    }
    result += tg_part;
    return cb * cb * result;
};

auto E = [](const auto &gp) {
    auto cb = cos(gp.bet); //TODO check if caching value cos(beta) tan(beta) cos(phi+gam)... improve speed
    return 2. * gp.omg * cb * cb * (-cos(gp.gam) * tan(gp.bet));
};

auto F = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    auto result = 0.;
    if (i > 0)
    {
        result = (cos(gp.bet) * sin(gp.gam + gp.phi));
        if (i > 1)
        {
            result *= D1_O2_i_bw(g, i, j, f_sqVmq2, f_m);
        }
        else
        {
            result *= D1_O1_i_bw(g, i, j, f_sqVmq2, f_m);
        }
    }
    return result;
};

auto eq_bet = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    const auto Vm = gp.Vm;
    return D(g, i, j) * Vm * Vm + E(gp) * Vm + F(g, i, j);
};

template <typename T, typename _Func>
inline void integrate_RK2_vm_sheet(T vmi, size_t i, MeridionalGrid<T> &g, _Func F)
{
    g(i, 0).Vm = vmi;
    size_t nj = g.nCols();

    for (auto j = 1; j < nj; j++)
    {
        const auto gp = g(i, j);
        const auto gp_prev = g(i, j - 1);
        auto sqVmq2 = f_sqVmq2(gp_prev);
        auto dl = gp.l - gp_prev.l;
        auto sqVmq2_1 = sqVmq2 + F(g, i, j - 1) * dl;
        g(i, j).Vm = sqrt(2. * sqVmq2_1);
        auto sqVmq2_2 = sqVmq2 + F(g, i, j) * dl;
        g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
    }
}

template <typename Iterator>
inline void compute_static_values(Iterator begin, Iterator end)
{
    std::for_each(
        std::execution::par,
        begin, end,
        [](auto &gp) {
            gp.Ts = gp.Tt - f_sqVmq2(gp) / gp.Cp;
            gp.Ps = gp.Pt * pow(gp.Ts / gp.Tt, gp.ga / (gp.ga - 1.));
            gp.rho = gp.Ps / (287.04 * gp.Ts);
        });
}

template <typename Iterator>
inline void compute_massflow(Iterator begin, Iterator end)
{
    std::transform(
        begin,
        std::next(end, -1),
        std::next(begin),
        std::next(begin),
        [](const auto &gp_prev, auto &gp) {
            gp.q = (f_mf(gp_prev) + f_mf(gp)) * (gp.l - gp_prev.l) * 0.5 + gp_prev.q;
            return gp;
        });
}

TEST(tests_eq, mass_flow)
{
    size_t ni = 50;
    size_t nj = 11;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1, r2, t1, t2, {0., 3.}, g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=10.;gp.rho=1.225; });
    compute_massflow(g.begin(0), g.end(0));
    auto S = 2 * PI * 3.;
    auto MF = S * 10. * 1.225;
    ASSERT_NEAR(MF, (*(g.end(0) - 1)).q, 1e-6);
}

TEST(tests_eq, functions)
{
    size_t ni = 50;
    size_t nj = 11;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1, r2, t1, t2, {0., 3.}, g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=10.;gp.rho=1.225; });
    ASSERT_NEAR(f_sqVmq2(g(0, 0)), 50., 1e-6);
}

TEST(tests_eq, eq_Vu)
{
    size_t ni = 50;
    size_t nj = 15;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1, r2, t1, t2, {0., 3.}, g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) { gp.Vm = 10.; });

    for (auto i = 0; i < ni; i++)
    {
        // eval tg part
        if (i == 0)
        {
            //boco
            std::for_each(g.begin(i), g.end(i), [](auto &gp) {
                gp.Vu = gp.Vm * tan(PI / 6.);
            });
        }
        else
        {
            std::transform(
                std::execution::par,
                g.begin(i - 1),
                g.end(i - 1),
                g.begin(i),
                g.begin(i),
                [](const auto &gp_prev, auto &gp) {
                    gp.Vu = gp_prev.Vu * gp_prev.y / gp.y;
                    return gp; // look for a way to avoid copying all
                });
        }

        integrate_RK2_vm_sheet(100., i, g,eq_vu);

        // p,t and rho a sufferer from a delay
        compute_static_values(g.begin(i), g.end(i));
    }
}

TEST(tests_eq, forced_vector_flow)
{
    size_t ni = 10;
    size_t nj = 500;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;

    make_uniform_grid(2., r2 - r1, g, 0, 0, r1);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    auto K_ = 3.;
    std::for_each(g.begin(), g.end(), [&K_](auto &gp) {gp.Vm=10.;gp.Vu=K_ * gp.y; });

    auto vmi = 10.;
    for (auto i = 0; i < ni; i++)
    {
        integrate_RK2_vm_sheet(vmi, i, g, eq_vu);
        std::for_each(
            g.begin(i),
            g.end(i),
            [&](const auto &gp) {
                auto vm_exact = sqrt(2 * K_ * K_ * (r1 - gp.y * gp.y) + vmi * vmi);
                ASSERT_LT((gp.Vm - vm_exact) / vm_exact * 100, 5e-5); //less than 0.00005 %
                //  std::cerr << 100* (gp.Vm -vm_exact) / vm_exact << std::endl;
            });
        compute_massflow(g.begin(i), g.end(i));
        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                ASSERT_LT(fabs(g(i,j).q-g(i-1,j).q),1e-5);
            }
        }
    }
}

TEST(tests_eq, constant_flow_angle)
{
    size_t ni = 10;
    size_t nj = 500;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;

    make_uniform_grid(2., r2 - r1, g, 0, 0, r1);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    auto a_ = PI / 6.;
    std::for_each(g.begin(), g.end(), [&a_](auto &gp) {gp.Vm=10.;gp.bet=PI / 6.; });

    auto vmi = 10.;
    for (auto i = 0; i < ni; i++)
    {
        integrate_RK2_vm_sheet(vmi, i, g,eq_bet);
        std::for_each(
            g.begin(i),
            g.end(i),
            [&](const auto &gp) {
                auto vm_exact = vmi * pow(r1/gp.y,sin(PI/6.)*sin(PI/6.));
                ASSERT_LT((gp.Vm -vm_exact) / vm_exact * 100 , 5e-5); //less than 0.00005 %
                // std::cerr << 100* (gp.Vm -vm_exact) / vm_exact << std::endl;
            });
        compute_massflow(g.begin(i), g.end(i));
        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                ASSERT_LT(fabs(g(i,j).q-g(i-1,j).q),1e-5);
            }
        }
    }

    
    // auto structuredGrid = make_vtkStructuredGrid(g);
    // add_value(g,structuredGrid,"Vm",[](const auto &gp){return gp.Vm;});
    // structuredGrid->GetPointData()->SetActiveScalars("Vm");
    // plot_vtkStructuredGrid(structuredGrid);
}