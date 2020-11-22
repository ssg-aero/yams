#include <gtest/gtest.h>
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>
#include <gridrender.h>
#include <gbs/bscinterp.h>
#include <gbs/extrema.h>

const double PI = acos(-1.);

const bool TESTS_USE_PLOT = false;
const double c_r = 287.04;

using namespace quiss;
using gbs::operator*;
using gbs::operator+;
using gbs::operator-;

auto f_sqVmq2 = [](const auto &gp) { return 0.5 * gp.Vm * gp.Vm; };

auto f_mf = [](const auto &gp) {
    return gp.Vm * gp.rho * cos(gp.phi + gp.gam) * 2 * PI * gp.y; //Schobeiri p. 273 adapted to Novak 1977 angle convention
};

auto f_rVu = [](const auto &gp) { return gp.y * gp.Vu; };
auto f_rTanBeta = [](const auto &gp) { return gp.y * tan(gp.bet); };
auto f_l = [](const auto &gp) { return gp.l; };
auto f_m = [](const auto &gp) { return gp.m; };
auto f_Mm = [](const auto &gp) { return gp.Vm / sqrt(gp.ga * c_r * gp.Ts); };
auto f_Vm = [](const auto &gp) { return gp.Vm; };

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
            m_part *=  D1_O2_i_bw(g, i, j, f_sqVmq2, f_m) / cos(beta);
        }
        else
        {
            m_part *=  D1_O1_i_bw(g, i, j, f_sqVmq2, f_m) / cos(beta);
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
            result *= D1_O2_i_bw(g, i, j, f_sqVmq2, f_m) / cos(gp.bet);
        }
        else
        {
            result *= D1_O1_i_bw(g, i, j, f_sqVmq2, f_m) / cos(gp.bet);
        }
    }
    return result;
};

auto eq_bet = [](const auto &g, size_t i, size_t j) {
    const auto gp = g(i, j);
    const auto Vm = gp.Vm;
    return D(g, i, j) * Vm * Vm + E(gp) * Vm + F(g, i, j);
};

auto eval_H_s = [](const auto &g, size_t i, size_t j) {
    if(i>0)
    {
        auto g1 = g(i - 1, j);
        auto g2 = g(i, j);
        g2.H = g1.H + g2.omg * (g2.y * g2.Vu - g1.y * g1.Vu); // So it works if g1 is not a blade
        g2.s = g1.s; // TODO modify entropy
    }
};

template <typename T, typename _Func>
inline void integrate_RK2_vm_sheet(T vmi, size_t i, MeridionalGrid<T> &g, _Func F)
{
    g(i, 0).Vm = vmi;
    size_t nj = g.nCols();

    eval_H_s(g,i,0);

    for (auto j = 1; j < nj; j++)
    {
        const auto gp = g(i, j);
        const auto gp_prev = g(i, j - 1);
        auto sqVmq2 = f_sqVmq2(gp_prev);
        auto dl = gp.l - gp_prev.l;
        auto sqVmq2_1 = sqVmq2 + F(g, i, j - 1) * dl;
        g(i, j).Vm = sqrt(2. * sqVmq2_1);// TODO add H and s to eq
        eval_H_s(g,i,j); // equations are using H and s (enthalpy and entropy)
        auto sqVmq2_2 = sqVmq2 + F(g, i, j) * dl; 
        g(i, j).Vm = 0.5 * (g(i, j).Vm + sqrt(2. * sqVmq2_2));
        eval_H_s(g,i,j);
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
            gp.rho = gp.Ps / (c_r * gp.Ts);
        });
}

template <typename Iterator>
inline void compute_massflow_distribution(Iterator begin, Iterator end)
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

template <typename T>
inline auto compute_massflow(const MeridionalGrid<T> &g, int i)
{
    auto nj = g.nCols();
    auto mf = 0.;
    for(auto j = 1 ; j < nj ; j++)
    {
        mf+=(f_mf(g(i,j)) + f_mf(g(i,j-1))) * (g(i,j).l - g(i,j-1).l) * 0.5;
    }
    return mf;
}

template <typename T, typename _Func>
inline auto eq_massflow(T vmi, MeridionalGrid<T> &g, int i, _Func F)
{
    auto nj = g.nCols();
    if (i > 0)
    {
        if (g(i, 0).iB == -1)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).Vu = g(i - 1, j).y * g(i - 1, j).Vu / g(i, j).y;
                g(i, j).bet = atan2(g(i, j).Vu, g(i, j).Vm);
            }
        }
        else
        {
            // TODO ecart flux-profile
        }
            
    }
    integrate_RK2_vm_sheet(vmi, i, g, F);
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
    return compute_massflow(g, i);
}

template <typename T>
auto newton_solve = [](const auto &crv,auto p, T u0, T tol_f = 1.e-3, T tol_u = 1.e-4, size_t it_max=100)
{
    auto delta = tol_u * 10., d0 = tol_f * 10.;
    auto u = u0;
    auto count =0;
    while (delta>tol_u && d0 > tol_f && count < it_max)
    {
        auto d0 = crv.value(u)-p;
        auto d1 = crv.value(u,1);
        auto d2 = crv.value(u,2);
        delta = d1*d0 / (d2*d0+d1*d1);
        u -= delta;
        count++;
    }
    return u;
};

template <typename T, typename _Func>
inline auto streamsheet_value_vector(const MeridionalGrid<T> &g, size_t i, size_t stride, _Func f)
{
    std::vector<T> vec(g.nCols());
    std::transform(
        std::execution::par,
        g.begin(i),
        g.end(i),
        g.begin(i + stride),
        vec.begin(),
        [&f](const auto &gp2,const auto &gp1) { return f(gp2)-f(gp1); });
        return vec;
}

template <typename T, typename _Func>
inline auto streamsheet_value_vector(const MeridionalGrid<T> &g, size_t i, _Func f)
{
    std::vector<T> vec(g.nCols());
    std::transform(
        std::execution::par,
        g.begin(i),
        g.end(i),
        vec.begin(),
        [&f](const auto &gp) { return f(gp); });
        return vec;
}

template <typename T, typename _Func>
inline auto find_streamsheet_max(const MeridionalGrid<T> &g, size_t i, _Func f)
{
    auto vec = streamsheet_value_vector(g,i,f);
    return *std::max_element(
        std::execution::par,
        vec.begin(),
        vec.end());
}

template <typename T, typename _Func>
inline auto find_streamsheet_max(const MeridionalGrid<T> &g, size_t i, size_t stride, _Func f)
{
    auto vec = streamsheet_value_vector(g,i,stride,f);
    return *std::max_element(
        std::execution::par,
        vec.begin(),
        vec.end());
}

template <typename T>
inline auto eval_RF(const MeridionalGrid<T> &g, int i, T B_)
{

    T dm_max = 0.;
    if (i > 0)
    {
        dm_max = find_streamsheet_max(g, i, -1, f_m);
    }
    if (i < g.nRows() - 1)
    {
        dm_max = fmax(dm_max, find_streamsheet_max(g, i + 1, -1, f_m));
    }
    auto Mm = fmax(0.95, find_streamsheet_max(g, i, f_Mm));
    auto l = (*(std::next(g.end(i), -1))).l;
    return 1. / (1. + (1 - Mm * Mm) * l * l / (B_ * dm_max * dm_max) );
}

template <typename T>
inline auto balance_massflow(MeridionalGrid<T> &g, int i, T tol_mf)
{
    auto nj = g.nCols();
    std::vector<gbs::constrType<T, 1, 1>> Q(nj);
    std::vector<gbs::constrType<T, 2, 1>> X(nj);
    std::vector<T> u(nj);
    auto l_tot = g(i, nj - 1).l;
    for (auto j = 0; j < nj; j++)
    {
        Q[j][0][0] = g(i, j).q;
        X[j][0][0] = g(i, j).x;
        X[j][0][1] = g(i, j).y;
        u[j] = g(i, j).l / l_tot;
    }
    auto f_Q = gbs::interpolate(Q, u, fmax(fmin(3, nj), 1), gbs::KnotsCalcMode::CHORD_LENGTH);
    auto f_X = gbs::interpolate(X, u, fmax(fmin(3, nj), 1), gbs::KnotsCalcMode::CHORD_LENGTH);
    auto delta_pos = 0.;
    // auto RF = 0.01;
    T B_ = fmin(8., 0.5 * 60. / g.nRows());
    auto RF = eval_RF(g,i,B_);
    
    for (auto j = 1; j < nj - 1; j++)
    {
        auto l = newton_solve<T>(f_Q, gbs::point<T,1>{g(0, j).q}, u[j]);
        auto X = f_X.value(l);
        auto dx = g(i, j).x - X[0];
        auto dy = g(i, j).y - X[1];
        delta_pos = fmax(fmax(fabs(dx), fabs(dy)), delta_pos);
        g(i, j).x += RF * (X[0] - g(i, j).x);
        g(i, j).y += RF * (X[1] - g(i, j).y);
    }
    return delta_pos;
}

template <typename T, typename _Func>
inline auto compute_vm_distribution(T mf,T &vmi,size_t i,MeridionalGrid<T> &g,_Func F,T tol_rel_mf,T eps)
{
    auto err_mf = tol_rel_mf * 10.;
    auto mf_ = 0., mf_pre = 0.; // mf shall allways be positive
    int count = 0;
    while (err_mf > tol_rel_mf && count < 100)
    {
        mf_pre = eq_massflow(vmi - eps, g, i, F);
        mf_ = eq_massflow(vmi, g, i, F);
        vmi = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
        err_mf = fabs(mf_ - mf) / mf;
        count++;
    }
}

template <typename T>
inline auto solve_grid(MeridionalGrid<T> &g,size_t max_geom = 100)
{
    compute_metrics(g);// TODO run in //
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    if(ni<3 && nj <3)
    {
        throw std::length_error("Grid must have dimensions >= 3");
    }
    auto vmi = g(0, 0).Vm;
    auto eps = 0.001;
    auto mf = compute_massflow(g, 0);
    auto tol_rel_mf = 0.01;
    auto tol_pos = 0.01 * g(0, nj - 1).l;
    int count_geom = 0;
    auto delta_pos_max = 0.;
    auto delta_pos = 0.;
    auto delta_pos_moy = 0.;
    auto converged = false;
    while (!converged)
    {
        for (auto i = 0; i < ni; i++)
        {
            if(g(i,0).iB==-1)
            {
                compute_vm_distribution(mf, vmi, i, g, eq_vu, tol_rel_mf, eps);
            }
            else
            {
                compute_vm_distribution(mf, vmi, i, g, eq_bet, tol_rel_mf, eps);
            }
            
        }
        delta_pos_max = 0.;
        delta_pos_moy = 0.;
        for (auto i = 0; i < ni; i++) // TODO run in //
        {
            compute_massflow_distribution(g.begin(i), g.end(i));
            if (i > 0)
            {
                delta_pos = balance_massflow(g, i, tol_rel_mf * mf);
                delta_pos_moy += delta_pos / (ni - 2.);
                delta_pos_max = fmax(delta_pos_max, delta_pos);
            }
        }
        compute_metrics(g);// TODO run in //
        count_geom++;
        converged = (delta_pos_moy < tol_pos) || (count_geom >= max_geom);
        // std::cerr << count_geom << " " << delta_pos_max << " " << delta_pos_moy << std::endl;
    }
    return converged;
}
TEST(tests_eq, mass_flow)
{
    size_t ni = 50;
    size_t nj = 11;
    MeridionalGrid<double> g(ni, nj);
    auto re = 1.;
    auto ri = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri, re, t1, t2, {0., 3.}, g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=10.;gp.rho=1.225; });
    compute_massflow_distribution(g.begin(0), g.end(0));
    auto S = 2 * PI * 3.;
    auto MF = S * 10. * 1.225;
    ASSERT_NEAR(MF, (*(g.end(0) - 1)).q, 1e-6);
}

TEST(tests_eq, functions)
{
    size_t ni = 50;
    size_t nj = 11;
    MeridionalGrid<double> g(ni, nj);
    auto re = 1.;
    auto ri = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri, re, t1, t2, {0., 3.}, g);
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
    auto re = 1.;
    auto ri = 2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri, re, t1, t2, {0., 3.}, g);
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
        compute_massflow_distribution(g.begin(i), g.end(i));
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
        compute_massflow_distribution(g.begin(i), g.end(i));
        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                ASSERT_LT(fabs(g(i,j).q-g(i-1,j).q),1e-5);
            }
        }
    }

    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, constant_flow_angle_clustered_grid)
{
    size_t ni = 10;
    size_t nj = 10;
    MeridionalGrid<double> g(ni, nj);
    auto r1 = 1.;
    auto r2 = 2.;

    make_uniform_clustered_grid(2,r2-r1,g,0.,0.,r1);
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
                ASSERT_LT((gp.Vm -vm_exact) / vm_exact * 100 , 1e-1); //less than 0.1 % on coarse grid
                // std::cerr << 100* (gp.Vm -vm_exact) / vm_exact << std::endl;
            });
        compute_massflow_distribution(g.begin(i), g.end(i));
        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                ASSERT_LT(fabs(g(i,j).q-g(i-1,j).q),1e-5);
            }
        }
    }
    
    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, constant_flow_vortex_circular)
{
    size_t ni = 50;
    size_t nj = 15;
    MeridionalGrid<double> g(ni,nj);
    auto re =  1.;
    auto ri =  2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri,re,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    compute_angles(g);
    compute_curvature(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=10.;gp.Vu=3; });

    auto vmi = 10.;
    for (auto i = 0; i < ni; i++)
    {

        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).Vu = g(i-1,j).y * g(i-1,j).Vu / g(i,j).y;
            }
        }
        integrate_RK2_vm_sheet(vmi, i, g, eq_vu);
        if (i > 0)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).bet = atan2(g(i,j).Vu,g(i,j).Vm);
            }
        }
        // std::cerr << compute_massflow(g,i) << std::endl;
    }

    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, constant_flow_vortex_circular_mass_flow_balance)
{
    size_t ni = 60;
    size_t nj = 15;
    MeridionalGrid<double> g(ni,nj);
    auto ri =  2.;
    auto re =  1.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri,re,t1,t2,{0.,3.},g);
    compute_metrics(g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=100.;gp.Vu=30;gp.H=gp.Cp*gp.Tt; });

    auto vmi = g(0,0).Vm;
    auto eps = 0.001;
    auto mf = compute_massflow(g, 0);
    auto tol_rel_mf = 0.01;
    int count_geom = 0;
    auto delta_pos_max = 0.;
    auto delta_pos=0.;
    auto delta_pos_moy =0.;
    do
    {
        delta_pos_max = 0.;
        delta_pos_moy = 0.;
        for (auto i = 0; i < ni; i++)
        {
            auto err_mf = tol_rel_mf * 10.;
            auto mf_ = 0., mf_pre = 0.; // mf shall allways be positive
            int count = 0;
            while (err_mf > tol_rel_mf && count < 100)
            {
                mf_pre = eq_massflow(vmi - eps, g, i, eq_vu);
                mf_    = eq_massflow(vmi, g, i, eq_vu);
                vmi    = vmi - eps * (mf_ - mf) / (mf_ - mf_pre);
                err_mf = fabs(mf_ - mf) / mf;
                count++;
            }
            ASSERT_LT(fabs(compute_massflow(g,i)-mf)/mf,tol_rel_mf);
        }
        for (auto i = 0; i < ni; i++) // TODO run in //
        {
            compute_massflow_distribution(g.begin(i), g.end(i));
            if (i > 0)
            {
                delta_pos = balance_massflow(g, i, tol_rel_mf * mf);
                delta_pos_moy += delta_pos / (ni-2.);
                delta_pos_max = fmax(delta_pos_max, delta_pos);
            }
        }
        compute_metrics(g);
        std::cerr << count_geom << " " << delta_pos_max << " " << delta_pos_moy << std::endl;
        count_geom++;
    }while (delta_pos_moy > 0.01 * g(0,nj-1).l && count_geom < 200);

    ASSERT_LT(count_geom , 200);

    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, solve_circular_grid)
{
    size_t ni = 60;
    size_t nj = 15;
    MeridionalGrid<double> g(ni,nj);
    auto re =  1.;
    auto ri =  2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(ri,re,t1,t2,{0.,3.},g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=100.;gp.Vu=30;gp.H=gp.Cp*gp.Tt; });

    ASSERT_TRUE(solve_grid(g));

    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, solve_igv)
{
    // size_t ni = 60;
    // size_t nj = 15;
    size_t ni = 21;
    size_t nj = 7;
    MeridionalGrid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto l  =  3.;
    auto b1 = PI/4;
    auto b2 = PI/3;
    make_straight_igv(r1,r2,l,b1,b2,g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=100.;gp.Vu=0.;gp.H=gp.Cp*gp.Tt; });

    ASSERT_TRUE(solve_grid(g));

    ASSERT_TRUE(g(ni-1,0).Vm>g(ni-1,nj-1).Vm);

    if (TESTS_USE_PLOT)
    {
        auto structuredGrid = make_vtkStructuredGrid(g);
        add_value(g, structuredGrid, "Vm", [](const auto &gp) { return gp.Vm; });
        structuredGrid->GetPointData()->SetActiveScalars("Vm");
        plot_vtkStructuredGrid(structuredGrid, true);
    }
}

TEST(tests_eq, streamsheet_value)
{
    size_t ni = 60;
    size_t nj = 15;
    MeridionalGrid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto l  =  3.;
    auto b1 = PI/4;
    auto b2 = PI/3;
    make_straight_igv(r1,r2,l,b1,b2,g);
    // init values
    std::for_each(g.begin(), g.end(), [](auto &gp) {gp.Vm=100.;gp.Vu=0.;gp.H=gp.Cp*gp.Tt; });
    auto Vm0 = streamsheet_value_vector(g,0,f_Vm);
    std::for_each(Vm0.begin(),Vm0.end(),[](const auto &v){ASSERT_LT(fabs(v-100.),1e-6);});
    auto dVm1 = streamsheet_value_vector(g,1,-1,f_Vm);
    std::for_each(dVm1.begin(),dVm1.end(),[](const auto &v){ASSERT_LT(fabs(v),1e-6);});
    for (auto j = 0; j < nj; j++)
    {
        g(0, j).Vm = 100. * j / (nj - 1.);
    }
    Vm0 = streamsheet_value_vector(g,0,f_Vm);
    for (auto j = 0; j < nj; j++)
    {
        ASSERT_LT(fabs(Vm0[j]-100.* j / (nj - 1.)),1e-6);
    }
    auto Vm_max = find_streamsheet_max(g,0,f_Vm);
    ASSERT_LT(fabs(Vm_max-100.),1e-6);
    auto dVm_max = find_streamsheet_max(g,1,-1,f_Vm);
    ASSERT_LT(fabs(dVm_max-100.),1e-6);
}