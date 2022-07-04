#include <gtest/gtest.h>
#include <datastorage.h>
#include <eqcurvaturesolver.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>


TEST(tests_solver_eq, direct_inv_eq)
{
    using namespace yams;
    using T = double;
    const T tol = 1e-6;

    size_t ni = 50;
    size_t nj = 15;
    MeridionalGrid<T> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto t1 = std::numbers::pi_v<T>;
    auto t2 = std::numbers::pi_v<T> * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);

    Array2d<Grid2dMetricsPoint<T>> g_metrics;
    compute_grid_metrics(g,g_metrics,fm,fl);
    
    T d_eth = T(1) / (ni - 1);
    T d_ksi = T(1) / (nj - 1);


    T Vm =100.;
    T Vu =100.;
    // T Vu =0.;
    T x  = 0.1;
    T y  = 0.2;
    // T U  = 200.;
    T U  = 0.;
    T omg= U / y;
    T Wu = Vu-U;
    T bet=atan2(Wu,Vm);
    T cur= 3.2;
    // T cur= 0.;

    std::for_each(
        g.begin(), g.end(),[=](auto &gp){
            gp.Vm = Vm;
            gp.Vu = Vu;
            gp.omg=omg;
            gp.cur=cur;
            gp.bet=bet;
        }
    );

    size_t i = ni /2, j = nj / 2;

    const auto &gp = g(i,j);

    auto d_bet = eq_bet(g,g_metrics,i,j,d_ksi, d_eth);
    auto D_    = D(g, g_metrics, i, j, d_ksi, d_eth);
    auto E_    = E(gp);
    auto F_    = F(g, g_metrics, i, j, d_ksi, d_eth);
    auto d_vu  = eq_vu(g,g_metrics,i,j,d_ksi, d_eth);
    auto G_    = G(gp);
    auto J_    = J(g, g_metrics, i, j, d_ksi, d_eth);
    auto K_    = K(g, g_metrics, i, j, d_ksi, d_eth);

    ASSERT_NEAR( 
        d_bet, 
        d_vu, 
        tol
    );

}