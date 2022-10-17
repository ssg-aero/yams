#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-mesh/tfi.h>
#include <gbs-mesh/smoothing.h>
#include <gbs-render/vtkgridrender.h>

#include <gas/gasModel.h>
#include <euler2d/euler2dSolver.h>

#include <chrono>

using namespace gbs;
using T = double;
using BSCurve2d = BSCurve<T,2>;
template <typename T>
using vector = std::vector<T>;
template <typename... _Types>
using tuple = std::tuple<_Types...>;



namespace yams{
    template <typename... _Types>
    inline auto make_tuple(_Types&&... _Args) -> tuple<_Types...>
    {
        return std::make_tuple(_Args...);
    }
    template <size_t i>
    auto get(auto &tuple)
    {
        return std::get<i>(tuple);
    }
    template <typename Tuple, typename Functor, size_t Index = 0>
    auto tuple_for_each(const Tuple &tpl, const Functor &f) -> void
    {
        constexpr auto tuple_size = std::tuple_size_v<Tuple>;
        if constexpr (Index < tuple_size)
        {
            f(std::get<Index>(tpl));
            tuple_for_each<Tuple, Functor, Index + 1>(tpl, f);
        }
    }
}

inline auto build_channel_mesh(size_t ni_ = 30,size_t nj_ = 30, T bump = 0.1, T l = 2., T h = 1.)
{
    auto south_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0,0.0},
            {0.2 * l,0.0},
            {0.4 * l,bump},
            {0.6 * l,bump},
            {0.8 * l,0.0},
            {1.0 * l,0.0}
        },
        {0.,0.25, 0.75, 1.},
        {4,1,1, 4},
        3
    });
    auto north_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0, h},
            {l, h}
        },
        {0.,1.},
        {2, 2},
        1
    });
    auto west_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {0.0,0.0},
            {0.0,h}
        },
        {0.,1.},
        {2, 2},
        1
    });
    auto east_crv = std::make_shared<BSCurve2d>(BSCurve2d{
        {
            {l,0.0},
            {l,h}
        },
        {0.,1.},
        {2, 2},
        1
    });

    return tfi_mesh_2d<T,2,1,1>({south_crv, north_crv},{west_crv, east_crv},{0., 1.},{0., 1.}, ni_, nj_);
}


struct FlowData
{
    tuple<vector<T>, vector<T>, vector<T>, vector<T>> d; // tuple can be used thanks to thrust in cuda
    FlowData(size_t n) {
        d = yams::make_tuple(
            vector<T>(n),
            vector<T>(n),
            vector<T>(n),
            vector<T>(n)
        );
    }
    void init(T v1, T v2, T v3, T v4)
    {
        std::fill(std::get<0>(d).begin(), std::get<0>(d).end(), v1);
        std::fill(std::get<1>(d).begin(), std::get<1>(d).end(), v2);
        std::fill(std::get<2>(d).begin(), std::get<2>(d).end(), v3);
        std::fill(std::get<3>(d).begin(), std::get<3>(d).end(), v4);
    }
};

// template <typename T>
// class rk_integrator
// {
// private:
//     T
//     const std::array<T,5> rk_coefs;
// public:
//     rk_integrator(const std::array<T,5> &rk_coefs = {0.0533,0.1263,0.2375, 0.4414,1.0}) rk_coefs{rk_coefs} {}
//     template <typename U, typename _Func>
//     U integrate(U W, T dt, _Func R)
//     {
//         for(auto alf : rk_coefs)
//         {

//         }
//     }
// };

// TEST(tests_euler2D, rk_integrator)
// {

// }



TEST(tests_euler2D, gris_access_speeds)
{

    size_t ni_{51}, nj_{25};
    // size_t ni_{5}, nj_{3};
    T bump = 0.3;
    auto [ pts, ni, nj, n_iso_ksi, n_iso_eth] = build_channel_mesh(ni_, nj_, bump);
    elliptic_structured_smoothing(pts, ni, 0, nj-1, 0, ni-1, 300);
    auto njm=nj-1;
    auto nim=ni-1;
    auto njp=nj+1;
    auto nip=ni+1;
    auto nm = nim * njm; // N cells without ghosts
    auto np = nip * njp; // N cells with ghosts
    auto n = ni * nj;    // N vertices

    // Define flow description
    FlowData U{np}; // U = (rho, rho*u, rho*v, rho*et)
    FlowData R{np};
    T Pt_inlet = 1.05e5;
    T Tt_inlet = 288;
    T Ps_outlet = 1e5;
    // Init flow
    yams::GasPerfectGammaConstant<T> air{287.04, 1.4};
    T g = air.g(Tt_inlet, Pt_inlet);
    T M_outlet = sqrt( 2/(g-1)*(pow(Pt_inlet/Ps_outlet,(g-1)/g)-1));
    T Ts_outlet= Tt_inlet * (1+(g-1)/2*M_outlet*M_outlet);
    T u_outlet = M_outlet * air.a(Ts_outlet, Ps_outlet);
    T rho_outlet = air.rho(Ts_outlet, Ps_outlet);
    T et_outlet  = air.et(Ts_outlet, Ps_outlet, u_outlet);
    U.init(
        rho_outlet,
        rho_outlet*u_outlet,
        0.,
        et_outlet
    );
    // compute gog and volume

    vector<T> V(nm),X(n),Y(n),Xc(nm),Yc(nm);
    std::transform(std::execution::par, pts.begin(), pts.end(),X.begin(),[](const auto &pt){return pt[0];});
    std::transform(std::execution::par, pts.begin(), pts.end(),Y.begin(),[](const auto &pt){return pt[1];});
    
    {    
        auto start = std::chrono::steady_clock::now();
        for (long j{}; j < njm; j++)
        {
            for (long i{}; i < nim; i++)
            {
                auto x1 = X[i + ni * j];
                auto y1 = Y[i + ni * j];
                auto x2 = X[i + 1 + ni * j];
                auto y2 = Y[i + 1 + ni * j];
                auto x3 = X[i + 1 + ni * (j + 1)];
                auto y3 = Y[i + 1 + ni * (j + 1)];
                auto x4 = X[i + ni * (j + 1)];
                auto y4 = Y[i + ni * (j + 1)];
                V[i + nim * j]  = T(0.5) * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));
                Xc[i + nim * j] = T(0.25)* (x1 + x2 + x3 + x4);
                Yc[i + nim * j] = T(0.25)* (y1 + y2 + y3 + y4); 
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout<< "Cells volume and center, double loop: " << elapsed << "mcs" << std::endl;
    }
    
    {
        auto start = std::chrono::steady_clock::now();
        // #pragma omp for
        for( int id{}; id< nm; id++)
        {
            size_t i = id % nim;
            size_t j = (id-i) / nim;
            auto id1 = i + ni * j;
            auto id2 = i + 1 + ni * j;
            auto id3 = i + 1 + ni * (j + 1);
            auto id4 = i + ni * (j + 1);

            auto x1 = X[id1];
            auto y1 = Y[id1];
            auto x2 = X[id2];
            auto y2 = Y[id2];
            auto x3 = X[id3];
            auto y3 = Y[id3];
            auto x4 = X[id4];
            auto y4 = Y[id4];
            V[ id]  = T(0.5) * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));
            Xc[id] = T(0.25)* (x1 + x2 + x3 + x4);
            Yc[id] = T(0.25)* (y1 + y2 + y3 + y4); 
            
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells volume and center single loop:" << elapsed << "mcs" << std::endl;
    }

    gbs::points_vector<T,2> pts_c(nm);
    for( int id{}; id< nm; id++)
        pts_c[id] = {Xc[id],Yc[id]};
    // Compute faces areas

    vector<T> S_EW_x((ni)*(njm)),S_EW_y((ni)*(njm)),S_SN_x((nim)*(nj)),S_SN_y((nim)*(nj));
    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < ni; i++)
            {
                auto id1 = i + ni * j;
                auto id2 = i + 1 + ni * j;
                S_EW_x[id1] = Y[id2] - Y[id2];
                S_EW_y[id1] = X[id1] - X[id2];
            }
        }
        for(size_t i{}; i < nim; i++)
        {
            for(size_t j{}; j < nj; j++)
            {
                auto id1 = j + nj * i;
                auto id2 = j + 1 + nj * i;
                S_SN_x[id1] = Y[id1] - Y[id2];
                S_SN_y[id1] = X[id2] - X[id2];
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Cells Surf/Normals:" << elapsed << "mcs" << std::endl;
    }

    FlowData E{np};
    FlowData F{np};
    FlowData E_EW{ni*njm}, E_SN{nim*nj};
    FlowData F_EW{ni*njm}, F_SN{nim*nj};

    auto & [U1, U2, U3, U4] = U.d;
    auto & [E1, E2, E3, E4] = E.d;
    auto & [F1, F2, F3, F4] = F.d;
    vector<T> gamma(np, 1.4);
    {
        auto start = std::chrono::steady_clock::now();
// #pragma omp for
        for( int id{}; id<np ; id++)
        {
            auto u = U2[id]/U1[id];
            auto v = U3[id]/U1[id];
            auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

            E1[id] = U2[id];
            E2[id] = U2[id]/U1[id] + p;
            E3[id] = U2[id]*U3[id]/U1[id];
            E4[id] = (U4[id] + p) * u; 

            F1[id] = U3[id];
            F2[id] = U2[id]*U3[id]/U1[id];
            F3[id] = U3[id]/U1[id] + p;
            F4[id] = (U4[id] + p) * v; 
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F cells centers single loop:" << elapsed << "mcs" << std::endl;
    }

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njp; j++)
        {
            for( size_t id{nip*j}; id < nip * (j+1); id++)
            {
                // auto id = i + nip * j;
                auto u = U2[id]/U1[id];
                auto v = U3[id]/U1[id];
                auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

                E1[id] = U2[id];
                E2[id] = U2[id]/U1[id] + p;
                E3[id] = U2[id]*U3[id]/U1[id];
                E4[id] = (U4[id] + p) * u; 

                F1[id] = U3[id];
                F2[id] = U2[id]*U3[id]/U1[id];
                F3[id] = U3[id]/U1[id] + p;
                F4[id] = (U4[id] + p) * v; 
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F cells centers double loop:" << elapsed << "mcs" << std::endl;
    }

    // vector<size_t> ids_jp(njp);
    // std::generate(ids_jp.begin(), ids_jp.end(), [j=0]() mutable {return j++;});
    // {
    //     auto start = std::chrono::steady_clock::now();
    //     // size_t id{};
    //     std::for_each(
    //         // std::execution::par,
    //         ids_jp.begin(), ids_jp.end(),
    //         [&](auto j)
    //         {
    //             for( size_t id{nip * j}; id < nip * (j+1); id++)
    //             {
    //                 // auto id = i + nip * j;
    //                 auto u = U2[id]/U1[id];
    //                 auto v = U3[id]/U1[id];
    //                 auto p = (gamma[id]-1) * ( U4[id] - 0.5 * (u*u+v*v));

    //                 E1[id] = U2[id];
    //                 E2[id] = U2[id]/U1[id] + p;
    //                 E3[id] = U2[id]*U3[id]/U1[id];
    //                 E4[id] = (U4[id] + p) * u; 

    //                 F1[id] = U3[id];
    //                 F2[id] = U2[id]*U3[id]/U1[id];
    //                 F3[id] = U3[id]/U1[id] + p;
    //                 F4[id] = (U4[id] + p) * v; 
    //             }
    //         }
    //     );
    //     auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
    //     std::cout << "E/F cells centers for_each vectorized loop:" << elapsed << "mcs" << std::endl;
    // }

    auto & [E1_EW, E2_EW, E3_EW, E4_EW] = E_EW.d;
    auto & [F1_EW, F2_EW, F3_EW, F4_EW] = F_EW.d;

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < ni; i++)
            {
                auto id  = i + ni * j;
                auto id1p= i + nip * j;
                auto id2p= i+1 + nip * j;
                E1_EW[id] = 0.5 * ( E1[id1p] + E1[id2p]);
                F1_EW[id] = 0.5 * ( F1[id1p] + F1[id2p]);
                E2_EW[id] = 0.5 * ( E2[id1p] + E2[id2p]);
                F2_EW[id] = 0.5 * ( F2[id1p] + F2[id2p]);
                E3_EW[id] = 0.5 * ( E3[id1p] + E3[id2p]);
                F3_EW[id] = 0.5 * ( F3[id1p] + F3[id2p]);
                E4_EW[id] = 0.5 * ( E4[id1p] + E4[id2p]);
                F4_EW[id] = 0.5 * ( F4[id1p] + F4[id2p]);
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F_EW cells face double loop:" << elapsed << "mcs" << std::endl;
    }

    auto & [E1_SN, E2_SN, E3_SN, E4_SN] = E_SN.d;
    auto & [F1_SN, F2_SN, F3_SN, F4_SN] = F_SN.d;

    {
        auto start = std::chrono::steady_clock::now();
        for(size_t i{}; i < nim; i++)
        {
            for(size_t j{}; j < nj; j++)
            {
                auto id  = j + nj * i;
                auto id1p= i + nip * j;
                auto id2p= i + nip * (j+1);
                E1_SN[id] = 0.5 * ( E1[id1p] + E1[id2p]);
                F1_SN[id] = 0.5 * ( F1[id1p] + F1[id2p]);
                E2_SN[id] = 0.5 * ( E2[id1p] + E2[id2p]);
                F2_SN[id] = 0.5 * ( F2[id1p] + F2[id2p]);
                E3_SN[id] = 0.5 * ( E3[id1p] + E3[id2p]);
                F3_SN[id] = 0.5 * ( F3[id1p] + F3[id2p]);
                E4_SN[id] = 0.5 * ( E4[id1p] + E4[id2p]);
                F4_SN[id] = 0.5 * ( F4[id1p] + F4[id2p]);
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "E/F_SN cells faces double loop:" << elapsed << "mcs" << std::endl;
    }

    auto & [R1, R2, R3, R4] = R.d;
    {
        auto start = std::chrono::steady_clock::now();
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < nim; i++)
            {
                auto id1 = i + ni*j;
                auto id2 = j + nj*i;
                auto id3 = i+1 + ni*j;
                auto id4 = j+1 + nj *i;
                auto S_EW_x_id1=S_EW_x[id1];
                auto S_SN_x_id2=S_SN_x[id2];
                auto S_EW_x_id3=S_EW_x[id3];
                auto S_SN_x_id4=S_SN_x[id4];
                auto S_EW_y_id1=S_EW_y[id1];
                auto S_SN_y_id2=S_SN_y[id2];
                auto S_EW_y_id3=S_EW_y[id3];
                auto S_SN_y_id4=S_SN_y[id4];
                R1[i + j *ni] = -(
                    -S_EW_x_id1 * E1_EW[id1] - S_EW_y_id1 * F1_EW[id1]
                    -S_SN_x_id2 * E1_SN[id2] - S_SN_y_id2 * F1_SN[id2]
                    +S_EW_x_id3 * E1_EW[id3] + S_EW_y_id3 * F1_EW[id3]
                    +S_SN_x_id4 * E1_SN[id4] + S_SN_y_id4 * F1_SN[id4]
                );
            }
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
        std::cout << "Residual cells  double loop:" << elapsed << "mcs" << std::endl;
    }

    const std::array<T,5> rk_coefs{0.0533,0.1263,0.2375, 0.4414,1.0};



    plot(
        make_structuredgrid_actor(pts, ni, nj)
        , pts_c
    );
}

TEST(tests_euler2D, base_solver)
{
    size_t ni_{51}, nj_{25};
    T bump = 0.;
    auto [ pts, ni, nj, n_iso_ksi, n_iso_eth] = build_channel_mesh(ni_, nj_, bump);

    yams::Euler2DSolver<T> solver{pts, ni, nj}; 

    {
        T u = 123;
        T v = 65;
        T p = 1.53e5;
        T t = 364;
        auto U = solver.U({u, v, t, p}, 1.39);
        auto P = solver.P(U, 1.39);
        ASSERT_NEAR( std::get<0>(P), u, 1e-7);
        ASSERT_NEAR( std::get<1>(P), v, 1e-7);
        ASSERT_NEAR( std::get<2>(P), t, 1e-7);
        ASSERT_NEAR( std::get<3>(P), p, 1e-7);
}
    solver.init();
    solver.evalResidual();

}

