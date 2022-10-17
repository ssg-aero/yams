#pragma once
#include <vector>
#include <gas/gasModel.h>

namespace yams{
    template <typename T>
    class Euler2DSolver
    {
        using vector = std::vector<T>;
        size_t ni;
        size_t nj;
        size_t njm;
        size_t nim;
        size_t njp;
        size_t nip;
        size_t nm; // N cells without ghosts
        size_t np; // N cells with ghosts
        size_t n;  // N vertices

        GasPerfectGammaConstant<T> gas{287.04, 1.4};
        T Pt_in{1.1e5};
        T Tt_in{340.};
        T Ps_out{1e5};

        vector V;
        vector Xc; // TODO suppress if not useful
        vector Yc; // TODO suppress if not useful
        vector X; //Storage useful only if moving mesh
        vector Y; //Storage useful only if moving mesh
        vector gamma;
        vector S_EW_x, S_EW_y, S_SN_x, S_SN_y;
        void computeCellsGeomInfo();
        virtual void computeCellsFlowInfo();
        virtual void computeCellFacesFlowInfo();
    protected:
        vector U1, U2, U3, U4;
        vector E1, E2, E3, E4;
        vector F1, F2, F3, F4;
        vector Q1, Q2, Q3, Q4;
        vector R1, R2, R3, R4;
        vector E1_EW, E2_EW, E3_EW, E4_EW;
        vector F1_EW, F2_EW, F3_EW, F4_EW;
        vector E1_SN, E2_SN, E3_SN, E4_SN;
        vector F1_SN, F2_SN, F3_SN, F4_SN;
    public:
        Euler2DSolver(std::vector<std::array<T,2>> &pts, size_t ni, size_t nj) :
            ni{ni}, nj{nj},
            nim{ni-1}, njm{nj-1},
            nip{ni+1}, njp{nj+1},
            nm{nim*njm}, n{ni*nj}, np{nip*njp},
            V(nm), Xc(nm), Yc(nm), 
            X(n), Y(n),
            S_EW_x((ni)*(njm)),
            S_EW_y((ni)*(njm)),
            S_SN_x((nim)*(nj)),
            S_SN_y((nim)*(nj)),
            U1(np), U2(np), U3(np), U4(np),
            E1(np), E2(np), E3(np), E4(np),
            F1(np), F2(np), F3(np), F4(np),
            Q1(np), Q2(np), Q3(np), Q4(np),
            R1(nm), R2(nm), R3(nm), R4(nm),
            E1_EW(ni*njm), E2_EW(ni*njm), E3_EW(ni*njm), E4_EW(ni*njm),
            F1_EW(ni*njm), F2_EW(ni*njm), F3_EW(ni*njm), F4_EW(ni*njm),
            E1_SN(nim*nj), E2_SN(nim*nj), E3_SN(nim*nj), E4_SN(nim*nj),
            F1_SN(nim*nj), F2_SN(nim*nj), F3_SN(nim*nj), F4_SN(nim*nj),
            gamma(np, 1.4)
        {
            // assert pts.size() == n;
            std::transform(std::execution::par, pts.begin(), pts.end(),X.begin(),[](const auto &pt){return pt[0];});
            std::transform(std::execution::par, pts.begin(), pts.end(),Y.begin(),[](const auto &pt){return pt[1];});
            computeCellsGeomInfo();
        }

        void init()
        {
            T gamma_loc = gas.g(Tt_in, Pt_in);
            T M_out = sqrt( 2/(gamma_loc-1)*(pow(Pt_in/Ps_out,(gamma_loc-1)/gamma_loc)-1));
            T Ts_out= Tt_in * (1+(gamma_loc-1)/2*M_out*M_out);
            T u_out = M_out * gas.a(Ts_out, Ps_out);
            T rho_out = gas.rho(Ts_out, Ps_out);
            T et_out  = gas.et(Ts_out, Ps_out, u_out);

            std::fill(U1.begin(), U1.end(), rho_out);
            std::fill(U2.begin(), U2.end(), rho_out*u_out);
            std::fill(U3.begin(), U3.end(), 0.);
            std::fill(U4.begin(), U4.end(), et_out);
            std::fill(gamma.begin(), gamma.end(), gamma_loc);

            computeCellsFlowInfo();
            computeCellFacesFlowInfo();
        }

        void evalResidual();

        auto sizes()
        {
            return std::make_pair(nim,njm);
        }
        auto idCell(size_t i, size_t j)
        {
            return i+nim*j;
        }
        auto idCellWithGhost(size_t i, size_t j)
        {
            return i+1+ni*(j+1);
        }
        auto R(size_t i, size_t j){
            auto id = idCell(i, j);
            return std::make_tuple(
                R1[id],
                R2[id],
                R3[id],
                R4[id]
            );
        }
        auto U(size_t i, size_t j){
            auto id = idCellWithGhost(i, j);
            return std::make_tuple(
                U1[id],
                U2[id],
                U3[id],
                U4[id]
            );
        }
        virtual auto U(const std::tuple<T,T,T,T> &tup, T g) -> std::tuple<T,T,T,T>
        {
            auto [ u, v, t, p ] = tup;
            auto u1 = p / gas.r() / t;
            auto u2 = u * u1;
            auto u3 = v * u1;
            auto u4 = p / (g-1) + 0.5 * u1 * (u*u+v*v);
            return std::make_tuple( u1, u2, u3, u4 );
        }
        auto P(size_t i, size_t j){
            auto id = idCellWithGhost(i, j);
            return P(std::make_tuple(
                U1[id],
                U2[id],
                U3[id],
                U4[id]
            ), gamma[id]);
        }
        virtual auto P(const std::tuple<T,T,T,T> &tup, T g) -> std::tuple<T,T,T,T>
        {
            auto [ u1, u2, u3, u4 ] = tup;
            auto ro=u1;
            auto u = u2/u1;
            auto v = u3/u1;
            auto p = (g-1) * ( u4 - 0.5 * ro*(u*u+v*v));
            auto t = p / gas.r() / ro;
            return std::make_tuple(
                u, v, t, p
            );
        }
    };


    template <typename T>
    void Euler2DSolver<T>::computeCellsGeomInfo()
    {
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

        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < ni; i++)
            {
                auto id  = i + ni * j;
                auto id1 = i + ni * j;
                auto id2 = i + ni * (j+1);
                S_EW_x[id] = Y[id2] - Y[id1];
                S_EW_y[id] = X[id1] - X[id2];
            }
        }

        for(size_t i{}; i < nim; i++)
        {
            for(size_t j{}; j < nj; j++)
            {
                auto id  = j + nj * i;
                auto id1 = i + ni * j;
                auto id2 = i+1 + ni * j;
                S_SN_x[id] = Y[id1] - Y[id2];
                S_SN_y[id] = X[id2] - X[id1];
            }
        }
    }
    
    template <typename T>
    void Euler2DSolver<T>::computeCellsFlowInfo()
    {
        for( size_t id{}; id<np ; id++)
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

            Q1[id] = 0.;
            Q2[id] = 0.;
            Q3[id] = 0.;
            Q4[id] = 0.;
        }
    }

    template <typename T>
    void Euler2DSolver<T>::computeCellFacesFlowInfo()
    {
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
    }

    template <typename T>
    void Euler2DSolver<T>::evalResidual()
    {
        for(size_t j{}; j < njm; j++)
        {
            for(size_t i{}; i < nim; i++)
            {
                auto id  = i + nim*j;
                auto id1 = i + ni*j;
                auto id2 = j + nj*i;
                auto id3 = i+1 + ni*j;
                auto id4 = j+1 + nj*i;

                auto S_EW_x_id1=S_EW_x[id1];
                auto S_SN_x_id2=S_SN_x[id2];
                auto S_EW_x_id3=S_EW_x[id3];
                auto S_SN_x_id4=S_SN_x[id4];

                auto S_EW_y_id1=S_EW_y[id1];
                auto S_SN_y_id2=S_SN_y[id2];
                auto S_EW_y_id3=S_EW_y[id3];
                auto S_SN_y_id4=S_SN_y[id4];

                R1[id] = -(
                    -S_EW_x_id1 * E1_EW[id1] - S_EW_y_id1 * F1_EW[id1]
                    -S_SN_x_id2 * E1_SN[id2] - S_SN_y_id2 * F1_SN[id2]
                    +S_EW_x_id3 * E1_EW[id3] + S_EW_y_id3 * F1_EW[id3]
                    +S_SN_x_id4 * E1_SN[id4] + S_SN_y_id4 * F1_SN[id4]
                );
                R2[id] = -(
                    -S_EW_x_id1 * E2_EW[id1] - S_EW_y_id1 * F2_EW[id1]
                    -S_SN_x_id2 * E2_SN[id2] - S_SN_y_id2 * F2_SN[id2]
                    +S_EW_x_id3 * E2_EW[id3] + S_EW_y_id3 * F2_EW[id3]
                    +S_SN_x_id4 * E2_SN[id4] + S_SN_y_id4 * F2_SN[id4]
                );
                R3[id] = -(
                    -S_EW_x_id1 * E3_EW[id1] - S_EW_y_id1 * F3_EW[id1]
                    -S_SN_x_id2 * E3_SN[id2] - S_SN_y_id2 * F3_SN[id2]
                    +S_EW_x_id3 * E3_EW[id3] + S_EW_y_id3 * F3_EW[id3]
                    +S_SN_x_id4 * E3_SN[id4] + S_SN_y_id4 * F3_SN[id4]
                );
                R4[id] = -(
                    -S_EW_x_id1 * E4_EW[id1] - S_EW_y_id1 * F4_EW[id1]
                    -S_SN_x_id2 * E4_SN[id2] - S_SN_y_id2 * F4_SN[id2]
                    +S_EW_x_id3 * E4_EW[id3] + S_EW_y_id3 * F4_EW[id3]
                    +S_SN_x_id4 * E4_SN[id4] + S_SN_y_id4 * F4_SN[id4]
                );
            }
        }
    }
}