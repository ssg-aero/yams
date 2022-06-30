#pragma once

#include "datastorage.h"


namespace yams
{
    template <typename T,typename fX,typename fY>
    auto D1_O1_i_bw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  ( Y(g(i,j)) - Y(g(i-1,j)) ) / ( X(g(i,j)) - X(g(i-1,j)) );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O1_j_bw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  ( Y(g(i,j)) - Y(g(i,j-1)) ) / ( X(g(i,j)) - X(g(i,j-1)) );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O1_i_fw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  ( Y(g(i,j)) - Y(g(i+1,j)) ) / ( X(g(i,j)) - X(g(i+1,j)) );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O1_j_fw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  ( Y(g(i,j)) - Y(g(i,j+1)) ) / ( X(g(i,j)) - X(g(i,j+1)) );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_i_fw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return 1 / ( X(g(i+2,j)) - X(g(i+1,j)) ) * (
                ( X(g(i+2,j)) - X(g(i,j)) ) / ( X(g(i+1,j)) - X(g(i,j)) ) * ( Y(g(i+1,j)) - Y(g(i,j)) )
                -
                ( X(g(i+1,j)) - X(g(i,j)) ) / ( X(g(i+2,j)) - X(g(i,j)) ) * ( Y(g(i+2,j)) - Y(g(i,j)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_i_bw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  1 / ( X(g(i-1,j)) - X(g(i-2,j)) ) * (
                ( X(g(i,j)) - X(g(i-2,j)) ) / ( X(g(i,j)) - X(g(i-1,j)) ) * ( Y(g(i,j)) - Y(g(i-1,j)) )
                -
                ( X(g(i,j)) - X(g(i-1,j)) ) / ( X(g(i,j)) - X(g(i-2,j)) ) * ( Y(g(i,j)) - Y(g(i-2,j)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_i_ct( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return 1 / ( X(g(i+1,j)) - X(g(i-1,j)) ) * (
                ( X(g(i+1,j)) - X(g(i,j)) ) / ( X(g(i,j)) - X(g(i-1,j)) ) * ( Y(g(i,j)) - Y(g(i-1,j)) )
                +
                ( X(g(i,j)) - X(g(i-1,j)) ) / ( X(g(i+1,j)) - X(g(i,j)) ) * ( Y(g(i+1,j)) - Y(g(i,j)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_j_fw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return 1 / ( X(g(i,j+2)) - X(g(i,j+1)) ) * (
                ( X(g(i,j+2)) - X(g(i,j)) ) / ( X(g(i,j+1)) - X(g(i,j)) ) * ( Y(g(i,j+1)) - Y(g(i,j)) )
                -
                ( X(g(i,j+1)) - X(g(i,j)) ) / ( X(g(i,j+2)) - X(g(i,j)) ) * ( Y(g(i,j+2)) - Y(g(i,j)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_j_bw( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return  1 / ( X(g(i,j-1)) - X(g(i,j-2)) ) * (
                ( X(g(i,j)) - X(g(i,j-2)) ) / ( X(g(i,j)) - X(g(i,j-1)) ) * ( Y(g(i,j)) - Y(g(i,j-1)) )
                -
                ( X(g(i,j)) - X(g(i,j-1)) ) / ( X(g(i,j)) - X(g(i,j-2)) ) * ( Y(g(i,j)) - Y(g(i,j-2)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_j_ct( T &g,size_t i, size_t j,fY Y,fX X)
    {
        return 1 / ( X(g(i,j+1)) - X(g(i,j-1)) ) * (
                ( X(g(i,j+1)) - X(g(i,j)) ) / ( X(g(i,j)) - X(g(i,j-1)) ) * ( Y(g(i,j)) - Y(g(i,j-1)) )
                +
                ( X(g(i,j)) - X(g(i,j-1)) ) / ( X(g(i,j+1)) - X(g(i,j)) ) * ( Y(g(i,j+1)) - Y(g(i,j)) )
            );
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_i( T &g,size_t i, size_t j,fY Y,fX X)
    {
        if(i==0)
        {
            return D1_O2_i_fw(g,i,j,Y,X);
        }
        else if (i==g.nRows()-1)
        {
            return D1_O2_i_bw(g,i,j,Y,X);
        }
        else{
            return D1_O2_i_ct(g,i,j,Y,X);
        }
    }

    template <typename T,typename fX,typename fY>
    auto D1_O2_j( T &g,size_t i, size_t j,fY Y,fX X)
    {
        if(j==0)
        {
            return D1_O2_j_fw(g,i,j,Y,X);
        }
        else if (j==g.nCols()-1)
        {
            return D1_O2_j_bw(g,i,j,Y,X);
        }
        else{
            return D1_O2_j_ct(g,i,j,Y,X);
        }
    }
    //////////////// with metrics   ///////////////////
    template <typename T, typename fX,typename L>
    auto D1_O1_ksi_bw(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        return (-1. * X(g(i - 1, j)) + 1. * X(g(i , j)) ) /  d_ksi;
    }
    template <typename T, typename fX,typename L>
    auto D1_O1_eth_bw(T &g, size_t i, size_t j, const fX &X, L d_eth)
    {
        return (-1. * X(g(i, j - 1)) + 1. * X(g(i , j)) ) /  d_eth;
    }
    template <typename T, typename fX,typename L>
    auto D1_O1_ksi_fw(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        return ( 1. * X(g(i + 1, j)) - 1. * X(g(i , j)) ) /  d_ksi;
    }
    template <typename T, typename fX,typename L>
    auto D1_O1_eth_fw(T &g, size_t i, size_t j, const fX &X, L d_eth)
    {
        return ( 1. * X(g(i, j + 1)) - 1. * X(g(i , j)) ) /  d_eth;
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_fw(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        return (-3. * X(g(i, j)) + 4. * X(g(i + 1, j)) - 1. * X(g(i + 2, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_ct(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        return (-1. * X(g(i - 1, j)) + 1. * X(g(i + 1, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_bw(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        return (1. * X(g(i - 2, j)) - 4. * X(g(i - 1, j)) + 3. * X(g(i, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_fw(T &g, size_t i, size_t j, const fX &X, L d_eth)
    {
        return (-3. * X(g(i, j)) + 4. * X(g(i, j + 1)) - 1. * X(g(i, j + 2)) ) / (2. * d_eth);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_ct(T &g, size_t i, size_t j, const fX &X, L d_eth)
    {
        return (-1. * X(g(i, j - 1)) + 1. * X(g(i, j + 1)) ) / (2. * d_eth);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_bw(T &g, size_t i, size_t j, const fX &X, L d_eth)
    {
        return (1. * X(g(i, j - 2)) - 4. * X(g(i, j - 1)) + 3. * X(g(i, j)) ) / (2. * d_eth);
    }
    template <typename T, typename fX, typename L>
    auto D1_O2_ksi(T &g, size_t i, size_t j, const fX &X, L d_ksi)
    {
        if (i == 0)
        {
            return D1_O2_ksi_fw(g, i, j, X, d_ksi);
            // return D1_O1_ksi_fw(g, i, j, X, d_ksi);
        }
        else if (i == g.nRows() - 1)
        {
            return D1_O2_ksi_bw(g, i, j, X, d_ksi);
            // return D1_O1_ksi_bw(g, i, j, X, d_ksi);
        }
        else
        {
            return D1_O2_ksi_ct(g, i, j, X, d_ksi);
        }
    }

    template <typename T,typename fX,typename L>
    auto D1_O2_eth( T &g,size_t i, size_t j,const fX &X,L d_eth)
    {
        if(j==0)
        {
            return D1_O2_eth_fw(g,i,j,X,d_eth);
            // return D1_O1_eth_fw(g,i,j,X,d_eth);
        }
        else if (j==g.nCols()-1)
        {
            return D1_O2_eth_bw(g,i,j,X,d_eth);
            // return D1_O1_eth_bw(g,i,j,X,d_eth);
        }
        else{
            return D1_O2_eth_ct(g, i, j, X, d_eth);
        }
    }

    auto compute_metrics(const auto &g, const auto &fx1, const auto &fx2, auto &g_metrics)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        g_metrics.resize(ni, nj);
        double d_ksi = 1. / (ni - 1.);
        double d_eth = 1. / (nj - 1.);
        for (auto j = 0; j < nj; j++)
        {
            for (auto i = 0; i < ni; i++)
            {
                auto &gp_m = g_metrics(i, j);
                gp_m.x1_ksi = yams::D1_O2_ksi(g, i, j, fx1, d_ksi);
                gp_m.x1_eth = yams::D1_O2_eth(g, i, j, fx1, d_eth);
                gp_m.x2_ksi = yams::D1_O2_ksi(g, i, j, fx2, d_ksi);
                gp_m.x2_eth = yams::D1_O2_eth(g, i, j, fx2, d_eth);
                gp_m.J = 1. / (gp_m.x1_ksi * gp_m.x2_eth - gp_m.x2_ksi * gp_m.x1_eth);
            }
        }
    }

    auto D1_O2_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    }

    auto D1_O2_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    }

    auto D1_O2_ct_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_ct(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_ct(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    }

    auto D1_O2_ct_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_ct(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_ct(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    }

    auto D1_O1_bw_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O1_ksi_bw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O1_eth_bw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    }

    auto D1_O1_bw_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O1_ksi_bw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O1_eth_bw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    }

    auto D1_O2_bw_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_bw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_bw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    }

    auto D1_O2_bw_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_bw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_bw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    }

    auto D1_O2_fw_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_fw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_fw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    }

    auto D1_O2_fw_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, auto d_ksi, auto d_eth, const auto &fv_)
    {
        auto v_ksi = yams::D1_O2_ksi_fw(g, i, j, fv_, d_ksi);
        auto v_eth = yams::D1_O2_eth_fw(g, i, j, fv_, d_eth);
        auto gp_m = g_metrics(i, j);
        return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    }

    // template <typename T>
    // auto D1_O2_so_dx1(const auto &g,const auto &g_metrics,size_t i, size_t j, T d_ksi, T d_eth, const auto &fv_)
    // {
    //     T v_ksi, v_eth;
    //     if (i == 0)
    //     {
    //         v_ksi = 0.;
    //     }
    //     else if (i == 1)
    //     {
    //         v_ksi = yams::D1_O1_ksi_bw(g, i, j, fv_, d_ksi);
    //     }
    //     else
    //     {
    //         v_ksi = yams::D1_O2_ksi_bw(g, i, j, fv_, d_ksi);
    //     }
    //     if (j == 0)
    //     {
    //         v_eth = 0.;
    //     }
    //     else if (j == 1)
    //     {
    //         v_eth = yams::D1_O1_eth_bw(g, i, j, fv_, d_eth);
    //     }
    //     else
    //     {
    //         v_eth = yams::D1_O2_eth_bw(g, i, j, fv_, d_eth);
    //     }
    //     auto gp_m = g_metrics(i, j);
    //     return gp_m.J * gp_m.x2_eth * v_ksi - gp_m.J * gp_m.x2_ksi * v_eth;
    // }

    // template <typename T>
    // auto D1_O2_so_dx2(const auto &g,const auto &g_metrics,size_t i, size_t j, T d_ksi, T d_eth, const auto &fv_)
    // {
    //     T v_ksi, v_eth;
    //     if (i == 0)
    //     {
    //         v_ksi = 0.;
    //     }
    //     else if (i == 1)
    //     {
    //         v_ksi = yams::D1_O1_ksi_bw(g, i, j, fv_, d_ksi);
    //     }
    //     else
    //     {
    //         v_ksi = yams::D1_O2_ksi_bw(g, i, j, fv_, d_ksi);
    //     }
    //     if (j == 0)
    //     {
    //         v_eth = 0.;
    //     }
    //     else if (j == 1)
    //     {
    //         v_eth = yams::D1_O1_eth_bw(g, i, j, fv_, d_eth);
    //     }
    //     else
    //     {
    //         v_eth = yams::D1_O2_eth_bw(g, i, j, fv_, d_eth);
    //     }
    //     auto gp_m = g_metrics(i, j);
    //     return -gp_m.J * gp_m.x1_eth * v_ksi + gp_m.J * gp_m.x1_ksi * v_eth;
    // }
}