#pragma once 
#include <datastorage.h>
namespace quiss
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

    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_fw(T &g, size_t i, size_t j, fX X, L d_ksi)
    {
        return (-3. * X(g(i, j)) + 4. * X(g(i + 1, j)) - 1. * X(g(i + 2, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_ct(T &g, size_t i, size_t j, fX X, L d_ksi)
    {
        return (-1. * X(g(i - 1, j)) + 1. * X(g(i + 1, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_ksi_bw(T &g, size_t i, size_t j, fX X, L d_ksi)
    {
        return (1. * X(g(i - 2, j)) - 4. * X(g(i - 1, j)) + 3. * X(g(i, j)) ) / (2. * d_ksi);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_fw(T &g, size_t i, size_t j, fX X, L d_eth)
    {
        return (-3. * X(g(i, j)) + 4. * X(g(i, j + 1)) - 1. * X(g(i, j + 2)) ) / (2. * d_eth);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_ct(T &g, size_t i, size_t j, fX X, L d_eth)
    {
        return (-1. * X(g(i, j - 1)) + 1. * X(g(i, j + 1)) ) / (2. * d_eth);
    }
    template <typename T, typename fX,typename L>
    auto D1_O2_eth_bw(T &g, size_t i, size_t j, fX X, L d_eth)
    {
        return (1. * X(g(i, j - 2)) - 4. * X(g(i, j - 1)) + 3. * X(g(i, j)) ) / (2. * d_eth);
    }
    template <typename T, typename fX, typename L>
    auto D1_O2_ksi(T &g, size_t i, size_t j, fX X, L d_ksi)
    {
        if (i == 0)
        {
            return D1_O2_ksi_fw(g, i, j, X, d_ksi);
        }
        else if (i == g.nRows() - 1)
        {
            return D1_O2_ksi_bw(g, i, j, X, d_ksi);
        }
        else
        {
            return D1_O2_ksi_ct(g, i, j, X, d_ksi);
        }
    }

    template <typename T,typename fX,typename L>
    auto D1_O2_eth( T &g,size_t i, size_t j,fX X,L d_eth)
    {
        if(j==0)
        {
            return D1_O2_eth_fw(g,i,j,X,d_eth);
        }
        else if (j==g.nCols()-1)
        {
            return D1_O2_eth_bw(g,i,j,X,d_eth);
        }
        else{
            return D1_O2_eth_ct(g, i, j, X, d_eth);
        }
    }

    auto compute_metrics(const auto &g, const auto &fx, const auto &fy, auto &g_metrics)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        g_metrics.resize(ni, nj);
        double ksi = 1. / (ni - 1.);
        double eth = 1. / (nj - 1.);
        for (auto j = 0; j < nj; j++)
        {
            for (auto i = 0; i < ni; i++)
            {
                g_metrics(i, j).x_ksi = quiss::D1_O2_ksi(g, i, j, fx, ksi);
                g_metrics(i, j).x_eth = quiss::D1_O2_eth(g, i, j, fx, eth);
                g_metrics(i, j).y_ksi = quiss::D1_O2_ksi(g, i, j, fy, ksi);
                g_metrics(i, j).y_eth = quiss::D1_O2_eth(g, i, j, fy, eth);
                g_metrics(i, j).J = 1. / (g_metrics(i, j).x_ksi * g_metrics(i, j).y_eth - g_metrics(i, j).y_ksi * g_metrics(i, j).x_eth);
            }
        }
    }
}