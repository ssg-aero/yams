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
}