# pragma once
#include <array>
template <typename T>
auto make_uniform_grid(double l,double h, T &g,double a=0, double b=0,double ri=0.) -> void
{
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        auto l_ = l * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto h_ = h * j / (nj - 1.);
            g(i, j).x = l_;
            g(i, j).y = h_ + ri;
            if(i>0) g(i, j).y += tan(a) * (g(i, j).x - g(i-1, j).x);
            if(j>0) g(i, j).x += tan(a) * (g(i, j).y - g(i, j-1).y);
        }
    }
}

template <typename T>
auto make_uniform_clustered_grid(double l,double h, T &g,double a=0, double b=0,double ri=0.) -> void
{
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    const double PI_2 = acos(-1.) / 2.;
    for (auto i = 0; i < ni; i++)
    {
        auto l_ = l * sin( i / (ni - 1.) * PI_2);
        for (auto j = 0; j < nj; j++)
        {
            auto h_ = h * sin( j / (nj - 1.) * PI_2 );
            g(i, j).x = l_;
            g(i, j).y = h_ + ri;
            if(i>0) g(i, j).y += tan(a) * (g(i, j).x - g(i-1, j).x);
            if(j>0) g(i, j).x += tan(a) * (g(i, j).y - g(i, j-1).y);
        }
    }
}

template <typename T>
auto make_circular_grid(double r1,double r2,double t1,double t2,std::array<double,2> O, T &g) -> void
{
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        auto th = t1 + (t2 - t1) * i / (ni - 1.);
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            g(i, j).x = r * cos(th) + O[0];
            g(i, j).y = r * sin(th) + O[1];
        }
    }
}

template <typename T, typename G>
auto make_straight_igv(T r1, T r2, T l, T b1, T b2, G &g) -> void
{
    size_t ni = g.nRows();
    size_t nj = g.nCols();
    for (auto i = 0; i < ni; i++)
    {
        auto z = l * i / (ni - 1.);
        auto z1 = l / 3.;
        auto z2 = 2. * l / 3.;
        for (auto j = 0; j < nj; j++)
        {
            auto r = r1 + (r2 - r1) * j / (nj - 1.);
            g(i, j).x = z;
            g(i, j).y = r;
            if (z >= z1 && z <= z2)
            {
                g(i, j).bet = (z - z1) / (z2 - z1) * (b1 + (b2 - b1) * j / (nj - 1.));
                g(i, j).iB=0;
            }
        }
    }
}