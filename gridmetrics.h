#pragma once
#include <datastorage.h>
namespace quiss
{
    template <typename T>
    inline auto distance(const MeridionalGridPoint<T> &gp1, const MeridionalGridPoint<T> &gp2) -> T
    {
        return sqrt((gp1.y - gp2.y) * (gp1.y - gp2.y) + (gp1.x - gp2.x) * (gp1.x - gp2.x));
    }

    template <typename T>
    inline auto compute_abscissas(MeridionalGrid<T> &g)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();
        for (auto i = 0; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i, j).m = i == 0 ? 0. : g(i - 1, j).m + distance(g(i, j), g(i - 1, j));
                g(i, j).l = j == 0 ? 0. : g(i, j - 1).l + distance(g(i, j), g(i, j - 1));
            }
        }
    }

    template <typename T>
    auto fz = [](MeridionalGridPoint<T> &gp) { return gp.x; };
    template <typename T>
    auto fr = [](MeridionalGridPoint<T> &gp) { return gp.y; };
    template <typename T>
    auto fm = [](MeridionalGridPoint<T> &gp) { return gp.m; };
    template <typename T>
    auto fl = [](MeridionalGridPoint<T> &gp) { return gp.l; };
    template <typename T>
    auto fphi = [](MeridionalGridPoint<T> &gp) { return gp.phi; };

    template <typename T>
    inline auto compute_angles(MeridionalGrid<T> &g)
    {
        size_t ni = g.nRows();
        size_t nj = g.nCols();

        T drqdm, dzqdm, drqdl, dzqdl;
        for (auto i = 0; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                // auto gam_ = -PI / 2 + PI * j / (nj - 1.);
                drqdm = D1_O2_i(g, i, j, fr<T>, fm<T>);
                dzqdm = D1_O2_i(g, i, j, fz<T>, fm<T>);
                g(i,j).phi = atan2(drqdm, dzqdm); // Stream line angle
                drqdl = D1_O2_j(g, i, j, fr<T>, fl<T>);
                dzqdl = D1_O2_j(g, i, j, fz<T>, fl<T>);
                g(i,j).gam = atan2(dzqdl, drqdl); // Span line angle
            }
        }
    }
    
    template <typename T>
    inline auto compute_curvature(MeridionalGrid<T> &g)
    {
        size_t ni = g.nRows()-1;
        size_t nj = g.nCols();

        for (auto i = 1; i < ni; i++)
        {
            for (auto j = 0; j < nj; j++)
            {
                g(i,j).cur = D1_O2_i(g, i, j, fphi<T>, fm<T>);// TODO check why Aungier put -DphiDm
            }
        }
        for (auto j = 0; j < nj; j++)
        {
            g(0, j).cur = 2. * g(1, j).cur - g(2, j).cur;
            g(ni, j).cur = 2. * g(ni - 1, j).cur - g(ni - 2, j).cur;
        }
    }

    template <typename T>
    inline auto compute_metrics(MeridionalGrid<T> &g)
    {
        compute_abscissas(g);
        compute_angles(g);
        compute_curvature(g);
    }
} // namespace quiss