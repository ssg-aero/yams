#pragma once

#include <tuple>
#include <vector>
#include <algorithm>

namespace yams{
    //////// CHANGE FOR CUDA //////////////////
    template <typename T>
    using vector = std::vector<T>;
    template <typename... _Types>
    using tuple = std::tuple<_Types...>;
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
    // template <class _FwdIt, class _Ty>
    // const auto fill = std::fill<_FwdIt, _Ty>; // Aliasing didn't support type deduction
    template <class _FwdIt, class _Ty>
    void fill( _FwdIt beg, _FwdIt end, _Ty val)
    {
        std::fill(beg, end, val);
    }
    /////////////////////////////////////////


    template<typename T>
    struct SolverData4
    {
        size_t ni;
        size_t nj;
        tuple<vector<T>, vector<T>, vector<T>, vector<T>> d; // tuple can be used thanks to thrust in cuda
        SolverData4(size_t ni, size_t nj) : ni{ni}, nj{nj} {
            auto n = ni * nj;
            d = yams::make_tuple(
                vector<T>(n),
                vector<T>(n),
                vector<T>(n),
                vector<T>(n)
            );
        }
        SolverData4(size_t ni, size_t nj, T v1, T v2, T v3, T v4) : SolverData4{ni, nj} {
            init(v1, v2, v3, v4);
        }
        SolverData4(size_t ni, size_t nj, const tuple<T, T, T, T> &v) : SolverData4{ni, nj} {
            init(v);
        }
        void init( const tuple<T, T, T, T> &v )
        {
            auto [ v1, v2, v3, v4 ] = v;
            init(v1, v2, v3, v4);
        }
        void init(T v1, T v2, T v3, T v4)
        {
            // TODO find a way to alias for cuda thrust
            std::fill(std::get<0>(d).begin(), std::get<0>(d).end(), v1);
            std::fill(std::get<1>(d).begin(), std::get<1>(d).end(), v2);
            std::fill(std::get<2>(d).begin(), std::get<2>(d).end(), v3);
            std::fill(std::get<3>(d).begin(), std::get<3>(d).end(), v4);
        }
    };
}