#pragma once
#include <vector>
#include <xtensor/xarray.hpp>
#include <execution>
namespace quiss
{
    template <typename Container,typename T>
    class _Array2d
    {
        size_t nj_ = 0;
        Container container_;

    public:
        _Array2d(size_t ni, size_t nj)
        {
            container_ = Container(ni * nj);
            nj_ = nj;
        }
        _Array2d(size_t ni, size_t nj, T v)
        {
            container_ = Container(ni * nj, v);
            nj_ = nj;
        }
        const auto &operator()(size_t i, size_t j) const noexcept
        {
            return container_[j + nj_ * i];
        }
        auto &operator()(size_t i, size_t j) noexcept
        {
            return container_[j + nj_ * i];
        }
        auto begin(size_t i)
        {
            return std::next(container_.begin(), nj_ * i);
        }
        auto end(size_t i)
        {
            return std::next(container_.begin(), nj_  + nj_ * i);
        }
        auto begin()
        {
            return container_.begin();
        }
        auto end()
        {
            return container_.end();
        }
        auto begin() const
        {
            return container_.begin();
        }
        auto end() const
        {
            return container_.end();
        }
        size_t nRows() const noexcept
        {
            return container_.size() / nj_;
        }
        size_t nCols() const noexcept
        {
            return nj_;
        }
        size_t size() const noexcept
        {
            return container_.size();
        }
    };

// xtensor specialization
    template <typename T>
    class _Array2d<xt::xarray<T>, T>
    {
        xt::xarray<T> container_;
        size_t nj_;
        public:
        _Array2d(size_t ni, size_t nj);
        _Array2d(size_t ni, size_t nj, T v);
        const T &operator()(size_t i, size_t j) const noexcept {return container_(i,j);}
        T &operator()(size_t i, size_t j) noexcept {return container_(i,j);}
        auto begin(size_t i) {return container_.begin() + (nj_ * i);}
        auto end(size_t i) {return container_.begin() + (nj_ - 1 + nj_ * i);}
        auto begin() const {return container_.cbegin();}
        auto end() const {return container_.cend();}
        auto begin() {return container_.begin();}
        auto end() {return container_.end();}
        size_t nRows() const noexcept { return container_.shape(0);}
        size_t nCols() const noexcept { return container_.shape(1);}
        size_t size() const noexcept  { return container_.size();  }
    };

    template <typename T>
    _Array2d<xt::xarray<T>, T>::_Array2d(size_t ni, size_t nj)
    {
        xt::xarray<T>::shape_type shape = {ni, nj};
        container_ = xt::xarray<T>(shape);
        nj_ = nj;
    }
    template <typename T>
    _Array2d<xt::xarray<T>, T>::_Array2d(size_t ni, size_t nj, T v)
    {
        xt::xarray<T>::shape_type shape = {ni, nj};
        container_ = xt::xarray<T>(shape, v);
        nj_ = nj;
    }
// Aliases
    template <typename T>
    using Array2d = _Array2d<std::vector<T>,T>;
    template <typename T>
    using ArrayX2d = _Array2d<xt::xarray<T>,T>;
    // template <typename T,size_t ni,size_t nj>
    // using Array2d = Array2dStdArrayBased<T,ni,nj>;

    template <typename T>
    struct MeridionalGridPoint
    {
        T x   = 0.;
        T y   = 0.;
        T l   = 0.;
        T m   = 0.; 
        T phi = 0.; // phi   = dz / dr
        T gam = 0.; // gamma = dr / dz 
        T bet = 0.; // atan2( Vu , Vm)
        T cur = 0.; // streamline curvature
        T cgp = 1.; // Cos( gamma + phi )
        T sgp = 0.; // Sin( gamma + phi )

        T Vm  = 0.;
        T Vu  = 0.;
        T rho = 1.225;
        T Pt  = 1.e5;
        T Tt  = 300;
        T Ps  = 1.e5;
        T Ts  = 300.;
        T Cp  = 1004.;
        T ga  = 1.4;
        T q   = 0.;
    };


    template <typename T>
    using MeridionalGrid = Array2d<MeridionalGridPoint<T>>;
    template <typename T>
    using MeridionalGridX = ArrayX2d<MeridionalGridPoint<T>>;

    template<typename T1,typename T2,template<typename> class S>
    auto copy(const S<T1> &a, S<T2> &b,size_t n) -> void
    {
        T1 *arrayT1 = (T1 *)&a;
        T2 *arrayT2 = (T2 *)&b;
        for (auto i = 0; i < n; i++)
        {
            arrayT2[i] = arrayT1[i];
        }
    }

    // template <typename T1,typename T2,template<typename> class S>
    // auto make_copy(const S<T1> &a,size_t n) -> S<T2>
    // {
    //     S<T2> b;
    //     copy(a,b,n);
    //     return b;
    // }

    template <typename Container1,typename T1,typename Container2,typename T2,template<typename> class S>
    auto copy(const _Array2d<Container1,S<T1>> &a, _Array2d<Container2,S<T2>> &b) -> void
    {
        if(a.size()==0) return;

        const size_t n =  sizeof(a(0,0)) / sizeof(T1);
        static_assert(sizeof(a(0,0)) == n*sizeof(T1));
        static_assert(sizeof(b(0,0)) == n*sizeof(T2));

        std::transform(
            std::execution::par,
            a.begin(),
            a.end(),
            b.begin(),
            b.begin(),
            [n](const auto &a_,const auto &b_){S<T2> bcp_ ;copy(a_,bcp_,n);return bcp_;}
        );
    }
} // namespace quiss