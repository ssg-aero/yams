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
            return std::next(container_.begin(), nj_ - 1 + nj_ * i);
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
    };

// xtensor specialization
    template <typename T>
    class _Array2d<xt::xarray<T>, T>
    {
        xt::xarray<T> container_;
        size_t nj_;
        size_t ni_;
        public:
        _Array2d(size_t ni, size_t nj);
        _Array2d(size_t ni, size_t nj, T v);
        const T &operator()(size_t i, size_t j) const noexcept {return container_(i,j);}
        T &operator()(size_t i, size_t j) noexcept {return container_(i,j);}
        auto begin(size_t i);
        auto end(size_t i);
        auto cbegin() const {return container_.cbegin();}
        auto cend() const {return container_.cend();}
        size_t nRows() const noexcept { return ni_;}
        size_t nCols() const noexcept { return nj_;}
    };

    template <typename T>
    _Array2d<xt::xarray<T>, T>::_Array2d(size_t ni, size_t nj)
    {
        xt::xarray<T>::shape_type shape = {ni, nj};
        container_ = xt::xarray<T>(shape);
        nj_ = nj;
        ni_ = ni;
    }
    template <typename T>
    _Array2d<xt::xarray<T>, T>::_Array2d(size_t ni, size_t nj, T v)
    {
        xt::xarray<T>::shape_type shape = {ni, nj};
        container_ = xt::xarray<T>(shape, v);
        nj_ = nj;
        ni_ = ni;
    }
// Aliases
    template <typename T>
    using Array2d = _Array2d<std::vector<T>,T>;
    template <typename T>
    using ArrayX2d = _Array2d<xt::xarray<T>,T>;
    // template <typename T,size_t ni,size_t nj>
    // using Array2d = Array2dStdArrayBased<T,ni,nj>;

    template <typename T>
    struct GridPoint
    {
        T x   = 0.;
        T y   = 0.;
        T l   = 0.;
        T m   = 0.;
        T phi = 0.;
        T gam = 0.;
        T bet = 0.;
        T cur = 0.;
    };


    template <typename T>
    using Grid = Array2d<GridPoint<T>>;
    template <typename T>
    using GridX = ArrayX2d<GridPoint<T>>;
    
    // template <typename Container1,typename T1,typename Container2,typename T2>
    // auto convert(_Array2d<Container1,T1> a) -> _Array2d<Container2,T2>
    // {
    //     _Array2d<Container2,T2> converted (a.nRows(),a.nCols());
    //     std::transform(
    //         std::execution::par,
    //         a.begin(),
    //         a.end(),
    //         converted.begin(),
    //         [](const auto &v_){return v_}
    //     )
    // }
} // namespace quiss