#pragma once
#include <vector>
namespace quiss
{
    template <typename Container,typename T>
    class _Array2d
    {
        Container container_;
        size_t nj_ = 0;

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
        size_t nRows() const noexcept
        {
            return container_.size() / nj_;
        }
        size_t nCols() const noexcept
        {
            return nj_;
        }
    };

    template <typename T,size_t ni,size_t nj>
    class Array2dStdArrayBased
    {
        std::array<T,ni*nj> container_;

    public:
        const T & operator() (size_t i, size_t j) const noexcept
        {
            return container_[j + nj * i];
        }
        T& operator() (size_t i, size_t j) noexcept
        {
            return container_[j + nj * i];
        }
        auto begin(size_t i)
        {
            return std::next(container_.begin(), nj * i);
        }
        auto end(size_t i)
        {
            return std::next(container_.begin(), nj - 1 + nj * i);
        }
        auto begin()
        {
            return container_.begin();
        }
        auto end()
        {
            return container_.end();
        }
        size_t nRows() const noexcept
        {
            return ni;
        }
        size_t nCols() const noexcept
        {
            return nj;
        }
    };

    template <typename T>
    using Array2d = _Array2d<std::vector<T>,T>;
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
    template <typename T,size_t ni,size_t nj>
    using Grid_ = Array2dStdArrayBased<GridPoint<T>,ni,nj>;

} // namespace quiss