#pragma once

template <typename T>
struct point2d
{
    T x;
    T y;
};

template <typename T>
point2d<T> operator+(const point2d<T> &p1,const point2d<T> &p2)
{
    return {
        p1.x+p2.x,
        p1.y+p2.y
    };
}

template <typename T>
point2d<T> operator-(const point2d<T> &p1,const point2d<T> &p2)
{
    return {
        p1.x-p2.x,
        p1.y-p2.y
    };
}

template <typename T>
point2d<T> operator*(T s,const point2d<T> &p)
{
    return {
        s*p.x,
        s*p.y
    };  
}

template <typename T>
point2d<T> rot90(const point2d<T> &p)
{
    return {-p.y, p.x};
}


template <typename T>
T sq_norm(const point2d<T> &p)
{
    return p.x * p.x + p.y * p.y;
}

template <typename T>
T sq_norm(const point2d<T> &p1,const point2d<T> &p2)
{
    return sq_norm(p2-p1);
}