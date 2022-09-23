#include <vector>

template <typename T>
void thomas_algorithm(const std::vector<T> &a,
                      const std::vector<T> &b,
                      std::vector<T> &c,
                      std::vector<T> &d)
{
    size_t n = d.size() - 1;

    c[0] = c[0] / b[0];
    d[0] = d[0] / b[0];

    for (size_t i = 1; i < n; i++)
    {
        auto m = (b[i] - a[i] * c[i - 1]);
        c[i] = c[i] / m;
        d[i] = (d[i] - a[i] * d[i - 1]) / m;
    }
    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);
    for (size_t i = n; i-- > 0;)
        d[i] -= c[i] * d[i + 1];
    // d[0] -= c[0] * d[1];
}

    // template <typename T>
    // void thomas_algorithm(const std::vector<T> &A,
    //                       const std::vector<T> &B,
    //                       const std::vector<T> &C,
    //                       const std::vector<T> &D,
    //                       std::vector<T> &F)
    // {
    //     size_t N = D.size();

    //     std::vector<T> c_star(N, T{});
    //     std::vector<T> d_star(N, T{});

    //     c_star[0] = c[0] / b[0];
    //     d_star[0] = d[0] / b[0];

    //     for (int i = 1; i < N; i++)
    //     {
    //         auto m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
    //         c_star[i] = c[i] * m;
    //         d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    //     }

    //     for (int i = N - 1; i-- > 0;)
    //     {
    //         f[i] = d_star[i] - c_star[i] * d[i + 1];
    //     }
    // }
