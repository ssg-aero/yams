#include <vector>

template <typename T>
void thomas_algorithm(const std::vector<T> &a,
                      const std::vector<T> &b,
                      const std::vector<T> &c,
                      const std::vector<T> &d,
                      std::vector<T> &f)
{
    size_t N = d.size();

    std::vector<T> c_star(N, T{});
    std::vector<T> d_star(N, T{});

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < N; i++)
    {
        auto m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    for (int i = N - 1; i-- > 0;)
    {
        f[i] = d_star[i] - c_star[i] * d[i + 1];
    }
}
