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

}
