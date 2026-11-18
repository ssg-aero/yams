#pragma once
#include <gbs/execution.h>

#include <algorithm>
#include <cstddef>
#include <exception>
#include <numeric>
#include <vector>

namespace yams::detail
{
    /// f(i) for i in [0, n), in parallel (GBS_PAR_EXEC). An exception escaping
    /// the function of a parallel algorithm calls std::terminate (MSVC: the
    /// process ends with 0xc0000409): each call is trapped, and the first
    /// failure (lowest i, as a serial loop would throw) is rethrown after. The
    /// other calls are all made, whether one fails or not.
    template <class F>
    void par_for(size_t n, F f)
    {
        std::vector<std::exception_ptr> errs(n);
        std::vector<size_t> idx(n);
        std::iota(idx.begin(), idx.end(), size_t{0});
        std::for_each(GBS_PAR_EXEC idx.begin(), idx.end(), [&](size_t i)
                      {
            try { f(i); }
            catch (...) { errs[i] = std::current_exception(); } });
        for (const auto &e : errs)
            if (e)
                std::rethrow_exception(e);
    }
}
