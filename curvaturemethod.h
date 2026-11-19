#pragma once
// How a streamline curvature solver computes the curvature of its stream lines
// from their points: shared by the meridional solver (compute_curvature,
// GridInfo::curvature) and the blade-to-blade one (BladeToBladeCurvatureSolver::
// setCurvature), so that both use the same methods.
#include <stdexcept>
#include <utility>

namespace yams
{
    /// How the stream lines' curvature is computed at a grid point, from the
    /// stream line's points at the stations (meridional solver, compute_curvature)
    /// or planes (blade-to-blade solver, BladeToBladeCurvatureSolver::setCurvature)
    /// around it. The descriptions below, and their stability analysis
    /// (max_relax_constant), are the meridional solver's; the blade-to-blade one
    /// takes a three-point second derivative of theta along its stream lines.
    enum class CurvatureMethod
    {
        /// The derivative along the stream line of its angle, itself a centred
        /// derivative of its points: the most stable (Wilkinson's k_min = -1), but
        /// blind to a zigzag from one station to the next, which can then stay
        /// in a converged solution (waves on the stream lines downstream of a
        /// leading edge).
        DoubleDifferentiation,
        /// A weighted mean of DoubleDifferentiation and ThreePoint, the share w of
        /// ThreePoint being GridInfo::curvature_blend: it sees the zigzag as soon
        /// as w > 0 (k = -4 w for it), at a small cost in stability while w is
        /// small: k_min = -1 / (1 - w) up to w = 1/2 (-2, half as stable as
        /// DoubleDifferentiation), -4 w beyond.
        Blended,
        /// The curvature of the circle through the point and its neighbours on
        /// the stream line: the most accurate, k_min = -4, a quarter as stable.
        ThreePoint,
    };

    /// The share of ThreePoint in a curvature method: 0 for DoubleDifferentiation,
    /// 1 for ThreePoint, blend for Blended (GridInfo::curvature_blend). The one
    /// place this mapping is made, for compute_curvature and max_relax_constant.
    /// @throw std::invalid_argument for Blended with a share outside [0, 1] (or NaN)
    constexpr double three_point_share(CurvatureMethod method, double blend)
    {
        switch (method)
        {
        case CurvatureMethod::DoubleDifferentiation: return 0.;
        case CurvatureMethod::ThreePoint: return 1.;
        case CurvatureMethod::Blended:
            if (!(blend >= 0. && blend <= 1.))
                throw std::invalid_argument("yams: curvature_blend must be between 0 and 1");
            return blend;
        }
        std::unreachable();
    }

    /// The largest relaxation constant (GridInfo::RF_B) a curvature method is
    /// stable with: Wilkinson's 96 / (5 |k_min|), the optimum of his analysis and
    /// about the limit found on an annulus and on NASA Rotor 37 (16 to 20, 8 to
    /// 9, 4 to 5 from 31 to 81 stations for DoubleDifferentiation, Blended with a
    /// share of 1/2 and ThreePoint).
    /// @param blend the share of ThreePoint in Blended (GridInfo::curvature_blend)
    /// @throw std::invalid_argument for Blended with a share outside [0, 1]
    constexpr double max_relax_constant(CurvatureMethod method, double blend = 0.5)
    {
        // k(s) = -4 s (1 - s (1 - w)), s = sin^2(theta / 2) in [0, 1]: its minimum is
        // at s = 1 / (2 (1 - w)) when w <= 1/2, else at s = 1 (the zigzag)
        const double w = three_point_share(method, blend);
        const double k_min = w <= 0.5 ? 1. / (1. - w) : 4. * w;
        return 96. / (5. * k_min);
    }

    /// The relaxation constant to use with a curvature method: 0.625 of
    /// max_relax_constant (12, 6 for Blended with a share of 1/2, 3), for a margin
    /// on other geometries.
    constexpr double default_relax_constant(CurvatureMethod method, double blend = 0.5)
    {
        return 0.625 * max_relax_constant(method, blend);
    }
}
