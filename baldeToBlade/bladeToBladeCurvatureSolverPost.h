#pragma once
#include "bladeToBladeCurvatureSolver.h"
#include <stdexcept>
#include <string>
namespace yams
{
    template <typename T>
    using arr3 = std::array<T,3>;

    /// The number of computation planes from j1 to j2 included.
    /// @throw std::invalid_argument when j2 < j1 (the size would wrap around)
    inline size_t side_size(size_t j1, size_t j2)
    {
        if (j2 < j1)
            throw std::invalid_argument("yams: the trailing edge plane (" + std::to_string(j2) +
                                        ") is before the leading edge plane (" + std::to_string(j1) + ")");
        return j2 - j1 + 1;
    }

    template <typename T>
    auto getValueStreamLine(const vector<T> &val, size_t i_stream, size_t j1, size_t j2, size_t ni)
    {
        // j1 to j2 included
        vector<T> values(side_size(j1, j2));
        for(auto j{j1}; j <= j2; j++)
        {
            auto id = i_stream + ni * j;
            values[j-j1] = val[id];
        }
        return values; 
    }

    template <typename T>
    auto getValuesSide1(const BladeToBladeCurvatureSolver<T> &solver,const vector<T> &val)
    {
        return getValueStreamLine<T>(val, 0, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.dimensions()[0]);
    }

    template <typename T>
    auto getValuesSide2(const BladeToBladeCurvatureSolver<T> &solver,const vector<T> &val)
    {
        return getValueStreamLine<T>(val, solver.dimensions()[0]-1, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.dimensions()[0]);
    }

    template <typename T>
    auto getPsSide1(const BladeToBladeCurvatureSolver<T> &solver)
    {
        return getValueStreamLine<T>(solver.data().PS, 0, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.dimensions()[0]);
    }
        
    template <typename T>
    auto getPsSide2(const BladeToBladeCurvatureSolver<T> &solver)
    {
        return getValueStreamLine<T>(solver.data().PS, solver.dimensions()[0]-1, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.dimensions()[0]);
    }

    template <typename T>
    auto getXYZStreamLine(const BladeToBladeCurvatureSolverMesh<T> &msh, size_t i_stream, size_t j1, size_t j2, const gbs::Curve<T, 2> &meridional_stream_line)
    {
        using namespace std;
        // j1 to j2 included, like getValueStreamLine
        vector<arr3<T>> XYZ(side_size(j1, j2));
        auto ni = msh.ni;
        for(auto j{j1}; j <= j2; j++)
        {
            auto id     = i_stream + ni * j;
            auto [z, r] = meridional_stream_line(msh.U[id]);
            auto th     = msh.TH[id];
            XYZ[j - j1] = {
                r * cos(th),
                r * sin(th),
                z
            };
        }
        return XYZ; 
    }

    template <typename T>
    auto getXYZSide1(const BladeToBladeCurvatureSolver<T> &solver)
    {
        return getXYZStreamLine(solver.mesh(), 0, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.meridionalStreamLine());
    }

    template <typename T>
    auto getXYZSide2(const BladeToBladeCurvatureSolver<T> &solver)
    {
        auto [ni,nj] = solver.dimensions();
        return getXYZStreamLine(solver.mesh(), ni-1, solver.leadingEdgeIndex(), solver.trailingEdgeIndex(), solver.meridionalStreamLine());
    }
}