#pragma once
#include <mesh2d>
#include <geom.h>
#include <vector>
#include <vtkNew.h>
#include <vtkStructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>


template <typename T>
auto build_grid(
    const point2d<T> &O,
    T lx, T ly,
    size_t ni, size_t nj
    )
{
    std::vector<point2d<T>> grid_points(ni*nj);
    const T dx = lx / (ni - 1);
    const T dy = ly / (nj - 1);

    for (size_t j{}; j < nj; j++)
    {
        for (size_t i{}; i < ni; i++)
        {
            grid_points[i + ni * j] = {O.x + dx * i, O.y + dy * j};
        }
    }

    return grid_points;
}

template <typename T>
auto build_grid(
    const point2d<T> &O,
    T ri, T re,
    T a1, T d,
    size_t ni, size_t nj
    )
{
    std::vector<point2d<T>> grid_points(ni*nj);
    const T dr = (re-ri) / (nj - 1);
    const T da = d / (ni - 1);

    for (size_t j{}; j < nj; j++)
    {
        auto r = ri + j * dr;
        for (size_t i{}; i < ni; i++)
        {
            auto a = a1 + da * i;
            grid_points[i + ni * j] = {O.x + r * std::cos(a), O.y + r * std::sin(a)};
        }
    }

    return grid_points;
}

template <typename T>
auto build_vtk_sgrid(const std::vector<point2d<T>> &grid_points, size_t ni, size_t nj)
{
    const int dims[3] = {static_cast<int>(ni), static_cast<int>(nj), 1};

    // Create the structured grid.
    vtkNew<vtkStructuredGrid> sgrid;
    sgrid->SetDimensions(dims);
    vtkNew<vtkPoints> points;

    for (size_t j{}; j < nj; j++)
    {
            for (size_t i{}; i < ni; i++)
            {
                auto pt = grid_points[i + ni * j];
                points->InsertNextPoint(pt.x,pt.y, 0);
            }
    }
    sgrid->SetPoints(points);
    return sgrid;
}

template <typename T, typename U, typename F>
void add_values(vtkStructuredGrid *sgrid,const CellCenteredStructuredGrid2D<T,U> &msh, const F &getter, const char *name)
{
    vtkNew<vtkDoubleArray> t;
    t->SetName(name);
    auto ni = msh.ni();
    auto nj = msh.nj();
    for (size_t j{}; j < nj; j++)
    {
        for (size_t i{}; i < ni; i++)
        {
            t->InsertNextValue(getter(msh.values()[msh.id(i, j)]));
        }
    }
    sgrid->GetCellData()->AddArray(t);
}

template <typename T>
auto build_vtk_blockSet(const std::vector<std::vector<point2d<T>>> &grid_points_lst, const std::vector<size_t> &ni_lst, const std::vector<size_t> &nj_lst)
{
    vtkNew<vtkMultiBlockDataSet> set;
    int index{};
    for(const auto &grid_points : grid_points_lst)
    {
        auto ni = ni_lst[index];
        auto nj = nj_lst[index];
        set->SetBlock(index, build_vtk_sgrid(grid_points, ni, nj));
        index++;
    }
    return set;
}

void WriteMultiBlockData(vtkMultiBlockDataSet* mbds, const std::string& prefix)
{
    // Sanity check
    vtkNew<vtkXMLMultiBlockDataWriter> writer;

    std::ostringstream oss;
    oss.str("");
    oss.clear();
    oss << prefix << "." << writer->GetDefaultFileExtension();
    writer->SetFileName(oss.str().c_str());
    writer->SetInputData(mbds);
    writer->Write();
}