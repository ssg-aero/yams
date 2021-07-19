#pragma once
#include <vtkNew.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <datastorage.h>
#include <numbers>

namespace quiss
{
    template<typename T>
    auto read_vtk_grid(const char *fname)
    {
        vtkNew<vtkXMLStructuredGridReader> reader;
        reader->SetFileName(fname);
        reader->Update();

        auto sgrid = reader->GetOutput();
        auto points= sgrid->GetPoints();
        auto dims  =sgrid->GetDimensions();
        size_t ni = dims[0];
        size_t nj = dims[1];

        MeridionalGrid<T> g(ni,nj);

        vtkIdType id {};
        auto iB_Array = sgrid->GetPointData()->GetAbstractArray("iB");
        auto k_Array = sgrid->GetPointData()->GetAbstractArray("k");
        for(size_t j {} ; j < nj ; j++)
        {
            for(size_t i {} ; i < ni ; i++)
            {
                auto pt = points->GetPoint(id);
                g(i,j).x = pt[0];
                g(i,j).y = pt[1];
                if(iB_Array) g(i,j).iB= iB_Array->GetVariantValue(id).ToInt();
                if(k_Array)  g(i,j).k= k_Array->GetVariantValue(id).ToDouble();
                id++;
            }
        }
        return g;

    }

    auto read_vtk_grid(auto &g, const char *fname)
    {
        vtkNew<vtkXMLStructuredGridReader> reader;
        reader->SetFileName(fname);
        reader->Update();

        auto sgrid = reader->GetOutput();
        auto points= sgrid->GetPoints();
        auto dims  =sgrid->GetDimensions();
        size_t ni = dims[0];
        size_t nj = dims[1];

        g.resize(ni,nj);

        vtkIdType id {};
        for(size_t j {} ; j < nj ; j++)
        {
            for(size_t i {} ; i < ni ; i++)
            {
                auto pt = points->GetPoint(id);
                g(i,j).x = pt[0];
                g(i,j).y = pt[1];
                id++;
            }
        }
        
    }

    template <typename T>
    auto write_vtk_grid(const MeridionalGrid<T> &g, const char *f_name) -> vtkSmartPointer<vtkStructuredGrid>
    {
            auto structuredGrid = make_vtkStructuredGrid(g);
            add_value(g, structuredGrid, "Vm", [](const auto &gp)
                      { return gp.Vm; });
            add_value(g, structuredGrid, "Vu", [](const auto &gp)
                      { return gp.Vu; });
            add_value(g, structuredGrid, "bet", [](const auto &gp)
                      { return gp.bet * 180 / std::numbers::pi; });
            add_value(g, structuredGrid, "Tt", [](const auto &gp)
                      { return gp.Tt; });
            add_value(g, structuredGrid, "Pt", [](const auto &gp)
                      { return gp.Pt; });
            add_value(g, structuredGrid, "Ts", [](const auto &gp)
                      { return gp.Ts; });
            add_value(g, structuredGrid, "Ps", [](const auto &gp)
                      { return gp.Ps; });
            add_value(g, structuredGrid, "s", [](const auto &gp)
                      { return gp.s; });
            add_value(g, structuredGrid, "rho", [](const auto &gp)
                      { return gp.rho; });
            add_value(g, structuredGrid, "H", [](const auto &gp)
                      { return gp.H; });
            add_value(g, structuredGrid, "q", [](const auto &gp)
                      { return gp.q; });
            add_value(g, structuredGrid, "cur", [](const auto &gp)
                      { return gp.cur; });
            add_value(g, structuredGrid, "gam", [](const auto &gp)
                      { return gp.gam; });
            add_value(g, structuredGrid, "phi", [](const auto &gp)
                      { return gp.phi; });
            add_value(g, structuredGrid, "r", [](const auto &gp)
                      { return gp.y; });
            add_value(g, structuredGrid, "z", [](const auto &gp)
                      { return gp.x; });
            add_value(g, structuredGrid, "cp", [](const auto &gp)
                      { return gp.Cp; });
            add_value(g, structuredGrid, "ga", [](const auto &gp)
                      { return gp.ga; });

            vtkNew<vtkXMLStructuredGridWriter> writer;
            writer->SetFileName(f_name);
            writer->SetInputData(structuredGrid);
            writer->Write();

            return structuredGrid;
    }
}