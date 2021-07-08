#pragma once
#include <vtkNew.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <datastorage.h>
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
}