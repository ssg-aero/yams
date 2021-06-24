#pragma once
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>

namespace quiss
{
    template <typename T>
    inline auto make_vtkStructuredGrid(const MeridionalGrid<T> &g)
    {
        vtkSmartPointer<vtkStructuredGrid> structuredGrid =
            vtkSmartPointer<vtkStructuredGrid>::New();

        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

        std::for_each(
            // std::execution::par,
            g.begin(),
            g.end(),
            [&](const auto &gp) {
                points->InsertNextPoint(gp.x, gp.y, 0);
            });
        auto ni = g.nRows();
        auto nj = g.nCols();
        structuredGrid->SetDimensions(nj, ni, 1); // i j inverted Grid inner loop is i
        structuredGrid->SetPoints(points);
        return structuredGrid;
    }

    template <typename T, typename _Func>
    inline auto add_value(const MeridionalGrid<T> &g, vtkStructuredGrid *structuredGrid, const char *name, _Func f)
    {
        vtkSmartPointer<vtkDoubleArray> phi =
            vtkSmartPointer<vtkDoubleArray>::New();

        phi->SetName(name);
        // phi->Allocate(g.size());

        std::for_each(
            // std::execution::par,
            g.begin(),
            g.end(),
            [&](const auto &gp) {
                phi->InsertNextValue(f(gp));
            });

        structuredGrid->GetPointData()->AddArray(phi);
    }

    inline void plot_vtkStructuredGrid(vtkStructuredGrid *structuredGrid,bool edges_on=false)
    {
        // Create a mapper and actor
        vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(structuredGrid);
        mapper->SetScalarRange(structuredGrid->GetScalarRange());

        vtkSmartPointer<vtkLookupTable> rainbowBlueRedLut =
            vtkSmartPointer<vtkLookupTable>::New();
        rainbowBlueRedLut->SetNumberOfColors(256);
        rainbowBlueRedLut->SetHueRange(0.667, 0.0);
        rainbowBlueRedLut->Build();
        mapper->SetLookupTable(rainbowBlueRedLut);

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        if(edges_on)
            actor->GetProperty()->EdgeVisibilityOn();

        // Create a renderer, render window, and interactor
        vtkSmartPointer<vtkRenderer> renderer =
            vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
                vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview


        // Add the actor to the scene
        renderer->AddActor(actor);
        vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(mapper->GetLookupTable());
        scalarBar->SetTitle(structuredGrid->GetPointData()->GetArrayName(0));
        scalarBar->SetNumberOfLabels(4);
        renderer->AddActor(scalarBar);
        renderer->SetBackground(.3, .6, .3); // Background color green

        renderWindowInteractor->SetInteractorStyle(style);

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
    }
} // namespace quiss