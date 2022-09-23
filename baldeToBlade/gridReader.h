#pragma once
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkDoubleArray.h>
#include <vtkProperty.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <gbs/gbslib.h>
#include "bladeToBladeCurvatureSolver.h"

namespace yams
{

    template <typename T>
    auto make_vtkStructuredGrid(const BladeToBladeCurvatureSolverMesh<T> &msh, T dth =0) -> vtkSmartPointer<vtkStructuredGrid>
    {
        // Build M, R*TH grid
        vtkSmartPointer<vtkStructuredGrid> structuredGrid =
            vtkSmartPointer<vtkStructuredGrid>::New();

        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

        auto nj = msh.nj;
        auto M = msh.M;
        auto R = msh.R;
        auto TH = msh.TH;
        auto n = M.size();
        for( size_t i{}; i<n; i++ )
        {
            points->InsertNextPoint(M[i], R[i]*( TH[i] + dth ), 0);
        }

        auto ni = M.size() / nj;
        structuredGrid->SetDimensions(ni, nj, 1);
        structuredGrid->SetPoints(points);

        return structuredGrid;
    }

    template <typename T>
    auto getVtkStructuredGridPoints(const char* f_name)
    {
        vtkNew<vtkXMLStructuredGridReader> reader;
        reader->SetFileName(f_name);
        reader->Update();

        vtkSmartPointer<vtkStructuredGrid> sGrid {reader->GetOutput()};

        auto points= sGrid->GetPoints();
        auto dims = sGrid->GetDimensions();
        size_t ni = dims[0];
        size_t nj = dims[1];

        gbs::points_vector<T,2> pts;

        for(int j {} ; j < nj ; j++)
        {
            auto jOffset = j * ni;
            for(int i {} ; i < ni ; i++)
            {
                auto offset = i + jOffset;
                auto X = points->GetPoint(offset);
                pts.push_back({X[0], X[1]});
            }
        }

        return std::make_pair(
            pts, nj
        );
    }

    template <typename T>
    void plot(const BladeToBladeCurvatureSolverMesh<T> &msh, const vector<T> &value, const char *name, bool edges_on = true, bool contour_on = true, T scale = 1.)
    {
        // Build grid with value
        auto nj = msh.nj;
        auto ni = msh.ni;
        auto structuredGrid     = yams::make_vtkStructuredGrid(msh);
        auto dth = msh.TH[ni-1]-msh.TH[0];
        auto structuredGrid_per = yams::make_vtkStructuredGrid(msh,dth);
        vtkSmartPointer<vtkDoubleArray> value_array =
            vtkSmartPointer<vtkDoubleArray>::New();
        value_array->SetName(name);

        std::for_each(
            value.begin(),
            value.end(),
            [scale, &value_array](auto v)
            {
                value_array->InsertNextValue(v * scale);
            });

        structuredGrid->GetPointData()->AddArray(value_array);
        structuredGrid_per->GetPointData()->AddArray(value_array);

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
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); // like paraview
        renderWindow->SetMultiSamples(4);

        renderer->SetBackground(0.7, 0.7, 0.7);

        renderWindowInteractor->SetInteractorStyle(style);

        // Create a mapper and actor
        structuredGrid->GetPointData()->SetActiveScalars(name);
        vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(structuredGrid);
        mapper->SetScalarRange(structuredGrid->GetScalarRange());

        structuredGrid_per->GetPointData()->SetActiveScalars(name);
        vtkSmartPointer<vtkDataSetMapper> mapper_per =
            vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_per->SetInputData(structuredGrid_per);
        mapper_per->SetScalarRange(structuredGrid_per->GetScalarRange());

        vtkSmartPointer<vtkLookupTable> rainbowBlueRedLut =
            vtkSmartPointer<vtkLookupTable>::New();
        rainbowBlueRedLut->SetNumberOfColors(256);
        rainbowBlueRedLut->SetHueRange(0.667, 0.0);
        rainbowBlueRedLut->Build();
        mapper->SetLookupTable(rainbowBlueRedLut);
        mapper_per->SetLookupTable(rainbowBlueRedLut);

        vtkSmartPointer<vtkActor> gridActor =
            vtkSmartPointer<vtkActor>::New();
        gridActor->SetMapper(mapper);
        vtkSmartPointer<vtkActor> gridActor_per =
            vtkSmartPointer<vtkActor>::New();
        gridActor_per->SetMapper(mapper_per);
        // add edges if required
        if (edges_on)
        {
            gridActor->GetProperty()->EdgeVisibilityOn();
            gridActor->GetProperty()->SetEdgeColor(0.9, 0.9, 0.9);
            gridActor_per->GetProperty()->EdgeVisibilityOn();
            gridActor_per->GetProperty()->SetEdgeColor(0.9, 0.9, 0.9);
        }
        // Add the actor to the scene
        renderer->AddActor(gridActor);
        renderer->AddActor(gridActor_per);
        if (contour_on)
        { // Add contour
            vtkNew<vtkContourFilter> contourFilter;
            contourFilter->SetInputData(structuredGrid);
            double range[2];
            structuredGrid->GetScalarRange(range);
            contourFilter->GenerateValues(15, range[0], range[1]);
            // Map the contours to graphical primitives
            vtkNew<vtkPolyDataMapper> contourMapper;
            contourMapper->SetInputConnection(contourFilter->GetOutputPort());
            contourMapper->ScalarVisibilityOff();

            // Create an actor for the contours
            vtkNew<vtkActor> contourActor;
            contourActor->SetMapper(contourMapper);
            contourActor->GetProperty()->SetLineWidth(1.5);
            contourActor->GetProperty()->SetColor(0., 0., 0.);
            contourActor->GetProperty()->SetOpacity(0.7);
            renderer->AddActor(contourActor);
/**************/
            vtkNew<vtkContourFilter> contourFilter_per;
            contourFilter_per->SetInputData(structuredGrid_per);
            contourFilter_per->GenerateValues(15, range[0], range[1]);
            // Map the contours to graphical primitives
            vtkNew<vtkPolyDataMapper> contourMapper_per;
            contourMapper_per->SetInputConnection(contourFilter_per->GetOutputPort());
            contourMapper_per->ScalarVisibilityOff();

            // Create an actor for the contours
            vtkNew<vtkActor> contourActor_per;
            contourActor_per->SetMapper(contourMapper_per);
            contourActor_per->GetProperty()->SetLineWidth(1.5);
            contourActor_per->GetProperty()->SetColor(0., 0., 0.);
            contourActor_per->GetProperty()->SetOpacity(0.7);
            renderer->AddActor(contourActor_per);
        }
        // Add scalar bar
        vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(mapper->GetLookupTable());
        scalarBar->SetTitle(name);
        scalarBar->SetNumberOfLabels(4);
        // scalarBar->GetTitleTextProperty()->SetColor(0.2, 0.2, 0.2);
        scalarBar->GetTitleTextProperty()->SetFontSize(24);
        // scalarBar->GetLabelTextProperty()->SetColor(0.2, 0.2, 0.2);
        scalarBar->GetLabelTextProperty()->SetFontSize(16);
        // scalarBar->SetDrawBackground(true);
        // scalarBar->GetBackgroundProperty()->SetOpacity(1.);
        // scalarBar->GetBackgroundProperty()->SetColor(0.3,0.3,0.3);
        scalarBar->SetUnconstrainedFontSize(true);
        renderer->AddActor(scalarBar);
        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
    }
}