#include <gtest/gtest.h>
#include <datastorage.h>
#include <diffop.h>
#include <gridsbuilders.h>
#include <gridmetrics.h>

#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

const double PI = acos(-1.);

using namespace quiss;

TEST(tests_visu, render_curvature)
{
    size_t ni = 50;
    size_t nj = 15;
    Grid<double> g(ni,nj);
    auto r1 =  1.;
    auto r2 =  2.;
    auto t1 = PI;
    auto t2 = PI * 2.;
    make_circular_grid(r1,r2,t1,t2,{0.,3.},g);
    compute_abscissas(g);
    compute_angles(g);

    // Create a grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
        vtkSmartPointer<vtkStructuredGrid>::New();

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> phi =
        vtkSmartPointer<vtkDoubleArray>::New();
    phi->SetName("Phi");
    vtkSmartPointer<vtkDoubleArray> gam =
        vtkSmartPointer<vtkDoubleArray>::New();
    gam->SetName("Gamma");

    std::for_each(
        // std::execution::par,
        g.begin(),
        g.end(),
        [&](const auto &gp){
            std::execution::par,
            points->InsertNextPoint(gp.x, gp.y, 0);
            gam->InsertNextValue(fmod(2*acos(-1)+gp.gam,2*acos(-1)));
            phi->InsertNextValue(gp.phi);
            }
    );

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(nj, ni, 1); // i j inverted Grid inner loop is i
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->AddArray(phi);
    structuredGrid->GetPointData()->AddArray(gam);
    structuredGrid->GetPointData()->SetActiveScalars("Phi");


    // Create a mapper and actor
    vtkSmartPointer<vtkDataSetMapper> mapper =
        vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(structuredGrid);
     mapper->SetScalarRange(structuredGrid->GetScalarRange());

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
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

    // Add the actor to the scene
    renderer->AddActor(actor);
    renderer->SetBackground(.3, .6, .3); // Background color green

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();
}