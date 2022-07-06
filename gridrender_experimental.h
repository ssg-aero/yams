#pragma once
#include "meridionalsolvercase_experiemental.h"

namespace yams
{
    template<typename T>
    void plot_vtkStructuredGrid(const SolverCaseSet<T> &set,const char* name,bool edges_on = false, bool countour_on = false)
    {

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

        std::for_each(
            set.cases().begin(), set.cases().end(),
            [&renderer, name, edges_on, countour_on](const auto &p_solver_case)
            {
                const auto &solver_case = *p_solver_case;
                auto structuredGrid = make_vtk_grid<T>(*(solver_case.gi->g));
                add_grid_to_renderer(renderer, structuredGrid, name, edges_on, countour_on);
                add_blades_to_renderer(renderer, solver_case);
            }
        );

        renderer->SetBackground(0.7, 0.7, 0.7); // Background color green

        renderWindowInteractor->SetInteractorStyle(style);

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
    }
}