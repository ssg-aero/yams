#pragma once
#include <gbs-render/vtkcurvesrender.h>
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
#include <vtkTextProperty.h>
#include <vtkProperty2D.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <algorithm>
#include <meridionalsolvercase.h>
namespace yams
{
    template <typename T>
    inline auto make_vtkStructuredGrid(const MeridionalGrid<T> &g) -> vtkSmartPointer<vtkStructuredGrid>
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

    void add_grid_to_renderer(vtkRenderer* renderer, vtkStructuredGrid *structuredGrid,const char* name,bool edges_on, bool countour_on)
    {
        structuredGrid->GetPointData()->SetActiveScalars(name);

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

        vtkSmartPointer<vtkActor> gridActor =
            vtkSmartPointer<vtkActor>::New();
        gridActor->SetMapper(mapper);
        if(edges_on)
        {
            gridActor->GetProperty()->EdgeVisibilityOn();
            gridActor->GetProperty()->SetEdgeColor(0.3, 0.3, 0.3);
        }
        // Add the actor to the scene
        renderer->AddActor(gridActor);
        if (countour_on)
        { // Add contour
            vtkNew<vtkContourFilter> contourFilter;
            contourFilter->SetInputData(structuredGrid);
            double range[2];
            structuredGrid->GetScalarRange(range);
            contourFilter->GenerateValues(8, range[0], range[1]);
            // Map the contours to graphical primitives
            vtkNew<vtkPolyDataMapper> contourMapper;
            contourMapper->SetInputConnection(contourFilter->GetOutputPort());
            contourMapper->ScalarVisibilityOff();

            // Create an actor for the contours
            vtkNew<vtkActor> contourActor;
            contourActor->SetMapper(contourMapper);
            contourActor->GetProperty()->SetLineWidth(1.5);
            contourActor->GetProperty()->SetColor(0., 0., 0.);
            // contourActor->GetProperty()->SetOpacity(0.7);
            renderer->AddActor(contourActor);
        }

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
    }

    inline void plot_vtkStructuredGrid(vtkStructuredGrid *structuredGrid,const char* name,bool edges_on = false, bool countour_on = false)
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

        add_grid_to_renderer(renderer, structuredGrid, name, edges_on,countour_on);

        renderer->SetBackground(0.7, 0.7, 0.7); // Background color green

        renderWindowInteractor->SetInteractorStyle(style);

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
    }

    template<typename T>
    void add_span_to_renderer(vtkRenderer* renderer,const GridInfo<T> &gi, size_t i, const std::array<double,3> &col, double width)
    {
        gbs::points_vector<T,3> pts(gi.nj);
        const auto &g = *(gi.g);
        std::transform(
            g.begin(i),g.end(i),
            pts.begin(),
            [](const auto &gp){return gbs::point<T,3>{gp.x,gp.y,T(0.)};}
        );
        auto pts_ = gbs::make_vtkPoints(pts.begin(), pts.end());
        std::array<double,3> col_{col};
        auto actor = gbs::make_polyline_(pts_, col_ );
        actor->GetProperty()->SetLineWidth(width);
        renderer->AddActor( actor );
    }  

    template<typename T>
    void add_blades_to_renderer(vtkRenderer* renderer,const SolverCase<T> &solver_case)
    {

        for( const auto &bld_info: solver_case.bld_info_lst)
        {
            add_span_to_renderer(renderer, *(solver_case.gi), bld_info.i1, {1., 0., 0.}, 2.5);
            add_span_to_renderer(renderer, *(solver_case.gi), bld_info.i2, {0., 0., 1.}, 2.5);
        }
    }

    template<typename T>
    void plot_vtkStructuredGrid(const SolverCase<T> &solver_case,const char* name,bool edges_on = false, bool countour_on = false)
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

        auto structuredGrid = make_vtk_grid<T>(*(solver_case.gi->g));
        add_grid_to_renderer(renderer, structuredGrid, name, edges_on, countour_on);
        add_blades_to_renderer(renderer, solver_case);

        renderer->SetBackground(0.7, 0.7, 0.7); // Background color green

        renderWindowInteractor->SetInteractorStyle(style);

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
    }

    // #include <vtkAppendFilter.h>
    // // #include <vtkNew.h>
    // template<typename T>
    // void plot_SolverCaseSet(const SolverCaseSet<T> &set)
    // {
    //     // vtkNew<vtkAppendFilter> appendFilter;
    //     vtkSmartPointer<vtkAppendFilter> appendFilter =
    //         vtkSmartPointer<vtkAppendFilter>::New();
    //     size_t count{};
    //     for(const auto &solver_case : set.cases())
    //     {
    //         auto sgrid = make_vtk_grid<T>(*(solver_case->gi->g));
    //         appendFilter->AddInputData(  sgrid);
    //         count++;
    //     }
    // }
} // namespace yams