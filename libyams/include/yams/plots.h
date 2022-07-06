#pragma once

#include "curvaturesolver.h"

#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vtkAxis.h>

namespace yams
{
    template <typename T>
    auto plot_residual(const SolverLog<T> &log)
    {
        vtkNew<vtkNamedColors> colors;

        // Create a table with some points in it.
        vtkNew<vtkTable> table;

        vtkNew<vtkFloatArray> arrX;
        arrX->SetName("X Axis");
        table->AddColumn(arrX);

        vtkNew<vtkFloatArray> arr_res_moy;
        arr_res_moy->SetName("Arg Residual");
        table->AddColumn(arr_res_moy);

        vtkNew<vtkFloatArray> arr_res_max;
        arr_res_max->SetName("Max Residual");
        table->AddColumn(arr_res_max);

        // Fill in the table with some example values.
        int numPoints = log.delta_pos_moy.size();

        table->SetNumberOfRows(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            table->SetValue(i, 0, i );
            table->SetValue(i, 1, std::log10(log.delta_pos_moy[i]));
            table->SetValue(i, 2, std::log10(log.delta_pos_max[i]));
        }

        // Set up the view
        vtkNew<vtkContextView> view;
        view->GetRenderWindow()->SetWindowName("LinePlot");
        view->GetRenderer()->SetBackground(colors->GetColor3d("SlateGray").GetData());

        // Add multiple line plots, setting the colors etc.
        vtkNew<vtkChartXY> chart;
        view->GetScene()->AddItem(chart);
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetInputData(table, 0, 1);
        line->SetColor(0, 255, 0, 255);
        line->SetWidth(1.0);
        line = chart->AddPlot(vtkChart::LINE);
        line->SetInputData(table, 0, 2);
        line->SetColor(255, 0, 0, 255);
        // line->SetWidth(5.0);
        chart->SetShowLegend(true);
        chart->GetAxis(vtkAxis::LEFT)->SetTitle("Residual magnitude");
        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Iterations");
        // For dotted line, the line type can be from 2 to 5 for different dash/dot
        // patterns (see enum in vtkPen containing DASH_LINE, value 2):
        // #ifndef WIN32
          line->GetPen()->SetLineType(vtkPen::DASH_LINE);
        // #endif
        // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
        //  machines with built-in Intel HD graphics card...)

        // view->GetRenderWindow()->SetMultiSamples(0);

        // Start interactor
        view->GetRenderWindow()->Render();
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();
    }
}