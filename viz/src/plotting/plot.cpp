#include <cstddef>
#include <filesystem>
#include <nane/viz/plotting/plot.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <vtkAxis.h>
#include <vtkChart.h>
#include <vtkChartMatrix.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkDoubleArray.h>
#include <vtkPNGWriter.h>
#include <vtkPlot.h>
#include <vtkPlotLine.h>
#include <vtkPlotPoints.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVector.h>
#include <vtkWindowToImageFilter.h>

namespace nane::viz
{
    struct plot::implementation
    {
        vtkSmartPointer<vtkContextView> view;
        vtkSmartPointer<vtkChartMatrix> matrix;
        vtkChart* active_chart = nullptr;
        std::vector<vtkSmartPointer<vtkTable>> tables;

        std::size_t rows = 1;
        std::size_t columns = 1;

        implementation() : view(vtkSmartPointer<vtkContextView>::New()), matrix(vtkSmartPointer<vtkChartMatrix>::New())
        {
            view->GetScene()->AddItem(matrix);
            view->GetRenderWindow()->SetWindowName("Nane");
            view->GetRenderWindow()->SetSize(900, 600);

            matrix->SetSize(vtkVector2i(1, 1));
            matrix->SetGutter(vtkVector2f(40.0F, 40.0F));

            active_chart = matrix->GetChart(vtkVector2i(0, 0));
            active_chart->SetShowLegend(true);
        }

        void select_subplot(std::size_t requested_rows, std::size_t requested_columns, std::size_t index)
        {
            if (requested_rows == 0 || requested_columns == 0)
                throw std::invalid_argument("subplot dimensions must be greater than zero");

            if (index >= requested_rows * requested_columns)
                throw std::invalid_argument("subplot index is outside the subplot grid");

            const bool layout_changed = requested_rows != rows || requested_columns != columns;

            if (layout_changed)
            {
                if (!tables.empty())
                    throw std::logic_error("subplot layout cannot be changed after adding data");

                rows = requested_rows;
                columns = requested_columns;
                matrix->SetSize(vtkVector2i(static_cast<int>(columns), static_cast<int>(rows)));
            }

            const std::size_t row_from_top = index / columns;
            const std::size_t column = index % columns;
            const std::size_t vtk_row = rows - 1 - row_from_top;

            active_chart = matrix->GetChart(vtkVector2i(static_cast<int>(column), static_cast<int>(vtk_row)));
            active_chart->SetShowLegend(true);
        }

        void add_series(std::span<const double> x, std::span<const double> y, std::string_view label, int plot_type)
        {
            if (x.size() != y.size())
                throw std::invalid_argument("plot coordinate sizes must match");

            auto table = vtkSmartPointer<vtkTable>::New();

            auto x_values = vtkSmartPointer<vtkDoubleArray>::New();
            auto y_values = vtkSmartPointer<vtkDoubleArray>::New();

            x_values->SetName("x");
            y_values->SetName("y");

            table->AddColumn(x_values);
            table->AddColumn(y_values);
            table->SetNumberOfRows(static_cast<vtkIdType>(x.size()));

            for (std::size_t i = 0; i < x.size(); ++i)
            {
                const auto index = static_cast<vtkIdType>(i);

                table->SetValue(index, 0, x[i]);
                table->SetValue(index, 1, y[i]);
            }

            vtkPlot* vtk_plot = active_chart->AddPlot(plot_type);

            vtk_plot->SetInputData(table, 0, 1);

            if (!label.empty())
                vtk_plot->SetLabel(std::string(label));

            if (plot_type == vtkChart::POINTS)
            {
                auto* points = vtkPlotPoints::SafeDownCast(vtk_plot);

                if (points != nullptr)
                {
                    points->SetMarkerStyle(vtkPlotPoints::PLUS);
                    points->SetMarkerSize(8.0);
                }
            }
            else if (plot_type == vtkChart::LINE)
            {
                auto* line = vtkPlotLine::SafeDownCast(vtk_plot);

                if (line != nullptr)
                    line->SetWidth(2.0);
            }

            tables.push_back(std::move(table));
        }

        void save(std::string_view filename)
        {
            if (filename.empty())
                throw std::invalid_argument("plot file name must not be empty");

            std::filesystem::path path{std::string(filename)};

            if (path.extension().empty())
                path.replace_extension(".png");

            else if (path.extension() != ".png")
                throw std::invalid_argument("plot save currently supports PNG files only");

            view->GetRenderWindow()->Render();

            auto image = vtkSmartPointer<vtkWindowToImageFilter>::New();

            image->SetInput(view->GetRenderWindow());
            image->SetInputBufferTypeToRGB();
            image->ReadFrontBufferOff();
            image->Update();

            auto writer = vtkSmartPointer<vtkPNGWriter>::New();

            const std::string output_path = path.string();

            writer->SetFileName(output_path.c_str());
            writer->SetInputConnection(image->GetOutputPort());
            writer->Write();
        }
    };

    plot::plot() : implementation_(std::make_unique<implementation>())
    {
    }

    plot::~plot() = default;

    plot::plot(plot&&) noexcept = default;

    plot& plot::operator=(plot&&) noexcept = default;

    void plot::subplot(std::size_t rows, std::size_t columns, std::size_t index)
    {
        implementation_->select_subplot(rows, columns, index);
    }

    void plot::add_points(std::span<const double> x, std::span<const double> y, std::string_view label)
    {
        implementation_->add_series(x, y, label, vtkChart::POINTS);
    }

    void plot::add_line(std::span<const double> x, std::span<const double> y, std::string_view label)
    {
        implementation_->add_series(x, y, label, vtkChart::LINE);
    }

    void plot::title(std::string_view value)
    {
        implementation_->active_chart->SetTitle(std::string(value));
    }

    void plot::x_label(std::string_view value)
    {
        implementation_->active_chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(std::string(value));
    }

    void plot::y_label(std::string_view value)
    {
        implementation_->active_chart->GetAxis(vtkAxis::LEFT)->SetTitle(std::string(value));
    }

    void plot::save(std::string_view filename)
    {
        implementation_->save(filename);
    }

    void plot::show()
    {
        implementation_->view->GetRenderWindow()->Render();
        implementation_->view->GetInteractor()->Initialize();
        implementation_->view->GetInteractor()->Start();
    }
} // namespace nane::viz
