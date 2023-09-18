#include "ParticleArea.h"

#include <Logging.h>
#include <vtkLookupTable.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>

#include "Libs/Optimize/Domain/VtkMeshWrapper.h"

namespace shapeworks {

//-----------------------------------------------------------------------------
void ParticleArea::assign_vertex_particles(vtkSmartPointer<vtkPolyData> poly_data,
                                           std::vector<itk::Point<double>> particles) {
  SW_DEBUG("Assigning vertex particles");
  // geodesics enabled mesh
  VtkMeshWrapper mesh(poly_data, true);

  // create "closest_particle" array
  auto closest_particle_array = vtkSmartPointer<vtkIntArray>::New();
  closest_particle_array->SetName("closest_particle");
  closest_particle_array->SetNumberOfComponents(1);
  closest_particle_array->SetNumberOfTuples(poly_data->GetNumberOfPoints());
  poly_data->GetPointData()->AddArray(closest_particle_array);

  // for each vertex
  for (int i = 0; i < poly_data->GetNumberOfPoints(); ++i) {
    // get the closest particle
    auto point = poly_data->GetPoint(i);

    // iterate over each particle to find the min
    int min_id = 0;
    double min_distance = std::numeric_limits<double>::max();
    for (int j = 0; j < particles.size(); ++j) {
      // get the geodesic distance
      auto distance = mesh.ComputeDistance(point, -1, particles[j], -1);
      if (distance < min_distance) {
        min_distance = distance;
        min_id = j;
      }
    }
    // assign the particle id
    closest_particle_array->SetTuple1(i, min_id);
  }
}

//-----------------------------------------------------------------------------
void ParticleArea::assign_vertex_colors(vtkSmartPointer<vtkPolyData> poly_data, std::vector<QColor> colors) {
  SW_DEBUG("Assigning vertex colors");
  // create rgb colors array
  auto colors_array = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors_array->SetNumberOfComponents(3);
  colors_array->SetName("colors");
  colors_array->SetNumberOfTuples(poly_data->GetNumberOfPoints());
  poly_data->GetPointData()->AddArray(colors_array);

  auto closest_particles = poly_data->GetPointData()->GetArray("closest_particle");

  // for each vertex
  for (int i = 0; i < poly_data->GetNumberOfPoints(); ++i) {
    // get the particle id from the "closest_point" array
    auto particle_id = closest_particles->GetTuple1(i);
    auto color = colors[particle_id];
    colors_array->SetTuple3(i, color.red(), color.green(), color.blue());
  }
  poly_data->GetPointData()->SetActiveScalars("colors");
}

//-----------------------------------------------------------------------------
std::vector<QColor> ParticleArea::colors_from_lut(vtkSmartPointer<vtkLookupTable> lut) {
  std::vector<QColor> colors;
  for (int i = 0; i < lut->GetNumberOfTableValues(); ++i) {
    double *color = lut->GetTableValue(i);
    colors.push_back(QColor(color[0] * 255, color[1] * 255, color[2] * 255));
  }
  return colors;
}

//-----------------------------------------------------------------------------
Eigen::VectorXd ParticleArea::compute_particle_areas(vtkSmartPointer<vtkPolyData> poly_data,
                                                     std::vector<itk::Point<double>> particles) {
  // The 'closest_particle' array should already be assigned and contains the particle id for each vertex
  auto closest_particles = poly_data->GetPointData()->GetArray("closest_particle");

  Eigen::VectorXd counts(particles.size());
  counts.setZero();
  // count how many vertices are assigned to each particle
  for (int i = 0; i < poly_data->GetNumberOfPoints(); ++i) {
    auto particle_id = closest_particles->GetTuple1(i);
    counts[particle_id] += 1;
  }

  // get total surface area using vtkMassProperties
  auto mass_properties = vtkSmartPointer<vtkMassProperties>::New();
  mass_properties->SetInputData(poly_data);
  mass_properties->Update();
  auto total_area = mass_properties->GetSurfaceArea();

  Eigen::VectorXd areas(particles.size());
  areas.setZero();

  // compute the area assigned to each particle
  for (int i = 0; i < particles.size(); ++i) {
    //    areas[i] = total_area * counts[i] / poly_data->GetNumberOfPoints();
    areas[i] = counts[i];
  }

  return areas;
}

//-----------------------------------------------------------------------------
Eigen::VectorXd ParticleArea::compute_particle_triangle_areas(vtkSmartPointer<vtkPolyData> poly_data,
                                                              std::vector<itk::Point<double>> particles) {
  auto closest_particles = poly_data->GetPointData()->GetArray("closest_particle");

  // for each cell
  Eigen::VectorXd areas(particles.size());
  areas.setZero();
  for (int i = 0; i < poly_data->GetNumberOfCells(); ++i) {
    // get the area of this cell
    auto cell = poly_data->GetCell(i);
    auto points = cell->GetPoints();
    double p0[3], p1[3], p2[3];
    points->GetPoint(0, p0);
    points->GetPoint(1, p1);
    points->GetPoint(2, p2);
    auto area = vtkTriangle::TriangleArea(p0, p1, p2);
  //    SW_TRACE(area);

    // for each vertex of the cell, give 1/3 of the area to the particle
    for (int j = 0; j < 3; ++j) {
      auto vertex_id = cell->GetPointId(j);
      auto particle_id = closest_particles->GetTuple1(vertex_id);
      areas[particle_id] += area / 3.0;
    }
  }

  return areas;
}

}  // namespace shapeworks
