#include "EigenUtils.h"

namespace shapeworks {

Eigen::MatrixXd itkTransformToEigen(TransformPtr itk_xform) {
  auto p = itk_xform->GetParameters();
  if (p.size() != 12) {
    std::ostringstream str;
    str << "itk transform is of incorrect dimensions (" << p.size() << ")";
    throw std::invalid_argument(str.str());
  }

  Eigen::MatrixXd eigen_xform = Eigen::MatrixXd::Identity(4, 4);
  for (size_t i = 0; i < p.size(); i++) {
    eigen_xform(i / 3, i % 3) = p[i];
  }

  return eigen_xform;
}

TransformPtr eigen44ToItkTransform(const Eigen::Matrix<double, 4, 4> &eigen_mat) {
  auto itk_xform = AffineTransform::New();

  auto p = itk::Array<double>(12);
  // the matrix
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      p[i * 3 + j] = eigen_mat(i, j);
    }
  }
  // the transform
  for (size_t j = 0; j < 3; j++) {
    p[3 * 3 + j] = eigen_mat(3, j);
  }
  itk_xform->SetParameters(p);

  return itk_xform;
}

using EigenMatrix = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;

EigenMatrix vtkTransformToEigen(MeshTransform vtk_xform) {
  EigenMatrix eigen = EigenMatrix::Identity(4, 4);
  eigen.block<4, 4>(0, 0) = Eigen::Map<EigenMatrix>(vtk_xform->GetMatrix()->GetData());
  return eigen;
}

MeshTransform eigen44ToVtkTransform(const EigenMatrix &eigen_mat) {
  vtkSmartPointer<vtkMatrix4x4> vtk_mat = vtkMatrix4x4::New();
  vtk_mat->DeepCopy(&eigen_mat(0, 0));
  MeshTransform xform = MeshTransform::New();
  xform->SetMatrix(vtk_mat);
  return xform;
}

std::pair<double, double> EigenUtils::compute_mean_and_stddev(const Eigen::MatrixXd &matrix, int row_index) {
  // Check if row_index is valid
  if (row_index < 0 || row_index >= matrix.rows()) {
    throw std::invalid_argument("Invalid row index");
  }

  // Extract the row as a vector
  Eigen::VectorXd row = matrix.row(row_index);

  // Calculate the mean of the row
  double mean = row.mean();

  // Calculate the sum of squared differences from the mean
  double sum_squared_diff = (row.array() - mean).square().sum();

  // Calculate the standard deviation
  double stddev = std::sqrt(sum_squared_diff / (matrix.cols() - 1));

  return std::make_pair(mean, stddev);
}

}  // namespace shapeworks
