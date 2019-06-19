#ifndef IGL_GRID_H
#define IGL_GRID_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Construct vertices of a regular grid, suitable for input to
  // `igl::marching_cubes`
  //
  // Inputs:
  //   res  3-long list of number of vertices along the x y and z dimensions
  // Outputs:
  //   GV  res(0)*res(1)*res(2) by 3 list of mesh vertex positions.
  //   
  IGL_INLINE void grid(const Eigen::RowVector3i & res, Eigen::MatrixXd & GV);
}
#ifndef IGL_STATIC_LIBRARY
#  include "grid.cpp"
#endif
#endif 
