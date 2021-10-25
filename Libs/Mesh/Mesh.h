#pragma once

#include "Shapeworks.h"
#include "ImageUtils.h"

class vtkCellLocator;
class vtkKdTreePointLocator;

namespace shapeworks {

class Mesh
{
public:
  enum AlignmentType { Rigid, Similarity, Affine };
  enum DistanceMethod { PointToPoint, PointToCell };
  enum CurvatureType { Principal, Gaussian, Mean };

  using MeshType = vtkSmartPointer<vtkPolyData>;

  Mesh(const std::string& pathname) : mesh(read(pathname)) {}
  Mesh(MeshType meshPtr) : mesh(meshPtr) { if (!mesh) throw std::invalid_argument("null meshPtr"); }
  Mesh(const Mesh& orig) : mesh(MeshType::New()) { mesh->DeepCopy(orig.mesh); }
  Mesh(Mesh&& orig) : mesh(orig.mesh) { orig.mesh = nullptr; }
  Mesh& operator=(const Mesh& orig) { mesh = MeshType::New(); mesh->DeepCopy(orig.mesh); return *this; }
  Mesh& operator=(Mesh&& orig) { mesh = orig.mesh; orig.mesh = nullptr; return *this; }

  /// append two meshes
  Mesh& operator+=(const Mesh& otherMesh);

  /// return the current mesh
  MeshType getVTKMesh() const { return this->mesh; }

  /// writes mesh, format specified by filename extension
  Mesh& write(const std::string &pathname);

  /// determines coverage between current mesh and another mesh (e.g. acetabular cup / femoral head)
  Mesh& coverage(const Mesh& otherMesh, bool allowBackIntersections = true,
                 double angleThreshold = 0, double backSearchRadius = 0);

  /// applies laplacian smoothing
  Mesh& smooth(int iterations = 0, double relaxation = 0.0);

  /// applies vtk windowed sinc smoothing
  Mesh& smoothSinc(int iterations = 0, double passband = 0.0);

  /// applies filter to reduce number of triangles in mesh
  Mesh& decimate(double reduction = 0.5, double angle = 15.0, bool preserveTopology = true);

  /// applies cvd (centroidal voronoi diagram) decimation filter
  Mesh& cvdDecimate(double percentage = 0.5);

  /// handle flipping normals
  Mesh& invertNormals();

  /// reflect meshes with respect to a specified center and specific axis
  Mesh& reflect(const Axis &axis, const Vector3 &origin = makeVector({ 0.0, 0.0, 0.0 }));

  /// creates a transform based on transform type
  MeshTransform createTransform(const Mesh &target, XFormType type = IterativeClosestPoint, AlignmentType align = Similarity, unsigned iterations = 10) const;

  /// applies the given transformation to the mesh
  Mesh& applyTransform(const MeshTransform transform);

  /// finds holes in a mesh and closes them
  Mesh& fillHoles();

  /// samples data values at specified point locations
  Mesh& probeVolume(const Image &image);

  /// clips a mesh using a cutting plane
  Mesh& clip(const Plane plane);

  /// helper to translate mesh
  Mesh& translate(const Vector3 &v);

  /// helper to scale mesh
  Mesh& scale(const Vector3 &v);

  /// computes bounding box of current mesh
  PhysicalRegion boundingBox() const;

  /// fix element winding of mesh
  Mesh& fixElement();

  /// computes distance from each vertex in this mesh to the closest location
  /// (PointToCell, the default) or closest vertex (PointToPoint) in target mesh
  Field distance(const Mesh &target, const DistanceMethod method = PointToCell) const;

  /// clips a mesh using a cutting plane resulting in a closed surface
  Mesh& clipClosedSurface(const Plane plane);

  /// computes and adds oriented point and cell normals
  Mesh& computeNormals();

  /// returns closest point on a face in the mesh to the given point in space
  /// - outside will be set to indicate if point is outside the mesh
  /// - face_id is the face on which this point is found
  /// non-const b/c determining if outside potentially requires computing cell normals
  Point3 closestPoint(const Point3 point, bool& outside, vtkIdType& face_id);

  /// returns closest point id in this mesh to the given point in space
  int closestPointId(const Point3 point) const;

  /// computes geodesic distance between two vertices (specified by their indices) on mesh
  double geodesicDistance(int source, int target) const;

  /// computes geodesic distance between a point (landmark) and each vertex on mesh
  Field geodesicDistance(const Point3 landmark) const;

  /// computes geodesic distance between a set of points (curve) and each vertex on mesh
  Field geodesicDistance(const std::vector<Point3> curve) const;

  /// computes curvature using principal (default) or gaussian or mean algorithms
  Field curvature(const CurvatureType type = Principal) const;

  /// rasterizes specified region to create binary image of desired dims (default: unit spacing)
  Image toImage(PhysicalRegion region = PhysicalRegion(), Point spacing = Point({1., 1., 1.})) const;

  /// converts specified region to distance transform image (default: unit spacing)
  Image toDistanceTransform(PhysicalRegion region = PhysicalRegion(), Point spacing = Point({1., 1., 1.})) const;

  // query functions //

  /// center of mesh
  Point3 center() const;

  /// center of mass of mesh
  Point3 centerOfMass() const;

  /// number of points
  int numPoints() const { return mesh->GetNumberOfPoints(); }

  /// number of faces
  int numFaces() const { return mesh->GetNumberOfCells(); }

  /// matrix with number of points with (x,y,z) coordinates of each point
  Eigen::MatrixXd points() const;

  /// matrix with number of faces with indices of the three points from which each face is composed
  Eigen::MatrixXi faces() const;

  /// (x,y,z) coordinates of vertex at given index
  Point3 getPoint(int id) const;

  /// return indices of the three points with which the face at the given index is composed
  IPoint3 getFace(int id) const;

  // fields of mesh points //

  /// print all field names in mesh
  std::vector<std::string> getFieldNames() const;

  /// sets the given field for points with array (*does not copy array's values)
  Mesh& setField(std::string name, Array array);

  /// sets the given field for faces with array (*does not copy array's values)
  Mesh& setFieldForFaces(std::string name, Array array);

  /// gets a pointer to the requested field, null if field doesn't exist
  Field getField(const std::string& name) const;

  /// gets a pointer to the requested field, null if field doesn't exist
  Field getFieldForFaces(const std::string& name) const;

  /// sets the given index of field to value
  void setFieldValue(const std::string& name, int idx, double value);

  /// gets the value at the given index of field (NOTE: returns first component of vector fields)
  double getFieldValue(const std::string& name, int idx) const;

  /// gets the multi value at the given index of [vertex] field
  Eigen::VectorXd getMultiFieldValue(const std::string& name, int idx) const;

  // mesh comparison //

  /// compare if values of the points in two (corresponding) meshes are (eps)equal
  bool compareAllPoints(const Mesh& other_mesh) const;

  /// compare if face indices in two (corresponding) meshes are equal
  bool compareAllFaces(const Mesh& other_mesh) const;

  /// compare if all fields in two meshes are (eps)equal
  bool compareAllFields(const Mesh& other_mesh, const double eps=-1.0) const;

  /// compare field of meshes to be (eps)equal (same field for both if only one specified)
  bool compareField(const Mesh& other_mesh, const std::string& name1, const std::string& name2="", const double eps=-1.0) const;

  // todo: add support for comparison of fields of mesh faces (ex: their normals)

  /// compare meshes
  bool compare(const Mesh& other_mesh, const double eps=-1.0) const;

  /// compare meshes
  bool operator==(const Mesh& other) const { return compare(other); }

  // public static functions //

  /// getSupportedTypes
  static std::vector<std::string> getSupportedTypes() { return {"vtk", "vtp", "ply", "stl", "obj"}; }

  /// Splits the mesh for FFCs by setting scalar and vector fields
  bool splitMesh(std::vector< std::vector< Eigen::Vector3d > > boundaries, Eigen::Vector3d query, size_t dom, size_t num);

  /// Gets values for FFCs
  double getFFCValue(Eigen::Vector3d query) const;

  /// Gets gradients for FFCs
  Eigen::Vector3d getFFCGradient(Eigen::Vector3d query) const;

  /// Formats mesh into an IGL format
  vtkSmartPointer<vtkPoints> getIGLMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const; // Copied directly from VtkMeshWrapper. this->poly_data_ becomes this->mesh. // WARNING: Copied directly from Meshwrapper. TODO: When refactoring, take this into account.

private:
  friend struct SharedCommandData;
  Mesh() : mesh(nullptr) {} // only for use by SharedCommandData since a Mesh should always be valid, never "empty"

  /// reads mesh (used only by constructor)
  static MeshType read(const std::string& pathname);

  /// Creates transform from source mesh to target using ICP registration
  MeshTransform createRegistrationTransform(const Mesh &target, AlignmentType align = Similarity, unsigned iterations = 10) const;

  MeshType mesh;

  /// invalidate cached cell and point locators (call when geometry changes)
  void invalidateLocators() const;

  /// Cell locator for functions that query for cells repeatedly
  // TODO: use vtkStaticCellLocator when vtk is upgraded to version 9
  mutable vtkSmartPointer<vtkCellLocator> cellLocator;
  void updateCellLocator() const;

  /// Point locator for functions that query for points repeatedly
  mutable vtkSmartPointer<vtkKdTreePointLocator> pointLocator;
  void updatePointLocator() const;

  /// Computes the gradient vector field for FFCs w.r.t the boundary
  std::vector<Eigen::Matrix3d> setGradientFieldForFFCs(vtkSmartPointer<vtkDoubleArray> absvalues, Eigen::MatrixXd V, Eigen::MatrixXi F) const;

  /// Computes scalar distance field w.r.t. the boundary
  vtkSmartPointer<vtkDoubleArray> setDistanceToBoundaryValueFieldForFFCs(vtkSmartPointer<vtkDoubleArray> values, vtkSmartPointer<vtkPoints> points, std::vector<size_t> boundaryVerts, vtkSmartPointer<vtkDoubleArray> inout, Eigen::MatrixXd V, Eigen::MatrixXi F, size_t dom);  // fixme: sets value, returns absvalues, not sure why

  /// Computes whether point is inside or outside the boundary
  vtkSmartPointer<vtkDoubleArray> computeInOutForFFCs(Eigen::Vector3d query, MeshType halfmesh); // similar issues to above

  /// Computes baricentric coordinates given a query point and a face number
  Eigen::Vector3d computeBarycentricCoordinates(const Eigen::Vector3d& pt, int face) const; // // WARNING: Copied directly from Meshwrapper. TODO: When refactoring, take this into account.

};

/// stream insertion operators for Mesh
std::ostream& operator<<(std::ostream &os, const Mesh& mesh);

} // shapeworks
