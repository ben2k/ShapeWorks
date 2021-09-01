#include <Eigen/Eigen>

#include <Libs/Optimize/Optimize.h>

Eigen::MatrixXd optimize_get_particle_system(shapeworks::Optimize *opt)
{
  shapeworks::MatrixContainer container = opt->GetParticleSystem();
  return container.matrix_;
}

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace pybind11::literals;

#include <sstream>
#include <itkImportImageFilter.h>
#include <vtkDoubleArray.h>

#include "Shapeworks.h"
#include "ShapeworksUtils.h"
#include "Image.h"
#include "VectorImage.h"
#include "ImageUtils.h"
#include "Mesh.h"
#include "MeshUtils.h"
#include "Optimize.h"
#include "ParticleSystem.h"
#include "ShapeEvaluation.h"
#include "ParticleShapeStatistics.h"
#include "EigenUtils.h"

using namespace shapeworks;

PYBIND11_MODULE(shapeworks_py, m)
{
  m.doc() = "ShapeWorks Python API";

  m.attr("Pi") = std::atan(1.0) * 4.0;

  // Transform
  py::class_<itk::SmartPointer<itk::Transform<double, 3u, 3u> >>(m, "Transform")
  .def("__repr__",
       [](const TransformPtr &transform) {
         std::stringstream stream;
         itk::Transform<double, 3, 3>::ParametersType p = transform->GetParameters();
         int r = 0;
         for (int i=0; i<4; i++)
         {
           for (int j=0; j<3; j++)
           {
             stream << p[r] << " ";
             r++;
           }
         }
         return stream.str();
       });

  // Constructs an itk::AffineTransform from the 3x3 scale-rotate-warp and 3x1 translation.
  m.def("createTransform",
        [](Eigen::Matrix<double, 3, 3> &mat, std::vector<double> v) -> decltype(auto) {
          Matrix33 mat33 = eigenToItk<double, 3, 3>(mat);
          return createTransform(mat33, makeVector({v[0], v[1], v[2]}));
        },
        "creates transform from 3x3 matrix and translation vector",
        "mat"_a, "translate"_a=std::vector<double>({0,0,0}));

  // Axis
  py::enum_<Axis>(m, "Axis")
  .value("invalid", Axis::invalid)
  .value("X", Axis::X)
  .value("Y", Axis::Y)
  .value("Z", Axis::Z)
  .export_values();

  m.def("axis_is_valid",
        [](const std::vector<double> &axis) -> decltype(auto) {
          return axis_is_valid(makeVector({axis[0], axis[1], axis[2]}));
        },
        "ensure an axis is valid",
        "axis"_a);

  m.def("degToRad",
        degToRad,
        "convert degrees to radians",
        "deg"_a);

  m.def("toAxis",
        toAxis,
        "convert to axis",
        "str"_a);

  // TransformType
  py::enum_<XFormType>(m, "TransformType")
  .value("CenterOfMass", XFormType::CenterOfMass)
  .value("IterativeClosestPoint", XFormType::IterativeClosestPoint)
  .export_values();

  // Image::InterpolationType (even though this is part of Image, it feels cleaner to keep it global in the module)
  py::enum_<Image::InterpolationType>(m, "InterpolationType")
  .value("Linear", Image::InterpolationType::Linear)
  .value("NearestNeighbor", Image::InterpolationType::NearestNeighbor)
  .export_values();

  // Image
  py::class_<Image>(m, "Image")
  .def(py::init<const std::string &>())
  .def(py::init<const Image&>())

  // Image constructor from numpy array (copies array, ensuring )
  .def(py::init
       ([](py::array& np_array) {
          // get input array info
          auto info = np_array.request();

#if 0
          /*
          struct buffer_info {
            void *ptr;
            py::ssize_t itemsize;
            std::string format;
            py::ssize_t ndim;
            std::vector<py::ssize_t> shape;
            std::vector<py::ssize_t> strides;
          };
          */

          std::cout << "buffer info: \n"
                    << "\tinfo.itemsize: " << info.itemsize << std::endl
                    << "\tinfo.format: " << info.format << std::endl
                    << "\tinfo.ndim: " << info.ndim << std::endl;
          std::cout << "\tinfo.shape: [ ";
          for (int i = 0; i < info.ndim; i++) {
            std::cout << info.shape[i] << " ";
          }
          std::cout << "]\n\tinfo.strides: [ ";
          for (int i = 0; i < info.ndim; i++) {
            std::cout << info.strides[i] << " ";
          }
          std::cout << "]\n";
          std::cout << "writeable: " << np_array.writeable() << std::endl
                    << "owns data: " << np_array.owndata() << std::endl;
#endif

          // Verify info type is same as Image::PixelType (currently hard-coded as float). The
          // reasons we don't simply specify py::array_t<float> as a parameter are:
          // - to show an error if another type is sent, and
          // - to ensure the array isn't silently cast, the default of pybind11.
          if (info.format != py::format_descriptor<Image::PixelType>::format()) {
            throw std::invalid_argument("array must be of dtype.float32");
          }

          // verify it's 2d or 3d
          if (info.ndim < 2 || info.ndim > 3) {
            throw std::invalid_argument(std::string("array must be 2d or 3d, but ndim = ") + std::to_string(info.ndim));
          }

          // verify data is densely packed by checking strides is same as shape
          const py::ssize_t scalar_size = sizeof(Image::PixelType);
          std::vector<py::ssize_t> strides{info.shape[2]*info.shape[1]*scalar_size,
                                           info.shape[2]*scalar_size,
                                           scalar_size};
          for (int i = 0; i < info.ndim; i++) {
            if (info.strides[i] != strides[i]) {
              std::cerr << "expected: " << strides[i] << ", actual: " << info.strides[i] << std::endl;
              throw std::invalid_argument("array must be densely packed");
            }
          }

          // if (dtype.float32) just steal array from python (#966):
          // - np_array owners data = false
          // - itk importer pass ownership = true

          // create itk importer to copy data into image
          using ImportType = itk::ImportImageFilter<Image::PixelType, 3>;
          auto importer = ImportType::New();
          importer->SetImportPointer(static_cast<Image::PixelType *>(info.ptr), np_array.size(), false /*pass ownership*/);

          ImportType::SizeType size;            // i.e., Dims
          size[0] = np_array.shape()[2];
          size[1] = np_array.shape()[1];
          size[2] = np_array.shape()[0];

          ImportType::IndexType start({0,0,0}); // i.e., Coord
          ImportType::RegionType region;
          region.SetIndex(start);
          region.SetSize(size);

          importer->SetRegion(region);
          importer->Update();
          return Image(importer->GetOutput());
        }),
       "initialize an image from a numpy array of dtype.float32")

  .def("__neg__", [](Image& img) -> decltype(auto) { return -img; })
  .def(py::self + py::self)
  .def(py::self += py::self)
  .def(py::self - py::self)
  .def(py::self -= py::self)
  .def(py::self * Image::PixelType())
  .def(py::self *= Image::PixelType())
  .def(py::self / Image::PixelType())
  .def(py::self /= Image::PixelType())
  .def(py::self + Image::PixelType())
  .def(py::self += Image::PixelType())
  .def(py::self - Image::PixelType())
  .def(py::self -= Image::PixelType())
  .def(py::self == py::self)

  .def("__repr__",
       [](const Image &img) -> decltype(auto) {
         std::stringstream stream;
         stream << img;
         return stream.str();
       })

  .def("copy",
       [](Image& image) -> decltype(auto) { return Image(image); })

  .def("write",
       &Image::write,
       "writes the current image (determines type by its extension)",
       "filename"_a, "compressed"_a=true)

  .def("antialias",
       &Image::antialias,
       "antialiases binary volumes (layers is set to 3 when not specified)",
       "iterations"_a=50, "maxRMSErr"_a=0.01f, "layers"_a=3)

  .def("resample",
       [](Image& image,
          const TransformPtr transform,
          const std::vector<double>& p,
          const std::vector<unsigned>& d,
          const std::vector<double>& v,
          const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &direction,
          Image::InterpolationType interp) -> decltype(auto) {
         return image.resample(transform,
                               Point({p[0], p[1], p[2]}),
                               Dims({d[0], d[1], d[2]}),
                               makeVector({v[0], v[1], v[2]}),
                               eigenToItk(direction),
                               interp);
       },
       "resamples by applying transform then sampling from given origin along direction axes at spacing physical units per pixel for dims pixels using specified interpolator",
       "transform"_a, "origin"_a, "dims"_a, "spacing"_a, "direction"_a,
       "interp"_a = Image::InterpolationType::NearestNeighbor)

  .def("resample",
       [](Image& image,
          const std::vector<double>& v,
          Image::InterpolationType interp) -> decltype(auto) {
         return image.resample(makeVector({v[0], v[1], v[2]}), interp);
       },
       "resamples image using new physical spacing, updating logical dims to keep all image data for this spacing",
       "physicalSpacing"_a, "interp"_a=Image::InterpolationType::Linear)

  .def("resample",
       py::overload_cast<double, Image::InterpolationType>(&Image::resample),
       "isotropically resamples image using giving isospacing",
       "isoSpacing"_a=1.0, "interp"_a=Image::InterpolationType::Linear)

  .def("resize",
       [](Image& image, std::vector<unsigned>& d, Image::InterpolationType interp) -> decltype(auto) {
         return image.resize(Dims({d[0], d[1], d[2]}), interp);
       },
       "change logical dims (computes new physical spacing)",
       "logicalDims"_a, "interp"_a=Image::InterpolationType::Linear)

  .def("recenter",
       &Image::recenter,
       "recenters an image by changing its origin in the image header to the physical coordinates of the center of the image")

  .def("pad",
       py::overload_cast<int, Image::PixelType>(&Image::pad),
       "pads an image by same number of pixels in all directions with constant value",
       "pad"_a, "value"_a=0.0)

  .def("pad",
       py::overload_cast<int, int, int, Image::PixelType>(&Image::pad),
       "pads an image by desired number of pixels in each direction with constant value",
       "padx"_a, "pady"_a, "padz"_a, "value"_a=0.0)

  .def("pad",
       py::overload_cast<IndexRegion&, Image::PixelType>(&Image::pad),
       "pads an image to include the given region with constant value",
       "region"_a, "value"_a=0.0)

  .def("translate",
       [](Image& image, const std::vector<double>& v) -> decltype(auto) {
         return image.translate(makeVector({v[0], v[1], v[2]}));
       },
       "translates image", "v"_a)

  .def("scale",
       [](Image& image, const std::vector<double>& scale_vec) -> decltype(auto) {
         return image.scale(makeVector({scale_vec[0], scale_vec[1], scale_vec[2]}));
       },
       "scale image by scale_vec around center (not origin)",
       "scale_vec"_a)

  .def("rotate",
       py::overload_cast<const double, const Vector3&>(&Image::rotate),
       "rotate around center (not origin) using axis (default z-axis) by angle (in radians)",
       "angle"_a, "axis"_a)

  .def("rotate",
       py::overload_cast<const double, Axis>(&Image::rotate),
       "rotate around center (not origin) using axis (default z-axis) by angle (in radians)",
       "angle"_a, "axis"_a)

  .def("rotate",
       [](Image& image, const double angle, const std::vector<double>& v) -> decltype(auto) {
         return image.rotate(angle, makeVector({v[0], v[1], v[2]}));
       },
       "rotate around center (not origin) using axis (default z-axis) by angle (in radians)",
       "angle"_a, "axis"_a)

  .def("applyTransform",
       py::overload_cast<TransformPtr, Image::InterpolationType>(&Image::applyTransform),
       "applies the given transformation to the image by using resampling filter",
       "transform"_a, "interp"_a=Image::InterpolationType::Linear)

  .def("applyTransform",
       [](Image& image, const TransformPtr transform,
          const std::vector<double>& p,
          const std::vector<unsigned>& d,
          const std::vector<double>& v,
          const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> &direction,
          Image::InterpolationType interp) {
         return image.applyTransform(transform,
                                     Point({p[0], p[1], p[2]}),
                                     Dims({d[0], d[1], d[2]}),
                                     makeVector({v[0], v[1], v[2]}),
                                     eigenToItk(direction),
                                     interp);
       },
       "applies the given transformation to the image by using resampling filter with new origin, dims, spacing and direction values",
       "transform"_a, "origin"_a, "dims"_a, "spacing"_a, "direction"_a,
       "interp"_a=Image::InterpolationType::NearestNeighbor)

  .def("extractLabel",
       &Image::extractLabel,
       "extracts/isolates a specific pixel label from a given multi-label volume and outputs the corresponding binary image",
       "label"_a=1.0)

  .def("closeHoles",
       &Image::closeHoles,
       "closes holes in a volume defined by values larger than specified value",
       "foreground"_a=0.0)

  .def("binarize",
       &Image::binarize,
       "sets portion of image greater than min and less than or equal to max to the specified value",
       "minVal"_a=0.0, "maxVal"_a=std::numeric_limits<Image::PixelType>::max(),
       "innerVal"_a=1.0, "outerVal"_a=0.0)

  .def("computeDT",
       &Image::computeDT,
       "computes signed distance transform volume from an image at the specified isovalue",
       "isovalue"_a=0.0)

  .def("applyCurvatureFilter",
       &Image::applyCurvatureFilter,
       "denoises an image using curvature driven flow using curvature flow image filter",
       "iterations"_a=10)

  .def("applyGradientFilter",
       &Image::applyGradientFilter,
       "computes gradient magnitude at each pixel using gradient magnitude filter")

  .def("applySigmoidFilter",
       &Image::applySigmoidFilter,
       "computes sigmoid function pixel-wise using sigmoid image filter",
       "alpha"_a=10.0, "beta"_a=10.0)

  .def("applyTPLevelSetFilter",
       &Image::applyTPLevelSetFilter,
       "segments structures in image using topology preserving geodesic active contour level set filter",
       "featureImage"_a, "scaling"_a=20.0)

  .def("applyIntensityFilter",
       &Image::applyIntensityFilter,
       "applies intensity windowing image filter",
       "min"_a=0.0, "max"_a=0.0)

  .def("gaussianBlur",
       &Image::gaussianBlur,
       "applies gaussian blur",
       "sigma"_a=0.0)

  .def("crop",
       &Image::crop,
       "crops the image down to the given (physica) region, with optional padding",
       "region"_a, "padding"_a=0)

  .def("clip",
       [](Image& image,
          const std::vector<double>& o,
          std::vector<double>& p1,
          std::vector<double>& p2,
          const Image::PixelType val) {
         return image.clip(makePlane(Point({o[0], o[1], o[2]}),
                                     Point({p1[0], p1[1], p1[2]}),
                                     Point({p2[0], p2[1], p2[2]})), val);
       },
       "sets values on the back side of cutting plane (containing three non-colinear points) to val (default 0.0)",
       "o"_a, "p1"_a, "p2"_a, "val"_a=0.0)

  .def("clip",
       [](Image& image,
          const std::vector<double>& n,
          std::vector<double>& q,
          const Image::PixelType val) {
         return image.clip(makePlane(Point({q[0], q[1], q[2]}), makeVector({n[0], n[1], n[2]})), val);
       },
       "sets values on the back side of cutting plane (normal n containing point p) to val (default 0.0)",
       "n"_a, "q"_a, "val"_a=0.0)

  .def("setSpacing",
       [](Image& self, std::vector<double>& v) -> decltype(auto) {
         return self.setSpacing(makeVector({v[0], v[1], v[2]}));
       },
       "set image spacing, the size of each pixel",
       "spacing"_a=std::vector<double>({1.0, 1.0, 1.0}))

  .def("setOrigin",
       [](Image& image, std::vector<double>& p) {
         return image.setOrigin(Point({p[0], p[1], p[2]}));
       },
       "sets the image origin in physical space to the given value",
       "origin"_a=std::vector<double>({0,0,0}))

  .def("reflect",
       &Image::reflect,
       "reflect image with respect to logical image center and the specified axis",
       "axis"_a)

  .def("dims",
       [](const Image& self) -> decltype(auto) { return py::array(3, self.dims().data()); },
       "logical dimensions of the image")

  .def("size",
       [](const Image& self) -> decltype(auto) { return py::array(3, self.size().GetDataPointer()); },
       "physical dimensions of the image (dims * spacing)")

  .def("spacing",
       [](const Image& self) -> decltype(auto) { return py::array(3, self.spacing().GetDataPointer()); },
       "physical spacing of the image")

  .def("origin",
       [](const Image& self) -> decltype(auto) { return py::array(3, self.origin().GetDataPointer()); },
       "physical coordinates of image origin")

  .def("center",
       [](const Image& self) -> decltype(auto) { return py::array(3, self.center().GetDataPointer()); },
       "physical coordinates of center of this image")

  .def("coordsys",
       [](const Image &self) -> decltype(auto) { return itkToEigen(self.coordsys()); },
       "return coordinate system in which this image lives in physical space")

  .def("centerOfMass",
       [](const Image& self, double minVal, double maxVal) -> decltype(auto) {
         return py::array(3, self.centerOfMass().GetDataPointer());
       },
       "returns average physical coordinate of pixels in range (minval, maxval]",
       "minVal"_a=0.0, "maxVal"_a=1.0)

  .def("min",
       &Image::min,
       "minimum of image")

  .def("max",
       &Image::max,
       "maximum of image")

  .def("mean",
       &Image::mean,
       "mean of image")

  .def("std",
       &Image::std,
       "standard deviation of image")

  .def("logicalBoundingBox",
       &Image::logicalBoundingBox,
       "returns the index coordinates of this image's region")

  // py::overload_cast doesn't work with const functions, so we directly static_cast for these two functions
  .def("physicalboundingBox",
       static_cast<PhysicalRegion (Image::*)() const>(&Image::physicalBoundingBox),
       "returns region of physical space occupied by this image")

  .def("physicalBoundingBox",
       static_cast<PhysicalRegion (Image::*)(Image::PixelType) const>(&Image::physicalBoundingBox),
       "returns region of physical space occupied by the region of data <= the given isoValue",
       "isovalue"_a=1.0)

  .def("logicalToPhysical",
       [](Image& self, IndexRegion region) -> decltype(auto) {
         return self.logicalToPhysical(region);
       },
       "converts from a logical region (index coordinates) to a physical region",
       "region"_a)

  .def("logicalToPhysical",
       [](Image& self, std::vector<long>& c) -> decltype(auto) {
         return py::array(3, self.logicalToPhysical(Coord({c[0], c[1], c[2]})).GetDataPointer());
       },
       "converts a logical (index) coordinate to physical space",
       "c"_a)

  .def("physicalToLogical",
       [](Image& self, PhysicalRegion region) -> decltype(auto) {
         return self.physicalToLogical(region);
       },
       "converts from a physical region to a logical region (index coordinates)",
       "region"_a)

  .def("physicalToLogical",
       [](Image& self, std::vector<double>& p) -> decltype(auto) {
         return py::array(3, self.physicalToLogical(Point({p[0], p[1], p[2]})).data());
       },
       "converts a physical coordinate to a logical (index) space",
       "p"_a)

  .def("compare",
       &Image::compare,
       "compares two images",
       "other"_a, "verifyall"_a=true, "tolerance"_a=0.0, "precision"_a=1e-12)

  .def("toArray",
       [](const Image &image) -> decltype(auto) {
         const auto dims = image.dims();
         const auto shape = std::vector<size_t>{dims[2], dims[1], dims[0]};
         return py::array(py::dtype::of<typename Image::ImageType::Pointer::ObjectType::PixelType>(),
                          shape, image.getITKImage()->GetBufferPointer());
       },
       "returns raw array of image data (note: spacing, origin, coordsys are not preserved)")

  .def("createTransform",
       py::overload_cast<XFormType>(&Image::createTransform),
       "creates a transform based on transform type",
       "type"_a=XFormType::CenterOfMass)

  .def("createTransform",
       py::overload_cast<const Image&, XFormType, float, unsigned>(&Image::createTransform),
       "creates a transform based on transform type",
       "target"_a, "type"_a=XFormType::IterativeClosestPoint, "isoValue"_a=0.0, "iterations"_a=20)

  .def("topologyPreservingSmooth",
       &Image::topologyPreservingSmooth,
       "creates a feature image (by applying gradient then sigmoid filters), then passes it to the TPLevelSet filter [curvature flow filter is often applied to the image before this filter]",
       "scaling"_a=20.0, "sigmoidAlpha"_a=10.5, "sigmoidBeta"_a=10.0)

  .def("compare",
       &Image::compare,
       "compares two images",
       "other"_a, "verifyall"_a=true, "tolerance"_a=0.0, "precision"_a=1e-12)

  .def("toMesh",
       [](Image &image, Image::PixelType isovalue) {
         return image.toMesh(isovalue);
       },
       "converts image to mesh at specified isovalue",
       "isovalue"_a)
  ;

  // PhysicalRegion
  py::class_<PhysicalRegion>(m, "PhysicalRegion")

  .def(py::init<>())

  .def(py::init
       ([](std::vector<double> min, std::vector<double> max) {
          return PhysicalRegion(Point({min[0], min[1], min[2]}),
                                Point({max[0], max[1], max[2]}));
        }))

  .def(py::self == py::self)

  .def("__repr__",
       [](const PhysicalRegion &region) {
         std::stringstream stream;
         stream << region;
         return stream.str();
       })

  .def_property("min",
                // read (return special array so members of returned point can be modified, e.g., min[0] = 1.0)
                py::cpp_function([](PhysicalRegion &region) -> py::array_t<double> {
                    py::str dummyDataOwner; // pretend the data has an owner and it won't be copied (pybind trick)
                    py::array arr(py::dtype::of<Point::ValueType>(),      // dtype
                                  std::vector<ssize_t>({3}),              // shape
                                  std::vector<ssize_t>(),                 // spacing
                                  region.min.GetDataPointer(),            // data ptr
                                  dummyDataOwner);     // "inspire" py::array not to copy data
                    return arr;
                  }, py::return_value_policy::move),

                // assign
                py::cpp_function([](PhysicalRegion &region, std::vector<double> min) -> decltype(auto) {
                                   region.min = Point({min[0], min[1], min[2]});
                                   return min;
                                 }, py::return_value_policy::reference),
                "min point of region")

  .def_property("max",
                // read (return special array so members of returned point can be modified, e.g., max[0] = 1.0)
                py::cpp_function([](PhysicalRegion &region) -> py::array_t<double> {
                    py::str dummyDataOwner; // pretend the data has an owner and it won't be copied (pybind trick)
                    py::array arr(py::dtype::of<Point::ValueType>(),      // dtype
                                  std::vector<ssize_t>({3}),              // shape
                                  std::vector<ssize_t>(),                 // spacing
                                  region.max.GetDataPointer(),            // data ptr
                                  dummyDataOwner);     // "inspire" py::array not to copy data
                    return arr;
                  }, py::return_value_policy::move),

                // assign
                py::cpp_function([](PhysicalRegion &region, std::vector<double> max) -> decltype(auto) {
                                   region.max = Point({max[0], max[1], max[2]});
                                   return max;
                                 }, py::return_value_policy::reference),
                "max point of region")

  .def("valid",
       &PhysicalRegion::valid,
       "ensure if region is valid")

  .def("origin",
       [](const PhysicalRegion &region) -> decltype(auto) {
         return py::array(3, region.origin().GetDataPointer());
       },
       "return origin of region")

  .def("size",
       [](const PhysicalRegion &region) -> decltype(auto) {
         return py::array(3, region.size().GetDataPointer());
       },
       "return size of region")

  .def("shrink",
       &PhysicalRegion::shrink,
       "shrink this region down to the smallest portions of both",
       "other"_a)

  .def("expand",
       py::overload_cast<const PhysicalRegion&>(&PhysicalRegion::expand),
       "expand this region up to the largest portions of both",
       "other"_a)

  .def("expand",
       py::overload_cast<const Point&>(&PhysicalRegion::expand),
       "expand this region to include this point",
       "point"_a)

  .def("pad",
       &PhysicalRegion::pad,
       "grows or shrinks the region by the specified amount",
       "padding"_a)
  ;

  // IndexRegion
  py::class_<IndexRegion>(m, "IndexRegion")

  .def(py::init<>())

  .def(py::init
       ([](std::vector<double> min, std::vector<double> max) {
          return IndexRegion(Coord({static_cast<long>(min[0]),
                                    static_cast<long>(min[1]),
                                    static_cast<long>(min[2])}),
                             Coord({static_cast<long>(max[0]),
                                    static_cast<long>(max[1]),
                                    static_cast<long>(max[2])}));
        }))

  .def(py::init
       ([](std::vector<double> dims) {
          if (dims[0] < 0.0 || dims[1] < 0.0 || dims[2] < 0.0)
            throw std::invalid_argument("cannot convert negative values to IndexRegion");
          return IndexRegion(Dims({static_cast<unsigned long>(dims[0]),
                                   static_cast<unsigned long>(dims[1]),
                                   static_cast<unsigned long>(dims[2])}));
        }))

  .def(py::self == py::self)

  .def("__repr__",
       [](const IndexRegion &region) {
         std::stringstream stream;
         stream << region;
         return stream.str();
       })

  .def_property("min",
                // read (return special array so members of returned point can be modified, e.g., min[0] = 1.0)
                py::cpp_function([](IndexRegion &region) -> py::array_t<Coord::IndexValueType> {
                    py::str dummyDataOwner; // pretend the data has an owner and it won't be copied (pybind trick)
                    py::array arr(py::dtype::of<Coord::IndexValueType>(), // dtype
                                  std::vector<ssize_t>({3}),              // shape
                                  std::vector<ssize_t>(),                 // spacing
                                  region.min.data(),                      // data ptr
                                  dummyDataOwner);     // "inspire" py::array not to copy data
                    return arr;
                  }, py::return_value_policy::move),

                // assign
                py::cpp_function([](IndexRegion &region, std::vector<double> min) -> decltype(auto) {
                                   region.min = Coord({static_cast<Coord::IndexValueType>(min[0]),
                                                       static_cast<Coord::IndexValueType>(min[1]),
                                                       static_cast<Coord::IndexValueType>(min[2])});
                                   return min;
                                 }, py::return_value_policy::reference),
                "min point of region")

  .def_property("max",
                // read (return special array so members of returned point can be modified, e.g., max[0] = 1.0)
                py::cpp_function([](IndexRegion &region) -> py::array_t<Coord::IndexValueType> {
                    py::str dummyDataOwner; // pretend the data has an owner and it won't be copied (pybind trick)
                    py::array arr(py::dtype::of<Coord::IndexValueType>(), // dtype
                                  std::vector<ssize_t>({3}),              // shape
                                  std::vector<ssize_t>(),                 // spacing
                                  region.max.data(),                      // data ptr
                                  dummyDataOwner);     // "inspire" py::array not to copy data
                    return arr;
                  }, py::return_value_policy::move),

                // assign
                py::cpp_function([](IndexRegion &region, std::vector<double> max) -> decltype(auto) {
                                   region.max = Coord({static_cast<Coord::IndexValueType>(max[0]),
                                                       static_cast<Coord::IndexValueType>(max[1]),
                                                       static_cast<Coord::IndexValueType>(max[2])});
                                   return max;
                                 }, py::return_value_policy::reference),
                "max point of region")

  .def("valid",
       &IndexRegion::valid,
       "ensure if region is valid")

  .def("origin",
       [](const IndexRegion &region) -> decltype(auto) {
         return py::array(3, region.origin().data());
       },
       "return origin of region")

  .def("size",
       [](const IndexRegion &region) -> decltype(auto) {
         return py::array(3, region.size().data());
       },
       "return size of region")

  .def("pad",
       &IndexRegion::pad,
       "grows or shrinks the region by the specified amount",
       "padding"_a)
  ;

  // VectorImage
  py::class_<VectorImage>(m, "VectorImage")

  .def(py::init
       ([](const Image &dt) {
          return VectorImage(dt);
        }),
       "create a vector image from an image (usually a distance transform) that can be sampled at any point in space",
       "image"_a)

  .def("evaluate",
       [](VectorImage &image, std::vector<double> &pt) -> decltype(auto) {
         auto v = image.evaluate(Point({pt[0], pt[1], pt[2]}));
         return std::vector<double>({v[0], v[1], v[2]});
       },
       "evaluate the vector image at any given point in space",
       "pt"_a)
  ;

  // ImageUtils
  py::class_<ImageUtils>(m, "ImageUtils")

  .def_static("boundingBox",
              py::overload_cast<const std::vector<std::string>&, Image::PixelType>(&ImageUtils::boundingBox),
              "compute largest bounding box surrounding the specified isovalue of the specified set of filenames",
              "filenames"_a, "isoValue"_a=1.0)

  .def_static("boundingBox",
              py::overload_cast<const std::vector<std::reference_wrapper<const Image>>&, Image::PixelType>(&ImageUtils::boundingBox),
              "compute largest bounding box surrounding the specified isovalue of the specified set of images",
              "images"_a, "isoValue"_a=1.0)

  .def_static("createWarpTransform",
              [](const std::string &source_landmarks,
                 const std::string &target_landmarks,
                 const int stride) -> decltype(auto) {
                auto xform_ptr = ImageUtils::createWarpTransform(source_landmarks,
                                                                 target_landmarks,
                                                                 stride);
                return xform_ptr;
              },
              "computes a warp transform from the source to the target landmarks",
              "source_landmarks"_a, "target_landmarks"_a, "stride"_a=1)
  ;

  // Mesh
  py::class_<Mesh> mesh(m, "Mesh");

  // Mesh::AlignmentType
  py::enum_<Mesh::AlignmentType>(mesh, "AlignmentType")
  .value("Rigid", Mesh::AlignmentType::Rigid)
  .value("Similarity", Mesh::AlignmentType::Similarity)
  .value("Affine", Mesh::AlignmentType::Affine)
  .export_values();
  ;

  // Mesh::DistanceMethod
  py::enum_<Mesh::DistanceMethod>(mesh, "DistanceMethod")
  .value("POINT_TO_POINT", Mesh::DistanceMethod::POINT_TO_POINT)
  .value("POINT_TO_CELL", Mesh::DistanceMethod::POINT_TO_CELL)
  .export_values();
  ;

  // Mesh bindings
  mesh.def(py::init<const std::string &>())

  .def(py::init<vtkSmartPointer<vtkPolyData>>())
  .def(py::self == py::self)
  .def(py::self += py::self)

  .def("__repr__",
       [](const Mesh &mesh) -> decltype(auto) {
         std::stringstream stream;
         stream << mesh;
         return stream.str();
       })

  .def("copy",
       [](Mesh& mesh) -> decltype(auto) { return Mesh(mesh); })

  .def("write",
       &Mesh::write,
       "writes mesh, format specified by filename extension",
       "pathname"_a)

  .def("coverage",
       &Mesh::coverage,
       "determines coverage between current mesh and another mesh (e.g. acetabular cup / femoral head)",
       "otherMesh"_a, "allowBackIntersections"_a=true, "angleThreshold"_a=0, "backSearchRadius"_a=0)

  .def("smooth",
       &Mesh::smooth,
       "applies laplacian smoothing",
       "iterations"_a=0, "relaxation"_a=0.0)

  .def("smoothSinc",
       &Mesh::smoothSinc,
       "applies windowed sinc smoothing",
       "iterations"_a=0, "passband"_a=0.0)

  .def("decimate",
       &Mesh::decimate,
       "applies filter to reduce number of triangles in mesh",
       "reduction"_a=0.5, "angle"_a=15.0, "preserveTopology"_a=true)

   .def("cvdDecimate",
       &Mesh::cvdDecimate,
       "applies cvd (centroidal voronoi diagram) decimation filter",
       "percentage"_a=0.5)

  .def("invertNormals",
       &Mesh::invertNormals,
       "handle flipping normals")

  .def("reflect",
       [](Mesh& mesh, const Axis &axis, std::vector<double>& v) -> decltype(auto) {
         return mesh.reflect(axis, makeVector({v[0], v[1], v[2]}));
       },
       "reflect meshes with respect to a specified center and specific axis",
       "axis"_a, "origin"_a=std::vector<double>({0.0, 0.0, 0.0}))

  .def("createTransform",
       &Mesh::createTransform,
       "creates a transform based on transform type",
       "target"_a, "type"_a=XFormType::IterativeClosestPoint,
       "align"_a=Mesh::AlignmentType::Similarity, "iterations"_a=10)

  .def("applyTransform",
       &Mesh::applyTransform,
       "applies the given transformation to the mesh",
       "transform"_a)

  .def("fillHoles",
       &Mesh::fillHoles,
       "finds holes in a mesh and closes them")

  .def("probeVolume",
       &Mesh::probeVolume,
       "samples data values at specified point locations",
       "image"_a)

  .def("clip",
       [](Mesh& mesh, const std::vector<double>& p, const std::vector<double>& n) -> decltype(auto) {
         return mesh.clip(makePlane(Point({p[0], p[1], p[2]}), makeVector({n[0], n[1], n[2]})));
       },
       "clips a mesh using a cutting plane",
       "point"_a,
       "normal"_a)

  .def("clip",
       [](Mesh& mesh, const std::vector<double>& o, const std::vector<double>& p1, const std::vector<double>& p2) -> decltype(auto) {
         return mesh.clip(makePlane(Point({o[0], o[1], o[2]}), Point({p1[0], p1[1], p1[2]}), Point({p2[0], p2[1], p2[2]})));
       },
       "clips a mesh using a cutting plane",
       "o"_a,
       "p1"_a,
       "p2"_a)

  .def("translate",
       [](Mesh& mesh, const std::vector<double>& v) -> decltype(auto) {
         return mesh.translate(makeVector({v[0], v[1], v[2]}));
       },
       "translates mesh",
       "v"_a)

  .def("scale",
       [](Mesh& mesh, const std::vector<double>& v) -> decltype(auto) {
         return mesh.scale(makeVector({v[0], v[1], v[2]}));
       },
       "scale mesh",
       "v"_a)

  .def("boundingBox",
       &Mesh::boundingBox,
       "computes bounding box of current mesh")

  .def("fixElement",
       &Mesh::fixElement,
       "fix element winding of mesh")

  .def("clipClosedSurface",
       [](Mesh& mesh, const std::vector<double>& p, const std::vector<double>& n) -> decltype(auto) {
         return mesh.clipClosedSurface(makePlane(Point({p[0], p[1], p[2]}), makeVector({n[0], n[1], n[2]})));
       },
       "clips a mesh using a cutting plane resulting in a closed surface",
       "point"_a,
       "normal"_a)

  .def("computeNormals",
       &Mesh::computeNormals,
       "computes and adds oriented point and cell normals")

  .def("toImage",
       [](Mesh& mesh, PhysicalRegion &region, std::vector<double>& spacing) -> decltype(auto) {
         return mesh.toImage(region, Point({spacing[0], spacing[1], spacing[2]}));
       },
       "rasterizes specified region to create binary image of desired dims (default: unit spacing)",
       "region"_a=PhysicalRegion(),
       "spacing"_a=std::vector<double>({1.0, 1.0, 1.0}))

  .def("closestPoint",
       [](Mesh &mesh, std::vector<double> p) -> decltype(auto) {
         return py::array(3, mesh.closestPoint(Point({p[0], p[1], p[2]})).GetDataPointer());
       },
       "returns closest point to given point on mesh",
       "point"_a)

  .def("closestPointId",
       [](Mesh &mesh, std::vector<double> p) -> decltype(auto) {
         return mesh.closestPointId(Point({p[0], p[1], p[2]}));
       },
       "returns closest point id in this mesh to the given point in space",
       "point"_a)

  .def("geodesicDistance",
       &Mesh::geodesicDistance,
       "computes geodesic distance between two vertices (specified by their indices) on mesh",
       "source"_a, "target"_a)

  .def("distance",
       &Mesh::distance,
       "computes surface to surface distance",
       "target"_a, "method"_a=Mesh::DistanceMethod::POINT_TO_POINT)

  .def("toDistanceTransform",
       [](Mesh& mesh, PhysicalRegion &region, std::vector<double>& spacing) -> decltype(auto) {
         return mesh.toDistanceTransform(region, Point({spacing[0], spacing[1], spacing[2]}));
       },
       "converts specified region to distance transform image (default: unit spacing)",
       "region"_a=PhysicalRegion(),
       "spacing"_a=std::vector<double>({1.0, 1.0, 1.0}))

  .def("center",
       [](Mesh &mesh) -> decltype(auto) {
         return py::array(3, mesh.center().GetDataPointer());
       },
       "center of mesh")

  .def("centerOfMass",
       [](Mesh &mesh) -> decltype(auto) {
         return py::array(3, mesh.centerOfMass().GetDataPointer());
       },
       "center of mass of mesh")

  .def("numPoints",
       &Mesh::numPoints,
       "number of points")

  .def("numFaces",
       &Mesh::numFaces,
       "number of faces")

  .def("points",
       &Mesh::points,
       "matrix with number of points with (x,y,z) coordinates of each point")

  .def("faces",
       &Mesh::faces,
       "matrix with number of faces with indices of the three points from which each face is composed")

  .def("getPoint",
       [](Mesh &mesh, int id) -> decltype(auto) {
         return py::array(3, mesh.getPoint(id).GetDataPointer());
       },
       "(x,y,z) coordinates of vertex at given index",
       "id"_a)

  .def("getFace",
       [](Mesh &mesh, int id) -> decltype(auto) {
         return py::array(3, mesh.getFace(id).GetDataPointer());
       },
       "return indices of the three points with which the face at the given index is composed",
       "id"_a)

  .def("getFieldNames",
       &Mesh::getFieldNames,
       "print all field names in mesh")

  //TODO: See github issue #966
  .def("setField",
       [](Mesh &mesh, std::vector<double>& v, std::string name) -> decltype(auto) {
         vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
         arr->SetNumberOfValues(v.size());
         for (int i=0; i<v.size(); i++) {
           arr->SetTuple1(i, v[i]);
         }
         return mesh.setField(name, arr);
       },
       "sets the given field for points with array",
       "array"_a, "name"_a)

  //TODO: See github issue #966
  .def("setField",
       [](Mesh &mesh, std::vector<std::vector<double>>& v, std::string name) -> decltype(auto) {
         vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
         arr->SetNumberOfComponents(3);
         arr->SetNumberOfTuples(v.size());
         for (int i=0; i<v.size(); i++) {
           arr->SetTuple3(i, v[i][0], v[i][1], v[i][2]);
          }
         return mesh.setField(name, arr);
       },
       "sets the given field for points with array",
       "array"_a, "name"_a)

  .def("getField",
       [](const Mesh &mesh, std::string name) -> decltype(auto) {
         auto array = mesh.getField<vtkDataArray>(name);
         const auto shape = std::vector<size_t>{static_cast<unsigned long>(array->GetNumberOfTuples()),
                                                static_cast<unsigned long>(array->GetNumberOfComponents()),
                                                1};
         auto vtkarr = vtkSmartPointer<vtkDoubleArray>(vtkDoubleArray::New());
         vtkarr->SetNumberOfValues(array->GetNumberOfValues());

         // LOTS of copying going on here, see github #903
         array->GetData(0, array->GetNumberOfTuples()-1,
                        0, array->GetNumberOfComponents()-1,
                        vtkarr);                               // copy1
         return py::array(py::dtype::of<double>(),
                          shape,
                          vtkarr->GetVoidPointer(0));          // copy2
       },
       "gets the field",
       "name"_a)

  .def("setFieldValue",
       &Mesh::setFieldValue,
       "sets the given index of field to value",
       "idx"_a, "value"_a, "name"_a="")

  .def("getFieldValue",
       &Mesh::getFieldValue,
       "gets the value at the given index of field",
       "idx"_a, "name"_a)

  .def("getFieldRange",
       &Mesh::getFieldRange,
       "returns the range of the given field",
       "name"_a)

  .def("getFieldMean",
       &Mesh::getFieldMean,
       "returns the mean the given field",
       "name"_a)

  .def("getFieldStd",
       &Mesh::getFieldStd,
       "returns the standard deviation of the given field",
       "name"_a)

  .def("compareField",
       &Mesh::compareField,
       "compares two meshes based on fields",
       "other_mesh"_a, "name1"_a, "name2"_a="")
  ;

  // MeshUtils
  py::class_<MeshUtils>(m, "MeshUtils")

  .def_static("boundingBox",
              py::overload_cast<const std::vector<std::string>&, bool>(&MeshUtils::boundingBox),
              "calculate bounding box incrementally for meshes",
              "filenames"_a, "center"_a=false)

  .def_static("boundingBox",
              py::overload_cast<const std::vector<std::reference_wrapper<const Mesh>>&, bool>(&MeshUtils::boundingBox),
              "calculate bounding box incrementally for meshes",
              "meshes"_a, "center"_a=false)

  .def_static("findReferenceMesh",
              &MeshUtils::findReferenceMesh,
              "find reference mesh from a set of meshes",
              "meshes"_a)

  .def_static("generateNormals",
              &MeshUtils::generateNormals,
              "generates and adds normals for points and faces for each mesh in given set of meshes",
              "meshes"_a, "forceRegen"_a=false)

  .def_static("computeMeanNormals",
               [](const std::vector<std::string>& filenames, bool autoGenerateNormals) -> decltype(auto) {
                  auto array = MeshUtils::computeMeanNormals(filenames, autoGenerateNormals);
                  const auto shape = std::vector<size_t>{static_cast<unsigned long>(array->GetNumberOfTuples()),
                                                         static_cast<unsigned long>(array->GetNumberOfComponents()),
                                                         1};
                  auto vtkarr = vtkSmartPointer<vtkDoubleArray>(vtkDoubleArray::New());
                  vtkarr->SetNumberOfValues(array->GetNumberOfValues());

                  // LOTS of copying going on here, see github #903
                  array->GetData(0, array->GetNumberOfTuples()-1,
                                 0, array->GetNumberOfComponents()-1,
                                 vtkarr);                               // copy1
                  return py::array(py::dtype::of<double>(),
                                   shape,
                                   vtkarr->GetVoidPointer(0));          // copy2
               },
               "computes average normals for each point in given set of meshes",
               "filenames"_a, "autoGenerateNormals"_a=true)

  .def_static("computeMeanNormals",
               [](const std::vector<std::reference_wrapper<const Mesh>>& meshes) -> decltype(auto) {
                  auto array = MeshUtils::computeMeanNormals(meshes);
                  const auto shape = std::vector<size_t>{static_cast<unsigned long>(array->GetNumberOfTuples()),
                                                         static_cast<unsigned long>(array->GetNumberOfComponents()),
                                                         1};
                  auto vtkarr = vtkSmartPointer<vtkDoubleArray>(vtkDoubleArray::New());
                  vtkarr->SetNumberOfValues(array->GetNumberOfValues());

                  // LOTS of copying going on here, see github #903
                  array->GetData(0, array->GetNumberOfTuples()-1,
                                 0, array->GetNumberOfComponents()-1,
                                 vtkarr);                               // copy1
                  return py::array(py::dtype::of<double>(),
                                   shape,
                                   vtkarr->GetVoidPointer(0));          // copy2
               },
               "computes average normals for each point in given set of meshes",
               "meshes"_a)
  ;

  // ParticleSystem
  py::class_<ParticleSystem>(m, "ParticleSystem")

  .def(py::init<const std::vector<std::string> &>())

  .def("Particles",
       &ParticleSystem::Particles)

  .def("Paths",
       &ParticleSystem::Paths)

  .def("N",
       &ParticleSystem::N)

  .def("D",
       &ParticleSystem::D)
  ;

  // ShapeEvaluation
  py::class_<ShapeEvaluation>(m, "ShapeEvaluation")

  .def_static("ComputeCompactness",
              &ShapeEvaluation::ComputeCompactness,
              "particleSystem"_a, "nModes"_a, "saveTo"_a="")

  .def_static("ComputeGeneralization",
              &ShapeEvaluation::ComputeGeneralization,
              "particleSystem"_a, "nModes"_a, "saveTo"_a="")

  .def_static("ComputeSpecificity",
              &ShapeEvaluation::ComputeSpecificity,
              "particleSystem"_a, "nModes"_a, "saveTo"_a="")

  .def_static("ComputeFullCompactness",
              &ShapeEvaluation::ComputeFullCompactness,
              "particleSystem"_a,"progress_callback"_a=nullptr)

  .def_static("ComputeFullGeneralization",
              &ShapeEvaluation::ComputeFullGeneralization,
              "particleSystem"_a,"progress_callback"_a=nullptr)

  .def_static("ComputeFullSpecificity",
              &ShapeEvaluation::ComputeFullSpecificity,
              "particleSystem"_a,"progress_callback"_a=nullptr)
  ;

  py::class_<ParticleShapeStatistics>(m, "ParticleShapeStatistics")

  .def(py::init<>())

  .def("PCA",
       py::overload_cast<ParticleSystem, int>(&ParticleShapeStatistics::DoPCA),
       "calculates the eigen values and eigen vectors of the data",
       "particleSystem"_a, "domainsPerShape"_a=1)

  .def("principalComponentProjections",
       &ParticleShapeStatistics::PrincipalComponentProjections,
       "projects the original data on the calculated principal components")

  .def("sampleSize",
       &ParticleShapeStatistics::SampleSize,
       "returns the sample size of the particle system")

  .def("numDims",
       &ParticleShapeStatistics::NumberOfDimensions,
       "returns the number of features of the particle system")

  .def("eigenVectors",
       [](ParticleShapeStatistics &stats) -> decltype(auto) {
         return vnlToEigen(stats.Eigenvectors());
       })

  .def("eigenValues",
       &ParticleShapeStatistics::Eigenvalues,
       "return the eigen values")

  .def("pcaLoadings",
       &ParticleShapeStatistics::PCALoadings,
       "returns the  coefficients of the linear combination of \
       the original variables from which the principal \
       components are constructed")

  .def("percentVarByMode",
       &ParticleShapeStatistics::PercentVarByMode,
       "return the variance accounted for by the principal components")
  ;

  // Optimize (TODO)
  py::class_<Optimize>(m, "Optimize")

  .def(py::init<>())

  .def("LoadParameterFile",
       &Optimize::LoadParameterFile)

  .def("Run",
       &Optimize::Run)

  .def("SetIterationCallbackFunction",
       &Optimize::SetIterationCallbackFunction)

  .def("GetParticleSystem",
       &optimize_get_particle_system)
  ;

} // PYBIND11_MODULE(shapeworks_py)
