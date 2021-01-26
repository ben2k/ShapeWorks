#include "Testing.h"

#include "Mesh.h"
#include "MeshUtils.h"
#include "Image.h"

using namespace shapeworks;

TEST(MeshTests, readTest)
{
  try {
    Mesh mesh(std::string(TEST_DATA_DIR) + "/foo.vtk");
  } catch(...) {
    return;
  }

  ASSERT_TRUE(true);
}

TEST(MeshTests, smoothTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.smooth(10, 0.01);
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/smooth1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, smoothTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.smooth();
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/smooth2.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, decimateTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.decimate();
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/decimate1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, decimateTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.decimate(0.0, 0.0, true);
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/decimate2.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

// TODO: Karthik will add MeshTest for inverNormal

TEST(MeshTests, reflectTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.reflect(Axis::X);
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/reflect1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, reflectTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.reflect(Axis::Y);
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/reflect2.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, fillHolesTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.fillHoles();
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/fillholes.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, probeTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.probeVolume(std::string(TEST_DATA_DIR) + "/femurVtkDT.nrrd");
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/probe.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, clipTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.clip(MeshUtils::createPlane(makeVector({0.0,0.0,0.0}), Point3({0.0,0.0,0.0})));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/clip1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, translateTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.translate(makeVector({1.0,1.0,1.0}));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/translate1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, translateTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.translate(makeVector({-10.0,-10.0,-10.0}));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/translate2.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, translateTest3)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.translate(makeVector({0,0,1.0}));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/translate3.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, scaleTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.scale(makeVector({1.0,1.0,1.0}));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/scale1.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

TEST(MeshTests, scaleTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  femur.scale(makeVector({-1.0, 1.5, 1.0}));
  Mesh ground_truth(std::string(TEST_DATA_DIR) + "/scale2.vtk");

  ASSERT_TRUE(femur == ground_truth);
}

// TODO: PreviewCmd results in seg fault. WHY?
// TEST(MeshTests, fixTest)
// {
//   Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
//   femur.fix();
//   Mesh ground_truth(std::string(TEST_DATA_DIR) + "/scale1.vtk");

//   ASSERT_TRUE(femur == ground_truth);
// }

TEST(MeshTests, distanceTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  double distance = femur.hausdorffDistance(pelvis);

  ASSERT_TRUE(std::abs(distance - 32.2215) < 1e-4);
}

TEST(MeshTests, relativeDistanceABTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  double relativeDistanceAB = femur.relativeDistanceAtoB(pelvis);

  ASSERT_TRUE(std::abs(relativeDistanceAB - 32.2215) < 1e-4);
}

TEST(MeshTests, relativeDistanceBATest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  double relativeDistanceBA = femur.relativeDistanceBtoA(pelvis);

  ASSERT_TRUE(std::abs(relativeDistanceBA - 16.1937) < 1e-4);
}

TEST(MeshTests, rasterizationOriginTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Region region = femur.boundingBox();
  Point3 origin({-114, -22, 1209});

  ASSERT_TRUE(femur.rasterizationOrigin(region) == origin);
}

TEST(MeshTests, rasterizationOriginTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  std::vector <Mesh> meshes;
  meshes.push_back(femur);
  meshes.push_back(pelvis);
  Region region = MeshUtils::boundingBox(meshes);
  Point3 origin({-114, -22, 1202});

  ASSERT_TRUE(femur.rasterizationOrigin(region) == origin);
}

TEST(MeshTests, rasterizationSizeTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Region region = femur.boundingBox();
  Dims size({45, 45, 41});

  ASSERT_TRUE(femur.rasterizationSize(region) == size);
}

TEST(MeshTests, rasterizationSizeTest2)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  std::vector <Mesh> meshes;
  meshes.push_back(femur);
  meshes.push_back(pelvis);
  Region region = MeshUtils::boundingBox(meshes);
  Dims size({48, 51, 54});

  ASSERT_TRUE(femur.rasterizationSize(region) == size);
}

TEST(MeshTests, center)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.ply");
  Point3 center(femur.center());
  std::cout << "center: " << center << std::endl;

  ASSERT_TRUE(epsEqual(center, Point3({90.7541, -160.557, -673.572}), 1e-3));
}

TEST(MeshTests, centerofmass)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.ply");
  Point3 com(femur.centerOfMass());
  std::cout << "com: " << com << std::endl;

  ASSERT_TRUE(epsEqual(com, Point3({92.6201, -157.662, -666.611}), 1e-3));
}

TEST(MeshTests, toImageTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.ply");
  Region region = femur.boundingBox();
  Image image = femur.toImage(makeVector({1.0,1.0,1.0}), femur.rasterizationSize(region), femur.rasterizationOrigin(region));
  Image ground_truth(std::string(TEST_DATA_DIR) + "/femurImage.nrrd");

  ASSERT_TRUE(image == ground_truth);
}

TEST(MeshTests, toDistanceTransformTest1)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.ply");
  Image image = femur.toDistanceTransform();
  Image ground_truth(std::string(TEST_DATA_DIR) + "/femurDT.nrrd");

  ASSERT_TRUE(image == ground_truth);
}

// <ctc> add toImage and toDT tests that specify params
// also, think of a way to specify padding automatically computing origin and size

TEST(MeshTests, coverageTest)
{
  Mesh femur(std::string(TEST_DATA_DIR) + "/femur.vtk");
  Mesh pelvis(std::string(TEST_DATA_DIR) + "/pelvis.vtk");
  pelvis.coverage(femur);

  Mesh baseline(std::string(TEST_DATA_DIR) + "/coverage.vtk");
  ASSERT_TRUE(pelvis.comparePointsEqual(baseline));
  ASSERT_TRUE(pelvis.compareScalarsEqual(baseline));
  ASSERT_TRUE(true);
}
