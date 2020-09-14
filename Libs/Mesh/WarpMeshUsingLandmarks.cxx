#include "TriMesh.h"
#include "meshFIM.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileReader.h"

Eigen::MatrixXd Vref;
Eigen::MatrixXd Vtemp;
Eigen::MatrixXi Fref;
Eigen::MatrixXd Vcontrol_static;

// void convertVTKtoOBJ(std::string inputFileName){

//   // tested for VTK to obj but should work for any input format
//   const unsigned int Dimension = 3;
//   typedef float TCoordinate;

//   typedef itk::Mesh< TCoordinate, Dimension > TMesh;
//   typedef itk::MeshFileReader< TMesh > TReader;
//   typedef itk::MeshFileWriter< TMesh > TWriter;
//   TReader::Pointer reader = TReader::New();
//   reader->SetFileName(inputFileName.c_str());
//   TWriter::Pointer writer = TWriter::New();
//   writer->SetFileName( "TemplateMesh.obj" );
//   writer->SetInput( reader->GetOutput() );
//   try
//     {
//     writer->Update();
//     }
//   catch( itk::ExceptionObject & error )
//     {
//     std::cerr << "Error: " << error << std::endl;
//     }

// }

// vector<Eigen::MatrixXd> W_precomputation(Eigen::MatrixXd Vcontrol_static, Eigen::MatrixXd TV, Eigen::MatrixXi TT, Eigen::MatrixXi TF){ 

//     vector<Eigen::MatrixXd> v;
//     Eigen::MatrixXd W;
//     Eigen::VectorXi b;
//     {
//       Eigen::VectorXi J = Eigen::VectorXi::LinSpaced(TV.rows(),0,TV.rows()-1);
//       Eigen::VectorXd sqrD;
//       Eigen::MatrixXd _2;
//       std::cout<<"Finding closest points..."<<std::endl;
//       igl::point_mesh_squared_distance(Vcontrol_static,TV,J,sqrD,b,_2);
//       assert(sqrD.minCoeff() < 1e-7 && "low.V should exist in high.V");
//     }
//     // force perfect positioning, rather have popping in low-res than high-res.
//     // The correct/elaborate thing to do is express original low.V in terms of
//     // linear interpolation (or extrapolation) via elements in (high.V,high.F)
//     igl::slice(TV,b,1,Vcontrol_static);
//     // list of points --> list of singleton lists
//     std::vector<std::vector<int> > S;
//     igl::matrix_to_list(b,S);
//     std::cout<<"Computing weights for "<<b.size()<<
//       " handles at "<<TV.rows()<<" vertices..."<<std::endl;
//     // Technically k should equal 3 for smooth interpolation in 3d, but 2 is
//     // faster and looks OK
//     const int k = 2;
//     igl::biharmonic_coordinates(TV,TT,S,k,W);
//     std::cout<<"Reindexing..."<< std::endl;
//     std::cout << W.rows() << " " << W.cols() << std::endl;
//     // Throw away interior tet-vertices, keep weights and indices of boundary
//     Eigen::VectorXi I,J;
//     igl::remove_unreferenced(TV.rows(),TF,I,J);
//     for_each(TF.data(),TF.data()+TF.size(),[&I](int & a){a=I(a);});
//     for_each(b.data(),b.data()+b.size(),[&I](int & a){a=I(a);});
//     igl::slice(Eigen::MatrixXd(TV),J,1,TV);
//     igl::slice(Eigen::MatrixXd(W),J,1,W);
//     std::cout << "It's done!!" << std::endl;
//     v.push_back(W);
//     v.push_back(Vcontrol_static);
//     return v;
// }

// Eigen::MatrixXd pointReadFormat(std::string refPointPath, int numP){
//   Eigen::MatrixXd Vout(numP, 3);
//   std::ifstream inFile;
//   inFile.open(refPointPath.c_str());
//   int count = 0;
//   while(!inFile.eof()){
//     inFile >> Vout(count, 0) >> Vout(count, 1) >> Vout(count, 2);
//     count++;
//   }
//   inFile.close();
//   return Vout;
// }

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cout << "Usage: " << argv[0] << " plyMeshFilePath inPointsPath outPointsPath outMeshFilePath" << std::endl;
        return EXIT_FAILURE;
    }
    std::string inMesh  = std::string(argv[1]);
    std::string inPoint   = std::string(argv[2]);
	std::string outMesh  = std::string(argv[1]);
    std::string outPoint   = std::string(argv[2]);

	// read the points
	Vcontrol_static = pointReadFormat(repPointpath, numParticles);
	// Compute the Warp Matrix

	
	// Compute Transformation
	
	// Save Output Mesh

    return 0;
}
