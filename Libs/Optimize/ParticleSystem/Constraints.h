#pragma once

#include "Constraint.h"
#include "PlaneConstraint.h"
#include "SphereConstraint.h"
#include "FreeFormConstraint.h"
#include <vector>
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_cross.h"
#include "Eigen/Dense"
#include "itkPoint.h"

namespace itk
{

class Constraints{
public:
    Constraints(){
        planeConsts = new std::vector<PlaneConstraint>();
        sphereConsts = new std::vector<SphereConstraint>();
        freeFormConsts = new std::vector<FreeFormConstraint>();
    }

    ~Constraints(){
        delete planeConsts;
        delete sphereConsts;
        delete freeFormConsts;
    }

  // Set constraints
  void addPlane(const vnl_vector<double> &a, const vnl_vector<double> &b,const vnl_vector<double> &c);
  void addSphere(const vnl_vector_fixed<double, DIMENSION> &v, double r);
  void setFreeFormConstraint(std::string filename);

  // Transforms
  bool transformConstraints(const vnl_matrix_fixed<double, 4, 4> &Trans);
  bool transformPlanes(const vnl_matrix_fixed<double, 4, 4> &Trans);

  // Apply functions
  bool applyBoundaryConstraints(vnl_vector_fixed<double, 3> &gradE, const Point<double, 3> &pos);
  bool applyPlaneConstraints(vnl_vector_fixed<double, 3> &gradE, const Point<double, 3> &pos);

  // Set & get defined
  void setCuttingPlaneDefined(bool value){cuttingPlaneDefined = value;}
  bool getCuttingPlaneDefined(){return cuttingPlaneDefined;}

  // Write constraints
  bool writePlanes(std::string filename){return true;}
  bool writeSpheres(std::string filename){return true;}
  bool writeFreeFormConstraint(std::string filename){return true;}

  // Is defined? functions
  bool IsCuttingPlaneDefined() const {if(planeConsts->size() > 0) return true; return false;}
  bool IsCuttingSphereDefined() const {if(sphereConsts->size() > 0) return true; return false;}

  // Plane constraint
  std::vector<PlaneConstraint> *getPlaneConstraints(){return planeConsts;}

protected:
  std::vector<PlaneConstraint> *planeConsts;
  std::vector<SphereConstraint> *sphereConsts;
  std::vector<FreeFormConstraint> *freeFormConsts;

  bool cuttingPlaneDefined;

private:
  Eigen::Vector3d projectOntoLine(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d p);
  Eigen::Vector3d linePlaneIntersect(Eigen::Vector3d n, Eigen::Vector3d p0, Eigen::Vector3d l0, Eigen::Vector3d l);
};


}
