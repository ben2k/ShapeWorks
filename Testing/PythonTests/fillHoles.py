import os
import sys
from shapeworks import *

def fillholesTest():
  mesh = Mesh(os.environ["DATA"] + "/femur.vtk")
  mesh.fillHoles()

  compareMesh = Mesh(os.environ["DATA"] + "/fillholes.vtk")

  return mesh == compareMesh

try:
  fillholesTest()
except ValueError:
  print("fillholesTest failed")
  sys.exit(1)
