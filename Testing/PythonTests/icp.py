import os
import sys
from shapeworks import *

def icpTest():
  imgSource = Image(os.environ["DATA"] + "/smooth1.nrrd")
  imgTarget = Image(os.environ["DATA"] + "/smooth2.nrrd")
  xform = imgSource.createTransform(imgTarget, TransformType.IterativeClosestPoint, 1.0, 5)
  imgSource.applyTransform(xform, imgTarget.origin(), imgTarget.dims(), imgTarget.spacing(), imgTarget.coordsys(), InterpolationType.NearestNeighbor)

  compareImg = Image(os.environ["DATA"] + "/icp.nrrd")

  return imgSource == compareImg

try:
  icpTest()
except ValueError:
  print("icpTest failed")
  sys.exit(1)
