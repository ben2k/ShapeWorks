import os
import sys
from typing import Type
from shapeworks import *

def curvatureTest1():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter()

  compareImg = Image(os.environ["DATA"] + "/curvature1.nrrd")

  return img.compare(compareImg)

utils.test(curvatureTest1)

def curvatureTest2():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter(5)

  compareImg = Image(os.environ["DATA"] + "/curvature2.nrrd")

  return img.compare(compareImg)

utils.test(curvatureTest2)

def curvaturefailTest():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter(-5)

  compareImg = Image(os.environ["DATA"] + "/curvaturefail.nrrd")

  return img.compare(compareImg)

utils.expectException(curvaturefailTest, TypeError)
