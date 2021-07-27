import os
import sys
from shapeworks import *

def curvatureTest1():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter()

  compareImg = Image(os.environ["DATA"] + "/curvature1.nrrd")

  return img.compare(compareImg)

try:
  curvatureTest1()
except ValueError:
  print("curvatureTest1 failed")
  sys.exit(1)

def curvatureTest2():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter(5)

  compareImg = Image(os.environ["DATA"] + "/curvature2.nrrd")

  return img.compare(compareImg)

try:
  curvatureTest2()
except ValueError:
  print("curvatureTest2 failed")
  sys.exit(1)

def curvaturefailTest():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  img.applyCurvatureFilter(-5)

  compareImg = Image(os.environ["DATA"] + "/curvaturefail.nrrd")

  return img.compare(compareImg)

try:
  curvaturefailTest()
  print("curvaturefailTest failed. The function cannot compute negative iterations.")
  sys.exit(1)
except TypeError:
  sys.exit(0)
