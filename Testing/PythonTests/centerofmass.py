import os
import sys
from shapeworks import *

def centerofmassTest1():
  img = Image(os.environ["DATA"] + "/1x2x2.nrrd")
  xform = img.createTransform()
  img.applyTransform(xform)

  compareImg = Image(os.environ["DATA"] + "/centerofmass1.nrrd")

  return img.compare(compareImg)

try:
  centerofmassTest1()
except ValueError:
  print("centerofmassTest1 failed")
  sys.exit(1)

def centerofmassTest2():
  img = Image(os.environ["DATA"] + "/la-bin.nrrd")
  xform = img.createTransform()
  img.applyTransform(xform)

  compareImg = Image(os.environ["DATA"] + "/centerofmass2.nrrd")

  return img.compare(compareImg)

try:
  centerofmassTest2()
except ValueError:
  print("centerofmassTest2 failed")
  sys.exit(1)
