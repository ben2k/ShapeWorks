import os
import sys
from shapeworks import *

success = True

def distanceTest1():
  femur = Mesh(os.environ["DATA"] + "/femur.vtk")
  pelvis = Mesh(os.environ["DATA"] + "/pelvis.vtk")
  f2p_distance_and_ids = femur.vertexDistance(pelvis)
  p2f_distance_and_ids = pelvis.vertexDistance(femur)
  femur.setField("distance", f2p_distance_and_ids[0])
  pelvis.setField("distance", p2f_distance_and_ids[0])
  femur.setField("closestPoints", f2p_distance_and_ids[1])
  pelvis.setField("closestPoints", p2f_distance_and_ids[1])

  f2p = Mesh(os.environ["DATA"] + "/meshdistance_point_fwd.vtk")
  p2f = Mesh(os.environ["DATA"] + "/meshdistance_point_rev.vtk")

  return femur == f2p and pelvis == p2f

success &= utils.test(distanceTest1)

def distanceTest2():
  femur1 = Mesh(os.environ["DATA"] + "/m03_L_femur.ply")
  femur2 = Mesh(os.environ["DATA"] + "/m04_L_femur.ply")
  fwd_distance_and_ids = femur1.distance(femur2)
  rev_distance_and_ids = femur2.distance(femur1)
  femur1.setField("distance", fwd_distance_and_ids[0])
  femur2.setField("distance", rev_distance_and_ids[0])
  femur1.setField("closestCells", fwd_distance_and_ids[1])
  femur2.setField("closestCells", rev_distance_and_ids[1])

  fwd = Mesh(os.environ["DATA"] + "/meshdistance_cell_fwd.vtk")
  rev = Mesh(os.environ["DATA"] + "/meshdistance_cell_rev.vtk")

  return femur1 == fwd and femur2 == rev

success &= utils.test(distanceTest2)

sys.exit(not success)
