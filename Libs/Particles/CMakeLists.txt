# Particles library

set(Particles_sources
 ParticleSystem.cpp
 ParticleShapeStatistics.cpp
 ShapeEvaluation.cpp
 ReconstructSurface.cpp
 ParticleNormalEvaluation.cpp
 )

set(Particles_headers
 ParticleSystem.h
 ParticleShapeStatistics.h
 itkParticlePositionReader.h
 itkParticlePositionWriter.h
 EvaluationUtil.h
 ShapeEvaluation.h
 ReconstructSurface.h
 ParticleNormalEvaluation.h
 )

add_library(Particles STATIC
  ${Particles_sources}
  ${Particles_headers}
  )

target_include_directories(Particles PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(Particles PUBLIC
  Mesh
  tinyxml
  Eigen3::Eigen
  Project
  )

# set
set_target_properties(Particles PROPERTIES PUBLIC_HEADER
  "${Particles_headers}"
  )

install(TARGETS Particles
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
