ADD_EXECUTABLE ( SU2_SOL
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/output_structure.cpp
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/output_cgns.cpp
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/output_paraview.cpp
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/output_tecplot.cpp
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/solver_structure.cpp
  ${CMAKE_SOURCE_DIR}/SU2_CFD/src/variable_structure.cpp
  src/SU2_SOL.cpp
  )

TARGET_LINK_LIBRARIES ( SU2_SOL
  SU2
  ${JSONCPP_json_LIBRARY}
  ${CGNS_cgns_LIBRARY}
  ${Metis_metis_LIBRARY}
  ${MPI_EXTRA_LIBRARY}
  ${MPI_LIBRARY}
  )
