# Include directories for library code.
include_directories(.)
# Find sources.
file(GLOB parmetis_sources *.c)
# Create libparmetis
add_library(parmetis ${ParMETIS_LIBRARY_TYPE} ${parmetis_sources})
# Link with metis and MPI libraries.
target_link_libraries(parmetis metis ${MPI_LIBRARIES})
set_target_properties(parmetis PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

install(TARGETS parmetis
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib)
