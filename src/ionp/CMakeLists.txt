set( PROJNAME IONP )

set( SRCS main.cpp
	  Aggregate.cpp
	  Ionp.cpp
	  Particle.cpp
	  Point3D.cpp
	  Proton.cpp
	  SimMatrix.cpp
	  Vector3D.cpp
	  Data.cpp
    RandomList.cpp
    MRsequence.cpp
	  SimSpace.cpp
   )

#find_package( spdlog REQUIRED )
#find_package( Boost COMPONENTS filesystem date_time REQUIRED )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable( ${PROJNAME} ${SRCS} )

target_link_system_libraries(
  ${PROJNAME}
  PRIVATE
)
target_link_libraries(
  ${PROJNAME}
  PRIVATE project_options
          project_warnings
)

target_include_directories(
  ${PROJNAME}
  PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include"
          "${CMAKE_SOURCE_DIR}/include"
          "${CMAKE_SOURCE_DIR}/external"
)
