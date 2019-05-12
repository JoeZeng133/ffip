cmake_minimum_required(VERSION 3.13)
project(ffip)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(HDF5_PREFER_PARALLEL TRUE)

find_package(CGAL 4.11 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(MPI 2 REQUIRED)
find_package(HDF5 1.10 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(HighFive REQUIRED)
# find_package(Catch2 REQUIRED)

set(HEADER_FILES
"include/dft.hpp"
"include/fields.hpp"
"include/geometry.hpp"
"include/medium.hpp"
"include/simulation.hpp"
"include/structure.hpp"
"include/utility.hpp"
)

message(HDF5_VERSION is ${HDF5_VERSION})

add_subdirectory(src)
# add_subdirectory(test)
add_subdirectory(app)

install(TARGETS ffip run_sim_json
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
PUBLIC_HEADER DESTINATION include)

install(TARGETS _ffip DESTINATION ${PROJECT_SOURCE_DIR}/scripts/ffip)