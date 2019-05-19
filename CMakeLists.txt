cmake_minimum_required(VERSION 3.12)
project(ffip)

# rpaht configurations
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set(HDF5_PREFER_PARALLEL TRUE)

find_package(CGAL 4.11 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(MPI 2 REQUIRED)
find_package(HDF5 1.10 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(HighFive REQUIRED)

# get python package location from PYBIND11 generated PYTHON_EXECUTABLE
# set(cmd "from distutils.sysconfig import get_python_lib; print(get_python_lib(plat_specific=False,standard_lib=True))")
set(cmd "import sys; print(sys.path[-1])")

if(NOT Python_STDLIB)
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "${cmd}"
    OUTPUT_VARIABLE Python_STDLIB)
endif(NOT Python_STDLIB)

message("Python Library Location=" ${Python_STDLIB})
message("Python Executable="${PYTHON_EXECUTABLE})

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
add_subdirectory(test)
add_subdirectory(app)

install(TARGETS ffip run_sim_json
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
PUBLIC_HEADER DESTINATION include)

set(FFIP_LIB_FILES
"ffip/geom.py"
"ffip/__init__.py"
"ffip/materials.py"
"ffip/source.py"
"ffip/simulation.py"
)

# install python module into site-packages folder
install(TARGETS _ffip DESTINATION ${Python_STDLIB}/ffip)
install(FILES ${FFIP_LIB_FILES} DESTINATION ${Python_STDLIB}/ffip)