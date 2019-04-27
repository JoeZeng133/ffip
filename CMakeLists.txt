cmake_minimum_required(VERSION 3.13)
project(ffip)

find_package(CGAl 4.13 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(MPI REQUIRED)
find_package(Python3 REQUIRED Development)
find_package(HDF5 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(HighFive REQUIRED)
find_package(Catch2 REQUIRED)

set(HEADER_FILES
"include/dft.hpp"
"include/fields.hpp"
"include/geometry.hpp"
"include/medium.hpp"
"include/simulation.hpp"
"include/structure.hpp"
"include/utility.hpp"
)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(app)