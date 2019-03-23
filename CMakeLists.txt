cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

find_package(MPI COMPONENTS CXX REQUIRED)
find_package(HDF5 REQUIRED)
find_package(nlohmann_json REQUIRED)


file(GLOB HEADER_FILES
"include/*.hpp"
"highfive/*.hpp"
)

# message(This is CXX HD5F Include HDF5_CXX_INCLUDE_DIRS)
include_directories(${MPI_CXX_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/highfive)

add_subdirectory(src)
add_subdirectory(app)
