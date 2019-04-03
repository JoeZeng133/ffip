cmake_minimum_required(VERSION 3.13)
set(Boost_DEBUG=ON)
set(CMAKE_BUILD_TYPE Release)

find_package(CGAl 4.13 REQUIRED) #it has to go before boost I think
find_package(Boost 1.68 REQUIRED COMPONENTS system python37)
find_package(MPI REQUIRED)
find_package(Python3 REQUIRED Development)
find_package(HDF5 REQUIRED)


file(GLOB HEADER_FILES
"include/*.hpp"
"src/*.cpp"
)

include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_DIRS})
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(py)
add_subdirectory(test)