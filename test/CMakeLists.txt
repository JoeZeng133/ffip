# mpi practice
add_executable(mpi_practice mpi_practice.cpp)
set_property(TARGET mpi_practice PROPERTY CXX_STANDARD 17)
target_link_libraries(mpi_practice PUBLIC MPI::MPI_CXX)

#  hdf5 practice
add_executable(hdf5_practice hdf5_practice.cpp)

set_property(TARGET hdf5_practice PROPERTY CXX_STANDARD 17)

target_include_directories(hdf5_practice
PUBLIC ${HDF5_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

target_link_libraries(hdf5_practice 
PUBLIC MPI::MPI_CXX ${HDF5_LIBRARIES})

# catch2 practice
add_executable(catch_practice catch_practice.cpp)
set_property(TARGET catch_practice PROPERTY CXX_STANDARD 17)
target_link_libraries(catch_practice
PUBLIC Catch2::Catch2)

# json practice
add_executable(json_practice json_practice.cpp)
set_property(TARGET json_practice PROPERTY CXX_STANDARD 17)
target_link_libraries(json_practice PUBLIC nlohmann_json::nlohmann_json)

# simulation test
add_executable(simulation_test simulation_test.cpp)
set_property(TARGET simulation_test PROPERTY CXX_STANDARD 17)
target_link_libraries(simulation_test 
PUBLIC ffip)

# utility test
add_executable(utility_test utility_test.cpp)
set_property(TARGET utility_test PROPERTY CXX_STANDARD 17)