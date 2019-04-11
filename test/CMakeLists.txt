add_executable(utility_test utility_test.cpp)
set_property(TARGET utility_test PROPERTY CXX_STANDARD 17)

add_executable(mpi_practice mpi_practice.cpp)
set_property(TARGET mpi_practice PROPERTY CXX_STANDARD 17)
target_link_libraries(mpi_practice ${MPI_CXX_LIBRARIES})

add_executable(hdf5_practice hdf5_practice.cpp)
set_property(TARGET hdf5_practice PROPERTY CXX_STANDARD 17)
target_link_libraries(hdf5_practice ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(catch_practice catch_practice.cpp)
set_property(TARGET catch_practice PROPERTY CXX_STANDARD 17)