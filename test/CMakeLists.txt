add_executable(utility_test utility_test.cpp)
set_property(TARGET utility_test PROPERTY CXX_STANDARD 17)

add_executable(mpi_test mpi_test.cpp)
set_property(TARGET mpi_test PROPERTY CXX_STANDARD 17)
target_link_libraries(mpi_test ${MPI_CXX_LIBRARIES})