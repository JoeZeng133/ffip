project(script)

add_executable(Temp_Test Temp_Test.cpp)
target_link_libraries(Temp_Test ffip)

add_executable(std_config std_config.cpp)
target_link_libraries(std_config ffip)

add_executable(hd5_test hd5_test.cpp)
target_link_libraries(hd5_test ${HDF5_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_executable(json_test json_test.cpp)
target_link_libraries(json_test PRIVATE nlohmann_json::nlohmann_json)
