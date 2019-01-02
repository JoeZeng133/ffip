project(script)

add_executable(Temp_Test Temp_Test.cpp)
target_link_libraries(Temp_Test ffip)

add_executable(std_config std_config.cpp)
target_link_libraries(std_config ffip)
