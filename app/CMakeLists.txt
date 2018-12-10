project(script)

add_executable(Temp_Test Temp_Test.cpp)
target_link_libraries(Temp_Test ffip)

add_executable(Dipole_Test Dipole_Test.cpp)
target_link_libraries(Dipole_Test ffip)

add_executable(N2F_Test N2F_Test.cpp)
target_link_libraries(N2F_Test ffip)

add_executable(Plane_Wave_Test Plane_Wave_Test.cpp)
target_link_libraries(Plane_Wave_Test ffip)

add_executable(Mie_Test Mie_Test.cpp)
target_link_libraries(Mie_Test ffip)

add_executable(std_config std_config.cpp)
target_link_libraries(std_config ffip)
