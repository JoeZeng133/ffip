add_executable(run_sim_json run_sim_json.cpp)
set_property(TARGET run_sim_json PROPERTY CXX_STANDARD 17)
target_link_libraries(run_sim_json PUBLIC ffip)