set(SOURCE_FILES
"simulation.cpp"
"utility.cpp"
"medium.cpp"
"fields.cpp"
"geometry.cpp"
"dft.cpp"
"structure.cpp"
)

add_library(ffip SHARED ${SOURCE_FILES})

set_property(TARGET ffip PROPERTY CXX_STANDARD 17)

target_include_directories(ffip
PUBLIC
    ${HDF5_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR} 
    $<TARGET_PROPERTY:HighFive,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(ffip 
PUBLIC 
    nlohmann_json::nlohmann_json 
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
    CGAL::CGAL
    $<TARGET_PROPERTY:HighFive,INTERFACE_LINK_LIBRARIES>
)



