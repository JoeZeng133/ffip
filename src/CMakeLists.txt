project(ffip)

set(SOURCE_FILES
simulation.cpp
utility.cpp
medium.cpp
fields.cpp
geometry.cpp
dft.cpp
structure.cpp
)

add_library(ffip STATIC ${SOURCE_FILES} ${HEADER_FILES})
set_property(TARGET ffip PROPERTY CXX_STANDARD 17)
target_link_libraries(ffip PRIVATE nlohmann_json::nlohmann_json ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} CGAL::CGAL)

