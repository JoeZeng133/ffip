project(ffip)
set (CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
simulation.cpp
utility.cpp
medium.cpp
chunk.cpp
geometry.cpp
analysis.cpp
chunk.cpp
source.cpp
)

add_library(ffip ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(ffip ${Boost_LIBRARIES})
