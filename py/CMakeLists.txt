add_library(greet SHARED greet.cpp)

set_property(TARGET greet PROPERTY CXX_STANDARD 17)
target_link_libraries(greet PRIVATE ${Boost_LIBRARIES} ${Python3_LIBRARIES} CGAL::CGAL)

set_target_properties(greet PROPERTIES PREFIX "")
if(APPLE)
    set_target_properties(greet PROPERTIES SUFFIX ".so")
endif()

