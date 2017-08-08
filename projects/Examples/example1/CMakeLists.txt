cmake_minimum_required (VERSION 2.6)

include_directories ("${LEMONADE_DIR}/include/")

add_executable(Example1 ex1main.cpp)
message("project_build_static = ${project_build_static}")
if( ( NOT DEFINED project_build_static ) OR ( project_build_static ) )
    target_link_libraries(Example1 libLeMonADE-static)
else()
    link_directories ("${CMAKE_BINARY_DIR}/lib/")
    target_link_libraries(Example1 LeMonADE)
endif()
