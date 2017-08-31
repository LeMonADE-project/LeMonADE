cmake_minimum_required (VERSION 2.6)

include_directories ("${LEMONADE_DIR}/include/")

add_executable(Example1 ex1main.cpp)
# This condition was in principle copied from LeMonADE/src/CMakeLists.txt,
# because only if the static library was added as a compile target to CMake,
# it is possible to directly specify the library by its CMake project name
# instead of having to add a library path. By using the CMake target directly,
# CMake also automatically infers the depedency graph correctly when using
# make -j, i.e. target Example1 will only be compiled after target
# libLeMonADE-static finished successfully
if( ( NOT DEFINED project_build_static ) OR ( project_build_static ) )
    target_link_libraries(Example1 libLeMonADE-static)
else()
    link_directories ("${CMAKE_BINARY_DIR}/lib/")
    target_link_libraries(Example1 LeMonADE)
endif()
