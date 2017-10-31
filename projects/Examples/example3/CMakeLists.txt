cmake_minimum_required (VERSION 2.6)

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

if (NOT DEFINED LEMONADE_LIBRARY_DIR)
message("LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location")
endif()

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

add_executable(Example3 ex3main.cpp)
if( ( NOT DEFINED project_build_static ) OR ( project_build_static ) )
    # using the project name set in LeMonADE/src/CMakeLists.txt
    # instead of the library file name correctly watches out for dependencies
    # when trying to compile in parallel with make -j
    target_link_libraries(Example3 libLeMonADE-static)
else()
    link_directories ("${CMAKE_BINARY_DIR}/lib/")
    target_link_libraries(Example3 LeMonADE)
endif()
