cmake_minimum_required(VERSION 2.8)

include_directories ("${INSTALLDIR_LEMONADE}/include/")
link_directories ("${INSTALLDIR_LEMONADE}/lib/")

add_executable(SimpleSimulator main.cpp)

target_link_libraries(SimpleSimulator LeMonADE)
