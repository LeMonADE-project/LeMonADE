cmake_minimum_required(VERSION 2.8)

include_directories ("${LEMONADE_DIR}/include/")
link_directories ("${LEMONADE_DIR}/lib/")

add_executable(SimpleSimulator main.cpp)

target_link_libraries(SimpleSimulator LeMonADE)
