cmake_minimum_required (VERSION 2.6)

include_directories ("${LEMONADE_DIR}/include/")
link_directories ("${LEMONADE_DIR}/lib/")

add_executable(Example1 ex1main.cpp)
target_link_libraries(Example1 LeMonADE)
