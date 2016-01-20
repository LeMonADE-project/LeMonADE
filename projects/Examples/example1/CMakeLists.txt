cmake_minimum_required (VERSION 2.6)

include_directories ("${INSTALLDIR_LEMONADE}/include/")
link_directories ("${INSTALLDIR_LEMONADE}/lib/")

add_executable(Example1 ex1main.cpp)
target_link_libraries(Example1 LeMonADE)
