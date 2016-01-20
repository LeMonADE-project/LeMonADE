cmake_minimum_required (VERSION 2.6)

include_directories ("${INSTALLDIR_LEMONADE}/include/")
link_directories ("${INSTALLDIR_LEMONADE}/lib/")

add_executable(Example3 ex3main.cpp)
target_link_libraries(Example3 LeMonADE)