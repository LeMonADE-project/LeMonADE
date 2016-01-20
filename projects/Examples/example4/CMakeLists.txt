cmake_minimum_required (VERSION 2.6)


include_directories ("${INSTALLDIR_LEMONADE}/include/")
link_directories ("${INSTALLDIR_LEMONADE}/lib/")


add_executable(Example4 ex4main.cpp)
target_link_libraries(Example4 LeMonADE)
