cmake_minimum_required (VERSION 2.6)


include_directories ("${LEMONADE_DIR}/include/")
link_directories ("${LEMONADE_DIR}/lib/")


add_executable(Example4 ex4main.cpp)
target_link_libraries(Example4 LeMonADE)
