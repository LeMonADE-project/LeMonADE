cmake_minimum_required (VERSION 2.6)

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

include_directories (${LEMONADE_INCLUDE_DIR})

add_executable(Example0 ex0main.cpp)

