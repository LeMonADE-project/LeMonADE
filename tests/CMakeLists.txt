# Enable ExternalProject CMake module
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/third_party)

# Add gtest
ExternalProject_Add(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG release-1.6.0
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIR ${source_dir}/include)

# Library
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARY_PATH ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
set(GTEST_LIBRARY gtest)
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${GTEST_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${GTEST_LIBRARY_PATH} )
add_dependencies(${GTEST_LIBRARY} googletest)

# copy all .test to build directory
configure_file(attributesTest.test tests/attributesTest.test COPYONLY)
configure_file(featureBondsetTest.test tests/featureBondsetTest.test COPYONLY)
configure_file(featureBoxTest.test tests/featureBoxTest.test COPYONLY)
configure_file(fileImportTest.test tests/fileImportTest.test COPYONLY)
configure_file(fileImportTest2.test tests/fileImportTest2.test COPYONLY)
configure_file(fileImportTest3.test tests/fileImportTest3.test COPYONLY)
configure_file(molecules.test tests/molecules.test COPYONLY)
configure_file(parserTest.test tests/parserTest.test COPYONLY)
configure_file(readbfmfile.test tests/readbfmfile.test COPYONLY)
configure_file(sc64CorruptedBfm.test tests/sc64CorruptedBfm.test COPYONLY)
configure_file(WriteReadAdditionalBondsTest.bfm tests/WriteReadAdditionalBondsTest.bfm COPYONLY)
#configure_file(tmpSolventTest.sol tests/tmpSolventTest.sol COPYONLY)

# build and run LeMonaDE-tests
INCLUDE_DIRECTORIES("${source_dir}/include")
FILE (GLOB_RECURSE test_SRCS *.cpp *.cxx *.cc *.C *.c *.h *.hpp)
SET (test_LIBS ${PROJECT_NAME} ${GTEST_LIBRARY} pthread dl)
SET (test_BIN ${PROJECT_NAME}-tests)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/tests)


IF (NOT CMAKE_CROSSCOMPILING)
    LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
    ADD_EXECUTABLE(${test_BIN} ${test_SRCS})
    TARGET_LINK_LIBRARIES(${test_BIN} ${test_LIBS})
    MESSAGE("Running LeMonADE tests...")
    ADD_CUSTOM_TARGET(check ALL "${CMAKE_BINARY_DIR}/tests/${test_BIN}" DEPENDS ${test_BIN} COMMENT "Executing LeMonADE tests..." VERBATIM SOURCES ${test_SRCS})
    ADD_CUSTOM_TARGET(testrun "${CMAKE_BINARY_DIR}/tests/${test_BIN}" DEPENDS ${test_BIN} COMMENT "Executing LeMonADE tests..." VERBATIM SOURCES ${test_SRCS})
ENDIF (NOT CMAKE_CROSSCOMPILING)

