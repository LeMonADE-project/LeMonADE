# ---------------------------------------------------------------------------------
#     ooo      L   attice-based  |
#   o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
#  o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
# oo---0---oo  A   lgorithm and  |
#  o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
#   o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
#     ooo                        | 
# ---------------------------------------------------------------------------------
#
# This file is part of LeMonADE.
# 
# LeMonADE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# LeMonADE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------
#
# Project Properties
#
CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (LeMonADE)
SET (APPLICATION_NAME "LeMonADE")
SET (APPLICATION_CODENAME "${PROJECT_NAME}")
SET (APPLICATION_COPYRIGHT_YEARS "2015")
SET (APPLICATION_VERSION_MAJOR 1)
SET (APPLICATION_VERSION_MINOR 0)
SET (APPLICATION_VERSION_PATCH 0)
SET (APPLICATION_VERSION_TYPE SNAPSHOT)
SET (APPLICATION_VERSION_STRING "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}-${APPLICATION_VERSION_TYPE}")
SET (APPLICATION_ID "${APPLICATION_VENDOR_ID}.${PROJECT_NAME}")


#
# Debugging Options
#
SET (CMAKE_VERBOSE_MAKEFILE 0) # Use 1 for debugging, 0 for release

#
# Project Output Paths
#
SET (LEMONADE_DIR ${PROJECT_SOURCE_DIR})

#
# Check if LEMONADE_DIR is given
#
if (DEFINED INSTALLDIR_LEMONADE)
    message("INSTALLDIR_LEMONADE set to "${INSTALLDIR_LEMONADE})
else()
    SET (INSTALLDIR_LEMONADE "${LEMONADE_DIR}/build/")
    message("INSTALLDIR_LEMONADE set automatically to "${INSTALLDIR_LEMONADE}". Use -DINSTALLDIR_LEMONADE=/path/to/install/LeMonADE/ for other destination")
endif()

SET (EXECUTABLE_OUTPUT_PATH "${INSTALLDIR_LEMONADE}/bin")
SET (LIBRARY_OUTPUT_PATH "${INSTALLDIR_LEMONADE}/lib")

#
# Project Search Paths
#
LIST (APPEND CMAKE_PREFIX_PATH "${LEMONADE_DIR}")
INCLUDE_DIRECTORIES("${LEMONADE_DIR}/include")


INSTALL(DIRECTORY "${LEMONADE_DIR}/include/LeMonADE/" DESTINATION  "${INSTALLDIR_LEMONADE}/include")
#
# Add Build Targets
#
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(projects)

option(LEMONADE_TESTS "Build the test" OFF)

if(LEMONADE_TESTS)
    add_subdirectory(tests)
endif(LEMONADE_TESTS)

#
# Add Install Targets
#

#
# Add Documentation Targets
#
SET (DOC_FILE_PATH "${LEMONADE_DIR}/docs")

FIND_PACKAGE (Doxygen)
IF (DOXYGEN_FOUND)
    MESSAGE("Build documentation with: make docs")
    IF (EXISTS ${DOC_FILE_PATH})
        MESSAGE("Existing File documentation with doxygen")
        ADD_CUSTOM_TARGET(
            docs
            ${DOXYGEN_EXECUTABLE} doxygen.conf
            WORKING_DIRECTORY ${DOC_FILE_PATH}
            COMMENT "Generating doxygen project documentation." VERBATIM
        )
    ELSE (EXISTS ${DOC_FILE_PATH})
        ADD_CUSTOM_TARGET(docs COMMENT "Doxyfile not found. Please generate a doxygen configuration file to use this target." VERBATIM)
    ENDIF (EXISTS ${DOC_FILE_PATH})
ELSE (DOXYGEN_FOUND)
    ADD_CUSTOM_TARGET(docs COMMENT "Doxygen not found. Please install doxygen to use this target." VERBATIM)
ENDIF (DOXYGEN_FOUND)

