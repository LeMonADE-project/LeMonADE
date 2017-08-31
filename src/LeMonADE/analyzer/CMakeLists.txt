## ----------------------------------------------------------------------------------
##     ooo      L   attice-based  |
##   o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
##  o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
## oo---0---oo  A   lgorithm and  |
##  o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
##   o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
##     ooo                        |
## ----------------------------------------------------------------------------------
##
## This file is part of LeMonADE.
##
## LeMonADE is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LeMonADE is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/LeMonADE/analyzer/*.h
  )

ADD_LIBRARY(obj_analyzer OBJECT ${_src} ${_header})
