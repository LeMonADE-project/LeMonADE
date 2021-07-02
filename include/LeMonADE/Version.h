/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2020 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

/*
 * AUTO-GENERATION WARNING:
 *     "Version.h" has been automatically generated from "Version.h.in".
 *     DO NOT edit "Version.h", as any changes will be automatically discarded.
 */


#ifndef LEMONADE_VERSION_H
#define LEMONADE_VERSION_H

#define APPLICATION_NAME               "LeMonADE"
#define APPLICATION_CODENAME           "LeMonADE"
#define APPLICATION_COPYRIGHT_YEARS    "2013-2021"
#define APPLICATION_VERSION_MAJOR      2
#define APPLICATION_VERSION_MINOR      2
#define APPLICATION_VERSION_PATCH      2
#define APPLICATION_VERSION_TYPE       "Release"
#define APPLICATION_VERSION_STRING     "2.2.2-(C)2013-2021-Release"
#define APPLICATION_ID                 "LeMonADE.LeMonADE"

#ifndef APPLICATION_VERSION_MAJOR
#   define APPLICATION_VERSION_MAJOR 0
#endif

#ifndef APPLICATION_VERSION_MINOR
#   define APPLICATION_VERSION_MINOR 0
#endif

#ifndef APPLICATION_VERSION_PATCH
#   define APPLICATION_VERSION_PATCH 0
#endif

#ifndef APPLICATION_VERSION_TYPE
#   define APPLICATION_VERSION_TYPE "RELEASE"
#endif

#define APPLICATION_VERSION_MAJOR_MINOR      2.2

const double LEMONADE_VERSION = APPLICATION_VERSION_MAJOR_MINOR;

#endif
