/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
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

#include "gtest/gtest.h"

#include <fstream>

#include <LeMonADE/io/Parser.h>

//this test simply checks if the parser finds the
//commands from a file,returns them in the correct format
TEST(ParserTest,FindRead)
{
  std::ifstream file;
  file.open("tests/parserTest.test");
  Parser parser(file);

  EXPECT_EQ("!firstCommand",parser.findRead());
  EXPECT_EQ("!secondCommand",parser.findRead());
  EXPECT_EQ("#!thirdCommand",parser.findRead());
  EXPECT_EQ("endoffile",parser.findRead());
}