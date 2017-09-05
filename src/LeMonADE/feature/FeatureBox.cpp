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

#include <LeMonADE/feature/FeatureBox.h>

/***********************************************************************/
/**
 *@file
 *@brief Implementation of class FeatureBox
 * */
/***********************************************************************/

FeatureBox::FeatureBox() :boxX(0)
        ,boxY(0)
        ,boxZ(0)
        ,periodicX(true)
        ,periodicY(true)
        ,periodicZ(true)
	,periodicInitX(false)
	,periodicInitY(false)
	,periodicInitZ(false)
{}


/************************************************************************
 * getter and setter for box dimension and periodicity
 * *********************************************************************/
uint64_t FeatureBox::getNumberOfLatticeSites() const {
    return boxX*boxY*boxZ;
}

void FeatureBox::setBoxX(int x) {
    boxX=x;
}
void FeatureBox::setBoxY(int y) {
    boxY=y;
}
void FeatureBox::setBoxZ(int z) {
    boxZ=z;
}
int FeatureBox::getBoxX() const {
    this->assertBoxSizeSetX();
    return boxX;
}
int FeatureBox::getBoxY() const {
    this->assertBoxSizeSetY();
    return boxY;
}
int FeatureBox::getBoxZ() const {
    this->assertBoxSizeSetZ();
    return boxZ;
}
void FeatureBox::setPeriodicX(bool x) {
    periodicX=x;
    periodicInitX=true;
    //periodicity_set.setX(true);
}
void FeatureBox::setPeriodicY(bool y) {
    periodicY=y;
    periodicInitY=true;
    //periodicity_set.setY(true);
}
void FeatureBox::setPeriodicZ(bool z) {
    periodicZ=z;
    periodicInitZ=true;
    //periodicity_set.setZ(true);
}
bool FeatureBox::isPeriodicX() const {
    this->assertPeriodicitySetX();
    return periodicX;
}
bool FeatureBox::isPeriodicY() const {
    this->assertPeriodicitySetY();
    return periodicY;
}
bool FeatureBox::isPeriodicZ() const {
    this->assertPeriodicitySetZ();
    return periodicZ;
}

bool FeatureBox::isCubic() const {
    if ((boxX == boxY) && (boxY==boxZ)) return true;
    else return false;
}


/************************************************************************
 * functions for asserting that variables have been set in bfm file
 * *********************************************************************/

/**
 * @throw <std::runtime_error> if size in x-direction (Length) has not been set.
 */
void FeatureBox::assertBoxSizeSetX() const {
    if ( boxX <= 0 ) throw  std::runtime_error("FeatureBox :: Box size has not been set in X direction or smaller than 0.");
}

/**
 * @throw <std::runtime_error> if size in y-direction (Width) has not been set.
 */
void FeatureBox::assertBoxSizeSetY() const {
    if ( boxY <= 0 ) throw  std::runtime_error("FeatureBox :: Box size has not been set in Y direction or smaller than 0.");
}

/**
 * @throw <std::runtime_error> if size in z-direction (Height) has not been set.
 */
void FeatureBox::assertBoxSizeSetZ() const {
    if ( boxZ <= 0 ) throw std::runtime_error ("FeatureBox :: Box size has not been set in Z direction or smaller than 0.");
}

/**
 * @throw <std::runtime_error> if size in all direction (Length, Width, Height) has not been set.
 */
void FeatureBox::assertBoxSizeSet() const {
    assertBoxSizeSetX();
    assertBoxSizeSetY();
    assertBoxSizeSetZ();
}

/**
 * @throw <std::runtime_error> if p.b.c in x-direction has not been set.
 */
void FeatureBox::assertPeriodicitySetX() const {
    if (!periodicInitX) throw std::runtime_error("FeatureBox :: Periodicity has not been set in X direction.");
}

/**
 * @throw <std::runtime_error> if p.b.c in y-direction has not been set.
 */
void FeatureBox::assertPeriodicitySetY() const {
    if (!periodicInitY) throw std::runtime_error("FeatureBox :: Periodicity has not been set in Y direction.");
}

/**
 * @throw <std::runtime_error> if p.b.c in z-direction has not been set.
 */
void FeatureBox::assertPeriodicitySetZ() const {
    if (!periodicInitZ) throw std::runtime_error("FeatureBox :: Periodicity has not been set in Z direction.");
}

/**
 * @throw <std::runtime_error> if p.b.c in all direction has not been set.
 */
void FeatureBox::assertPeriodicitySet() const {
    assertPeriodicitySetX();
    assertPeriodicitySetY();
    assertPeriodicitySetZ();
}
