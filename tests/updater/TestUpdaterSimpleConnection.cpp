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


/*****************************************************************************/
/**
 * @file
 * @brief Tests for UpdaterAbstractCreate
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/UpdaterSimpleConnection.h>

class TestUpdaterSimpleConnectionSc: public ::testing::Test{
public:

  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureConnectionSc, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;
  
  IngredientsType ingredients;
  
  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };

private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;

};

TEST_F(TestUpdaterSimpleConnectionSc, Conversion)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);
  ingredients.modifyMolecules().addMonomer(10,8,8);
  ingredients.modifyMolecules().addMonomer(10,10,8);
  ingredients.modifyMolecules().addMonomer(8,4,8);
  ingredients.modifyMolecules().addMonomer(8,6,8);
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().connect(3,4);

  ingredients.modifyMolecules()[0].setReactive(true);
  ingredients.modifyMolecules()[0].setNumMaxLinks(4);
  
  ingredients.modifyMolecules()[1].setReactive(true);
  ingredients.modifyMolecules()[1].setNumMaxLinks(2);
  
  ingredients.modifyMolecules()[2].setReactive(false);
  ingredients.modifyMolecules()[2].setNumMaxLinks(2);
  
  ingredients.modifyMolecules()[3].setReactive(true);
  ingredients.modifyMolecules()[3].setNumMaxLinks(2);
  
  ingredients.modifyMolecules()[4].setReactive(true);
  ingredients.modifyMolecules()[4].setNumMaxLinks(1);
  EXPECT_NO_THROW(ingredients.synchronize());
  
  UpdaterSimpleConnection<IngredientsType,MoveLocalSc,MoveConnectSc> update(ingredients,10);
  update.initialize();
  // double getConversion(){return (1.*NReactedSites)/(1.*NReactiveSites);}; // NReactedSites=4, NReactiveSites=8
  EXPECT_EQ( 0.5,update.getConversion());
  
  
}