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
 * @brief Tests for the class UpdaterLipidsRandomInitializer
 *
 * @author Ankush Checkervarty
 * @date 01.08.2019
 * */
/*****************************************************************************/


#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/updater/UpdaterLipidsRandomInitializer.h>
#include "TesterFunctions.h"

using namespace std;




/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class TestUpdaterLipidsRandomInitializer: public ::testing::Test{
public:

  //redirect cout output
  virtual void SetUp(){
    originalBuffer=cout.rdbuf();
    cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    cout.rdbuf(originalBuffer);
  };
  
private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;

};

//This tests are similar to tests for UpdaterLipidsCreator, the only difference
//here is that Orientations would random instead only towards z axis.
TEST_F(TestUpdaterLipidsRandomInitializer, TestSimpleRandomisedSystem)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;
  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
  int32_t box_volume=64*64*64;
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  EXPECT_NO_THROW(ingredients.synchronize());

  UpdaterLipidsRandomInitializer<IngredientsType> first(ingredients,2,0.1);
  
  int32_t defaultLipidTail=5;
  int32_t lipidsLength=2*defaultLipidTail+3;
  int32_t numLipids=2;
  
  
  TesterFunctions<IngredientsType> TF;  
  
  // first execution
  EXPECT_NO_THROW(first.initialize());

  EXPECT_NO_THROW(ingredients.synchronize());
  
  EXPECT_EQ((numLipids*lipidsLength),TF.getTotalNonsolventMonomers(ingredients));

  EXPECT_EQ(int32_t((0.1*box_volume-numLipids*lipidsLength*8)/8),TF.getTotalsolventMonomers(ingredients));

  
  // check (default) monomer tags
  // first lipid.
  EXPECT_EQ(1,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(1,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[3].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[10].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[12].getAttributeTag());
  
  // second lipid.
  EXPECT_EQ(1,ingredients.getMolecules()[13].getAttributeTag());
  EXPECT_EQ(1,ingredients.getMolecules()[15].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[17].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[20].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[25].getAttributeTag());
  
  // solvent .
  EXPECT_EQ(5, ingredients.getMolecules()[30].getAttributeTag());
  EXPECT_EQ(5, ingredients.getMolecules()[32].getAttributeTag());
  EXPECT_EQ(5, ingredients.getMolecules()[27].getAttributeTag());
  EXPECT_EQ(5, ingredients.getMolecules()[28].getAttributeTag());

  // check connectivity
  //first lipid.
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(0));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(1));
  EXPECT_EQ(3, ingredients.getMolecules().getNumLinks(2));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(5));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(7));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(9));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(12));
  
  //second lipid.
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(13));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(14));
  EXPECT_EQ(3, ingredients.getMolecules().getNumLinks(15));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(18));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(20));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(22));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(25));
 
  //solvent. 
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(30));
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(32));
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(27));
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(28)); 

  //Orientations of lipids.
  std::vector<VectorFloat3> Orientations = TF.getOrientations(ingredients,defaultLipidTail);

  EXPECT_EQ(2, Orientations.size());
  
  //Different from Lipids creator.
  EXPECT_GE(fabs(Orientations[0].getX()),0);
  EXPECT_GE(fabs(Orientations[0].getY()),0);
  EXPECT_GE(fabs(Orientations[0].getZ()),0);
  
  EXPECT_GE(fabs(Orientations[1].getX()),0);
  EXPECT_GE(fabs(Orientations[1].getY()),0);
  EXPECT_GE(fabs(Orientations[1].getZ()),0);
  
}

TEST_F(TestUpdaterLipidsRandomInitializer, TestMixedlipids)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;
  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
  int32_t box_volume=64*64*64;
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  EXPECT_NO_THROW(ingredients.synchronize());

  UpdaterLipidsRandomInitializer<IngredientsType> first_1(ingredients,4,0.1,5,7,0.5);
  
  int32_t longerLipidsLength=2*5+3;
  int32_t shorterLipidsLength=2*7+3;
  int32_t longerTail=7;
  int32_t shorterTail=5;
  //number of longerlipids.
  int32_t n1=2;
  //number of shorterlipids.
  int32_t n2=2;
  
  
  TesterFunctions<IngredientsType> TF;
  // first execution
  EXPECT_NO_THROW(first_1.initialize());

  EXPECT_NO_THROW(ingredients.synchronize());
  
  
  EXPECT_EQ((n1*longerLipidsLength + n2*shorterLipidsLength),TF.getTotalNonsolventMonomers(ingredients));

  EXPECT_EQ(int32_t((0.1*box_volume-(n1*longerLipidsLength + n2*shorterLipidsLength)*8)/8),TF.getTotalsolventMonomers(ingredients));
  
  //check (default) monomer tags  
  //shorter lipids(3 head monomers here).
  int32_t numMonHead1=TF.getNumMonWithAttribute(ingredients,1);
  EXPECT_EQ(3*n2,numMonHead1);

  int32_t numMonTail1=TF.getNumMonWithAttribute(ingredients,2);
  EXPECT_EQ(shorterTail*2*n2,numMonTail1);

  //longer lipids(3 head monomers).
  int32_t numMonHead2=TF.getNumMonWithAttribute(ingredients,3);
  EXPECT_EQ(3*n1,numMonHead2);

  int32_t numMonTail2=TF.getNumMonWithAttribute(ingredients,4);
  EXPECT_EQ(longerTail*2*n1,numMonTail2);

  //check connectivity
  //All the heads with three connections.
  int32_t numHeads3links=TF.getNumMonWithLinks(ingredients,3);
  EXPECT_EQ(1*(n1+n2),numHeads3links);
  
  //Two tail ends connections plus one head monomer.  
  int32_t numTailsEnds=TF.getNumMonWithLinks(ingredients,1);
  EXPECT_EQ((2+1)*(n1+n2),numTailsEnds);
  
  
  //Orientations of shorter lipids.
  std::vector<VectorFloat3> Orientations = TF.getOrientations(ingredients,shorterTail,false);
  
  EXPECT_EQ(2, Orientations.size());

  EXPECT_GE(fabs(Orientations[0].getX()),0);
  EXPECT_GE(fabs(Orientations[0].getY()),0);
  EXPECT_GE(fabs(Orientations[0].getZ()),0);
  
  EXPECT_GE(fabs(Orientations[1].getX()),0);
  EXPECT_GE(fabs(Orientations[1].getY()),0);
  EXPECT_GE(fabs(Orientations[1].getZ()),0);
  
  //Orientations of longer lipids.
  Orientations = TF.getOrientations(ingredients,longerTail,false);
  
  EXPECT_EQ(2, Orientations.size());

  EXPECT_GE(fabs(Orientations[0].getX()),0);
  EXPECT_GE(fabs(Orientations[0].getY()),0);
  EXPECT_GE(fabs(Orientations[0].getZ()),0);
  
  EXPECT_GE(fabs(Orientations[1].getX()),0);
  EXPECT_GE(fabs(Orientations[1].getY()),0);
  EXPECT_GE(fabs(Orientations[1].getZ()),0);
}
