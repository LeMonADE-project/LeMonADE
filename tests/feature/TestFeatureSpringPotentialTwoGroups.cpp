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
 * @brief test for Feature Spring Potential
 * @author Martin
 * @date 08.05.2019
 * */
/*****************************************************************************/
#include "gtest/gtest.h"


#include <iostream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>

#include <LeMonADE/feature/FeatureSpringPotentialTwoGroups.h>


class TestFeatureSpringPotentialTwoGroups: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureSpringPotentialTwoGroups) Features;
  typedef ConfigureSystem<VectorInt3,Features> ConfigType;
  typedef Ingredients < ConfigType> IngredientsType;
  IngredientsType ingredients;

  //dummy move class used to check response to unknown move type
  class UnknownMove:public MoveBase
  {
	  public:
		template<class IngredientsType> bool check(IngredientsType& ingredients) const
		{
		return ingredients.checkMove(ingredients,*this);
		}

		template<class IngredientsType> void apply(IngredientsType& ingredients)
		{
		ingredients.applyMove(ingredients,*this);
		}

		template <class IngredientsType> void init(const IngredientsType& ingredients){};
  };
  
  /* suppress cout output for better readability -->un-/comment here: */
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
  /* ** */
};


TEST_F(TestFeatureSpringPotentialTwoGroups,Basics){
  // test default constructor
  FeatureSpringPotentialTwoGroups feature;
  EXPECT_DOUBLE_EQ(0.0,feature.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.0,feature.getSpringConstant());
  
  //test constructor in ingredients
  EXPECT_DOUBLE_EQ(0.0,ingredients.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.0,ingredients.getSpringConstant());

  //setup a small system
  ingredients.setSpringConstant(0.25);
  ingredients.setEquilibriumLength(4);
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(0,0,0);

  // check the spring parameters and the monomer extension
  EXPECT_DOUBLE_EQ(4.00,ingredients.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.25,ingredients.getSpringConstant());
  EXPECT_EQ(0,ingredients.getMolecules()[0].getMonomerGroupTag());

  // check the monomer extension setter and getter
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(1));
  EXPECT_EQ(1,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(2));
  EXPECT_EQ(2,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(0));
  EXPECT_EQ(0,ingredients.getMolecules()[0].getMonomerGroupTag());

  // other values than 0,1,2 are not allowed
  EXPECT_ANY_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(3));

  // equality operator
  ingredients.modifyMolecules().addMonomer(2,0,0);
  EXPECT_NO_THROW(ingredients.modifyMolecules()[1].setMonomerGroupTag(1));
  EXPECT_FALSE( ingredients.getMolecules()[0].getMonomerGroupTag() == ingredients.getMolecules()[1].getMonomerGroupTag() );
  ingredients.modifyMolecules()[0].setMonomerGroupTag(ingredients.getMolecules()[1].getMonomerGroupTag());
  EXPECT_TRUE( ingredients.getMolecules()[0].getMonomerGroupTag() == ingredients.getMolecules()[1].getMonomerGroupTag() );
  
}

TEST_F(TestFeatureSpringPotentialTwoGroups,MoveChecks){
  //setup a small system with a rather high spring constant
  ingredients.setSpringConstant(8.0);
  ingredients.setEquilibriumLength(2);
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(0,0,0);
  ingredients.modifyMolecules().addMonomer(8,0,0);

  // check the spring parameters and the monomer extension
  EXPECT_DOUBLE_EQ(2.0,ingredients.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(8.0,ingredients.getSpringConstant());
  EXPECT_EQ(0,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_EQ(0,ingredients.getMolecules()[1].getMonomerGroupTag());

  // check the monomer extension setter and getter
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(1));
  EXPECT_EQ(1,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_NO_THROW(ingredients.modifyMolecules()[1].setMonomerGroupTag(2));
  EXPECT_EQ(2,ingredients.getMolecules()[1].getMonomerGroupTag());

  // perform some moves and check the resulting center to center distance
  // **************   check base move   **************
  UnknownMove basemove;
  basemove.init(ingredients);
  EXPECT_TRUE(basemove.check(ingredients));
  //should change nothing
  basemove.apply(ingredients);
  EXPECT_EQ(VectorInt3(0,0,0),ingredients.getMolecules()[0]);
  EXPECT_EQ(VectorInt3(8,0,0),ingredients.getMolecules()[1]);
  EXPECT_EQ(1,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_EQ(2,ingredients.getMolecules()[1].getMonomerGroupTag());
  ingredients.synchronize();

  //MoveLocalSc scmove;
  MoveLocalSc scmove;
  
  for(uint32_t i=0;i<100;i++){
    scmove.init(ingredients);
    if(scmove.check(ingredients)){
      scmove.apply(ingredients);
    }
  }
  // monomers should be close
  std::cout << (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength() <<std::endl;
  EXPECT_TRUE(abs( abs( (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength()-2 ) ) < 2 );

  // add a nonaffected monomer and check the movement
  ingredients.modifyMolecules().addMonomer(8,8,8);
  int32_t counter(0);
  VectorInt3 referencePosition(ingredients.getMolecules()[2]);
  ingredients.synchronize();

  for(uint32_t i=0;i<100;i++){
    scmove.init(ingredients);
    if(scmove.check(ingredients)){
      scmove.apply(ingredients);
      if( scmove.getIndex() == 2 ){
        EXPECT_FALSE(referencePosition == ingredients.getMolecules()[2]);
	      referencePosition=ingredients.getMolecules()[2];
        counter++;
      }
    }
  }
  // monomers should still be close
  std::cout << (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength() <<std::endl
  <<", move counter = "<<counter<<std::endl;
  EXPECT_TRUE(abs( abs( (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength()-2 ) ) < 2 );
  // there should have been some moves of the unaffected monomer
  EXPECT_TRUE(counter>0);
  
  //MoveLocalScDiag scmovediag;
  MoveLocalScDiag scmovediag;
  //reste the positions
  ingredients.modifyMolecules()[0].modifyVector3D().setAllCoordinates(0,0,0);
  ingredients.modifyMolecules()[1].modifyVector3D().setAllCoordinates(13,8,9);
  
  for(uint32_t i=0;i<100;i++){
    scmovediag.init(ingredients);
    if(scmovediag.check(ingredients)){
      scmovediag.apply(ingredients);
    }
  }
  // monomers should be close
  std::cout << (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength() <<std::endl;
  EXPECT_TRUE(abs( abs( (ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength()-2 ) ) < 2 );
  
}

TEST_F(TestFeatureSpringPotentialTwoGroups,fileReadWrite){
  IngredientsType ingredientsWrite;
  IngredientsType ingredientsRead;

  ingredientsWrite.setBoxX(64);
  ingredientsWrite.setPeriodicX(true);
  ingredientsWrite.setBoxY(64);
  ingredientsWrite.setPeriodicY(true);
  ingredientsWrite.setBoxZ(64);
  ingredientsWrite.setPeriodicZ(true);
  ingredientsWrite.setSpringConstant(1.2);
  ingredientsWrite.setEquilibriumLength(4);
  ingredientsWrite.modifyMolecules().addMonomer(0,0,0);
  ingredientsWrite.modifyMolecules()[0].setMonomerGroupTag(1);
  ingredientsWrite.modifyMolecules().addMonomer(2,0,0);
  ingredientsWrite.modifyMolecules()[1].setMonomerGroupTag(1);
  ingredientsWrite.modifyMolecules().addMonomer(4,0,0);
  ingredientsWrite.modifyMolecules()[2].setMonomerGroupTag(2);
  ingredientsWrite.modifyMolecules().addMonomer(6,0,0);
  ingredientsWrite.modifyMolecules()[3].setMonomerGroupTag(2);
  ingredientsWrite.modifyMolecules().addMonomer(8,0,0);

  EXPECT_EQ(1,ingredientsWrite.getMolecules()[0].getMonomerGroupTag());
  EXPECT_EQ(1,ingredientsWrite.getMolecules()[1].getMonomerGroupTag());
  EXPECT_EQ(2,ingredientsWrite.getMolecules()[2].getMonomerGroupTag());
  EXPECT_EQ(2,ingredientsWrite.getMolecules()[3].getMonomerGroupTag());
  EXPECT_EQ(0,ingredientsWrite.getMolecules()[4].getMonomerGroupTag());

  //write to file and read back in
  AnalyzerWriteBfmFile<IngredientsType> outfile("tests/springPotentialTestOut.test",ingredientsWrite,AnalyzerWriteBfmFile<IngredientsType>::NEWFILE);
  outfile.initialize();
  outfile.execute();

  FileImport<IngredientsType> infile ("tests/springPotentialTestOut.test",ingredientsRead);

  //scan file for !mcs and read-in first frame
  infile.initialize();

  EXPECT_EQ(ingredientsRead.getMolecules()[0].getMonomerGroupTag(),1);
  EXPECT_EQ(ingredientsRead.getMolecules()[1].getMonomerGroupTag(),1);
  EXPECT_EQ(ingredientsRead.getMolecules()[2].getMonomerGroupTag(),2);
  EXPECT_EQ(ingredientsRead.getMolecules()[3].getMonomerGroupTag(),2);
  EXPECT_EQ(ingredientsRead.getMolecules()[4].getMonomerGroupTag(),0); /*has default value*/

  EXPECT_DOUBLE_EQ(1.2,ingredientsRead.getSpringConstant());
  EXPECT_DOUBLE_EQ(4.0,ingredientsRead.getEquilibriumLength());

    //remove the temporary file
  EXPECT_EQ(0,remove("tests/springPotentialTestOut.test"));
}
