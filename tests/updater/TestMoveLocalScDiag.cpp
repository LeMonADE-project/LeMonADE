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
 * @brief Tests for the class MoveLocalSc
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

class TestMoveLocalScDiag: public ::testing::Test{
public:
    typedef LOKI_TYPELIST_4(FeatureBoltzmann, FeatureMoleculesIO, FeatureFixedMonomers, FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<bool> >) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> IngredientsType;

    IngredientsType ingredients;
    const IngredientsType& getIngredients() const {return ingredients;}
    TestMoveLocalScDiag () {
     	//perpendicular moves
	steps[0]=VectorInt3(1,0,0);
	steps[1]=VectorInt3(-1,0,0);
	steps[2]=VectorInt3(0,1,0);
	steps[3]=VectorInt3(0,-1,0);
	steps[4]=VectorInt3(0,0,1);
	steps[5]=VectorInt3(0,0,-1);
	//diagonal moves
	steps[6]=VectorInt3(0,1,1);
	steps[7]=VectorInt3(0,-1,1);
	steps[8]=VectorInt3(0,1,-1);
	steps[9]=VectorInt3(0,-1,-1);
	steps[10]=VectorInt3(1,0,1);
	steps[11]=VectorInt3(1,0,-1);
	steps[12]=VectorInt3(-1,0,1);
	steps[13]=VectorInt3(-1,0,-1);
	steps[14]=VectorInt3(1,1,0);
	steps[15]=VectorInt3(1,-1,0);
	steps[16]=VectorInt3(-1,1,0);
	steps[17]=VectorInt3(-1,-1,0); 
    }
    // //redirect cout output
    // virtual void SetUp(){
    //     originalBuffer=std::cout.rdbuf();
    //     std::cout.rdbuf(tempStream.rdbuf());
    // };

    // //restore original output
    // virtual void TearDown(){
    //     std::cout.rdbuf(originalBuffer);
    // };
    VectorInt3 steps[18];
private:
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;
    
};

TEST_F(TestMoveLocalScDiag, initialiseSetterGetter)
{
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(true);
    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.modifyMolecules().addMonomer(8,8,8);

    EXPECT_NO_THROW(ingredients.synchronize());
    EXPECT_EQ(1,ingredients.getMolecules().size());

    MoveLocalScDiag move;

    // ######################################################################## //
    // use empty init interface: dice a random monomer and a random move direction
    move.init(getIngredients());
    EXPECT_EQ(1.0,move.getProbability());
    EXPECT_EQ(0,move.getIndex());
    EXPECT_TRUE (move.getDir()*move.getDir() == 1 
                || move.getDir()*move.getDir() == 2);

    //change probability
    move.multiplyProbability(0.5);

    //add a new monomer, check if init changes properties correctly
    ingredients.modifyMolecules().addMonomer(4,8,8);
    EXPECT_NO_THROW(ingredients.synchronize());
    EXPECT_EQ(2,ingredients.getMolecules().size());

    move.init(getIngredients());
    EXPECT_EQ(1.0,move.getProbability());
    EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
    EXPECT_TRUE (move.getDir()*move.getDir() == 1 
                || move.getDir()*move.getDir() == 2);
    

    // ######################################################################## //
    // use index init interface: set the index and dice a random move direction

    EXPECT_NO_THROW(move.init(getIngredients(),0));
    EXPECT_EQ(1.0,move.getProbability());
    EXPECT_EQ(0,move.getIndex());
    EXPECT_TRUE (move.getDir()*move.getDir() == 1 
                || move.getDir()*move.getDir() == 2);
    //change probability
    move.multiplyProbability(0.5);

    EXPECT_EQ(2,ingredients.getMolecules().size());
    EXPECT_NO_THROW(move.init(getIngredients(),1));
    EXPECT_EQ(1.0,move.getProbability());
    EXPECT_EQ(1,move.getIndex());
    EXPECT_TRUE (move.getDir()*move.getDir() == 1 
                || move.getDir()*move.getDir() == 2);

    //check wrong index by initialize
    EXPECT_ANY_THROW(move.init(getIngredients(),ingredients.getMolecules().size()));
    EXPECT_ANY_THROW(move.init(getIngredients(),2)); // ( =same as molecules.size() )
    EXPECT_ANY_THROW(move.init(getIngredients(),-1));
    EXPECT_ANY_THROW(move.init(getIngredients(),68468468));

    // ######################################################################## //
    // use direction init interface: dice an index and set the move direction
    for(auto direction : steps ){
        EXPECT_NO_THROW(move.init(getIngredients(),direction));
        EXPECT_EQ(1.0,move.getProbability());
        EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
        EXPECT_EQ(direction,move.getDir());
        //change probability
        move.multiplyProbability(0.5);
    }
    VectorInt3 direction(2,0,0);
    EXPECT_ANY_THROW(move.init(getIngredients(),direction));
    direction.setAllCoordinates(0,-2,0);
    EXPECT_ANY_THROW(move.init(getIngredients(),direction));

    // ######################################################################## //
    // use direction and index init interface: set index and the move direction
    //check all possible directions
    
    for(auto index=0; index<2;index++){
        for(auto direction : steps ){
            EXPECT_NO_THROW(move.init(getIngredients(),index,direction));
            EXPECT_EQ(1.0,move.getProbability());
            EXPECT_EQ(index,move.getIndex());
            EXPECT_EQ(direction,move.getDir());
            //change probability
            move.multiplyProbability(0.5);
        }
    }
        
    //check some forbidden directions and idicees
    EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));
    EXPECT_ANY_THROW(move.init(getIngredients(),-1,direction));
    direction.setAllCoordinates(2,0,0);
    EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
    direction.setAllCoordinates(1,-1,1);
    EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
    EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));

}

TEST_F(TestMoveLocalScDiag, checkAndApply)
{
    // starting here with a new instance of ingredients!
    // class member lifetime is only one testcase!
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(false);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(true);
    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.modifyMolecules().addMonomer(2,8,8);
    ingredients.modifyMolecules().addMonomer(0,8,8);

    EXPECT_NO_THROW(ingredients.synchronize());
    EXPECT_EQ(2,ingredients.getMolecules().size());

    MoveLocalScDiag move;

    // ################################# //
    //check move with FeatureExcludedVolumeSc:
    VectorInt3 direction(1,0,0);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_FALSE(move.check(ingredients));

    direction.setAllCoordinates(1,1,0);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_FALSE(move.check(ingredients));

    direction.setAllCoordinates(1,0,1);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_FALSE(move.check(ingredients));
    
    direction.setAllCoordinates(1,-1,0);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_FALSE(move.check(ingredients));

    direction.setAllCoordinates(1,1,0);
    EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);

    direction.setAllCoordinates(0,-1,1);
    EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);

    direction.setAllCoordinates(0,1,-1);
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);

    direction.setAllCoordinates(0,0,-1);
    EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
    EXPECT_TRUE(move.check(ingredients));

    direction.setAllCoordinates(1,0,1);
    EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    // now monomers are at positions 

  // ################################# //
  //check move with FeatureMoleculesIO:
  //bondset:
  ingredients.modifyMolecules().connect(0,1);
  direction.setAllCoordinates(1,-1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  EXPECT_EQ(VectorInt3(3,8,9),VectorInt3(ingredients.getMolecules()[0]));
  EXPECT_EQ(VectorInt3(1,9,8),VectorInt3(ingredients.getMolecules()[1]));
  // now monomers are at positions 2,9,9 and 0,9,9

  //nonperiodic wall
  direction.setAllCoordinates(-1,0,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  direction.setAllCoordinates(-1,0,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_FALSE(move.check(ingredients));
}
TEST_F(TestMoveLocalScDiag, unentangledRings)
{
    // starting here with a new instance of ingredients!
    // class member lifetime is only one testcase!
    ingredients.setBoxX(128);
    ingredients.setBoxY(32);
    ingredients.setBoxZ(32);
    ingredients.setPeriodicX(false);
    ingredients.setPeriodicY(false);
    ingredients.setPeriodicZ(false);
    ingredients.modifyBondset().addBFMclassicBondset();
    
    ingredients.modifyMolecules().addMonomer(4,8,8);
    ingredients.modifyMolecules().addMonomer(7,9,8);
    ingredients.modifyMolecules().addMonomer(6,12,8);
    ingredients.modifyMolecules().addMonomer(3,11,8);
    ingredients.modifyMolecules().connect(0,1);
    ingredients.modifyMolecules().connect(1,2);
    ingredients.modifyMolecules().connect(2,3);
    ingredients.modifyMolecules().connect(3,0);

    ingredients.modifyMolecules().addMonomer(5,10,8);
    ingredients.modifyMolecules().addMonomer(3,9,10);
    ingredients.modifyMolecules().addMonomer(1,10,8);
    ingredients.modifyMolecules().addMonomer(3,9,6);
    ingredients.modifyMolecules().connect(4,5);
    ingredients.modifyMolecules().connect(5,6);
    ingredients.modifyMolecules().connect(6,7);
    ingredients.modifyMolecules().connect(7,4);
    ingredients.synchronize();
    MoveLocalScDiag move; 
    // std::string filename("TestMoveLocalScDiag.bfm");
    // AnalyzerWriteBfmFile<IngredientsType> BfmWriter(filename, ingredients, AnalyzerWriteBfmFile<IngredientsType>::NEWFILE);
    // BfmWriter.initialize();
    for (uint32_t j=0; j < 100; j++){
        for(uint32_t n=0;n<100;n++){
            for(size_t m=0;m<ingredients.getMolecules().size();m++) {
                move.init(ingredients);
                if( move.getDir().getX() == -1 && move.getIndex() >=4 ) //prefer positive x direction
                {
                    move.multiplyProbability(0.2);
                }
                else if( move.getDir().getX() == 1 && move.getIndex() <4 ) //prefer positive x direction
                {
                    move.multiplyProbability(0.2);
                }
                if(move.check(ingredients)==true)
                    move.apply(ingredients);
            }
            ingredients.modifyMolecules().setAge(ingredients.getMolecules().getAge()+100);
        }
        ingredients.synchronize();
        // BfmWriter.execute();
    }
    // BfmWriter.cleanup();
    VectorDouble3 COM1(0.,0.,0.),COM2(0.,0.,0.);
    for (auto i=0; i < 4 ; i++)
        COM1+=ingredients.getMolecules()[i].getVector3D();
    for (auto i=4; i < 8 ; i++)
        COM2+=ingredients.getMolecules()[i].getVector3D();
    COM1.setX(COM1.getX()/4.);
    COM2.setX(COM2.getX()/4.);
    EXPECT_GE(COM2.getX()-COM1.getX(),64);
    // remove(filename.c_str());
}
