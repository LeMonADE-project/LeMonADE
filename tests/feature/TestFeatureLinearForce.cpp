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
#include <cstdio>
#include <sstream>
#include <string>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureLinearForce.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>

class TestFeatureLinearForce : public ::testing::Test{
public:
    typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureLinearForce ) Features;
    typedef ConfigureSystem<VectorInt3,Features,17> Config;
    typedef Ingredients<Config> Ing;
    Ing ingredients;

    //redirect cout output
    // virtual void SetUp(){
    //     originalBuffer=std::cout.rdbuf();
    //     std::cout.rdbuf(tempStream.rdbuf());
    // };

    // //restore original output
    // virtual void TearDown(){
    //     std::cout.rdbuf(originalBuffer);
    // };
    //tension force
    const double force=1.212;
    //! filename for output 
    const  std::string  filename="TestFeatureLinearForce.bfm";
    //set some basic thins 
    void initIng(){
        //prepare ingredients
        ingredients.setBoxX(12);
        ingredients.setBoxY(12);
        ingredients.setBoxZ(12);
        ingredients.setPeriodicX(0);
        ingredients.setPeriodicY(0);
        ingredients.setPeriodicZ(0);
        ingredients.synchronize();
    }
private:
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;

};
TEST_F(TestFeatureLinearForce,SetterGetter){
    ingredients.setAmplitudeForce(force);
    EXPECT_EQ(ingredients.getAmplitudeForce(),force);
    EXPECT_FALSE(ingredients.isForceOn());
    ingredients.setForceOn(true);
    EXPECT_TRUE(ingredients.isForceOn());
}
TEST_F(TestFeatureLinearForce,WriterReader){
    initIng();
    ingredients.setAmplitudeForce(force);
    ingredients.setForceOn(true);
    ingredients.modifyMolecules().resize(1);
    ingredients.modifyMolecules()[0].setAllCoordinates(6,6,6);
    ingredients.synchronize();
    AnalyzerWriteBfmFile<Ing> writer(filename, ingredients, AnalyzerWriteBfmFile<Ing>::NEWFILE);
    writer.initialize();
    writer.execute();

    Ing ingredients2;
    UpdaterReadBfmFile<Ing> reader(filename, ingredients2, UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE);
    reader.initialize();
    EXPECT_EQ(ingredients2.getAmplitudeForce(),force);
    EXPECT_TRUE(ingredients2.isForceOn());
    EXPECT_EQ(0,remove(filename.c_str()));

}
TEST_F(TestFeatureLinearForce,checkMoveLocalSc){
    ingredients.setAmplitudeForce(2);
    ingredients.setForceOn(true);
    initIng();
    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(6,6,6);
    ingredients.modifyMolecules()[1].setAllCoordinates(6,6,6);
    ingredients.modifyMolecules()[2].setAllCoordinates(6,6,6);
    ingredients.modifyMolecules()[0].setAttributeTag(1);
    ingredients.modifyMolecules()[1].setAttributeTag(4);
    ingredients.modifyMolecules()[2].setAttributeTag(5);
    ingredients.synchronize();
    MoveLocalSc move; 
    //monomer without tag
    move.init(ingredients,0,VectorInt3(1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,0,VectorInt3(-1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);
    

    move.init(ingredients,0,VectorInt3(0,1,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);
    
    move.init(ingredients,1,VectorInt3(0,-1,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,0,VectorInt3(0,0,1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,0,VectorInt3(0,0,-1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);
    
    //direction where no force is applied on for a taged monomer 
    //second monomer 
    move.init(ingredients,1,VectorInt3(0,1,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,1,VectorInt3(0,-1,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,1,VectorInt3(0,0,1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,1,VectorInt3(0,0,-1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    //third monomer 
    move.init(ingredients,2,VectorInt3(0,1,0));    
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,2,VectorInt3(0,1,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,2,VectorInt3(0,0,1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,2,VectorInt3(0,0,-1)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    //direction where force is applied on for a taged monomer 

    move.init(ingredients,1,VectorInt3(-1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);

    move.init(ingredients,1,VectorInt3( 1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , exp(-2));
    
    move.init(ingredients,2,VectorInt3(1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , 1.0);
    
    move.init(ingredients,2,VectorInt3(-1,0,0)); 
    EXPECT_TRUE(move.check(ingredients));
    EXPECT_EQ(move.getProbability() , exp(-2));
    
}



