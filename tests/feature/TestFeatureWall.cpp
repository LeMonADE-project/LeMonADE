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

#include <iostream>
#include <cstdio>

#include "gtest/gtest.h"

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureWall.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>


class TestFeatureWall: public ::testing::Test{
public:
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

TEST(TestFeatureWall,Moves)
{

    typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureWall) Features;

    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients;

    //prepare ingredients
    ingredients.setBoxX(32);
    ingredients.setBoxY(32);
    ingredients.setBoxZ(32);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);

    //one move of every type
    MoveLocalSc scmove;
    MoveAddMonomerSc<> addmove;

    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(2,2,0);

    EXPECT_NO_THROW(ingredients.synchronize(ingredients));

    //************  scmove  ***************

    //build wall in y-z-plane at x=2 -> x=1 should be forbidden
    Wall wall1;
    wall1.setBase(2,0,0);
    // check forbidden nomal vector: must be a unit vector! P(1,0,0)
    EXPECT_ANY_THROW(wall1.setNormal(2,0,0));
    EXPECT_ANY_THROW(wall1.setNormal(0,1,1));

    //add correct normal vector
    wall1.setNormal(1,0,0);
    ingredients.addWall(wall1);
/*
    std::cout << "getWalls(): " << ingredients.getWalls()[0].getBase().getX() << " " << ingredients.getWalls()[0].getBase().getY() << " " << ingredients.getWalls()[0].getBase().getZ() << " " << ingredients.getWalls()[0].getNormal().getX() << " " << ingredients.getWalls()[0].getNormal().getY() << " " << ingredients.getWalls()[0].getNormal().getZ() << std::endl;
*/

    //check monomer 0 movement to backside of wall
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients); //init(ingredients) searches randomly one monomer in ingredients and one direction for the move -> while you have not the right monomer and the prefered direction, keep on searching randomly
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check monomer 1 movement to frontside of wall
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check monomer 2 movement to frontside of wall
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 0 movements
    while((scmove.getDir().getX()==1) || (scmove.getIndex()!=0)) scmove.init(ingredients); //search randomly for any direction except direction to the wall (x=1) for monomer 0
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 1 movements
    while((scmove.getDir().getX()==-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 2 movements
    while((scmove.getDir().getX()==-1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //build additional wall in x-z-plane at y=2 -> x=1 and y=1 should be forbidden now
    Wall wall2;
    wall2.setBase(0,2,0);
    wall2.setNormal(0,1,0);
    ingredients.addWall(wall2);
/*
    for(size_t i=0; i< ingredients.getWalls().size(); i++){
        std::cout << "getWalls()[" << i << "]: " << ingredients.getWalls()[i].getBase().getX() << " " << ingredients.getWalls()[i].getBase().getY() << " " << ingredients.getWalls()[i].getBase().getZ() << " " << ingredients.getWalls()[i].getNormal().getX() << " " << ingredients.getWalls()[i].getNormal().getY() << " " << ingredients.getWalls()[i].getNormal().getZ() << std::endl;
    }
*/
    //check monomer 0 movement to backside of both walls
    while((scmove.getDir().getY()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check monomer 1 movement to frontside of wall1 and backside of wall2
    while((scmove.getDir().getY()!=1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check monomer 2 movement to frontside of both walls
    while((scmove.getDir().getY()!=-1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 0 movements
    while((scmove.getDir().getX()==1) || (scmove.getDir().getY()==1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 1 movements
    while((scmove.getDir().getX()==-1) || (scmove.getDir().getY()==1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 2 movements
    while((scmove.getDir().getX()==-1) || (scmove.getDir().getY()==-1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    ingredients.clearAllWalls();

    //add former wall1 again
    ingredients.addWall(wall1);

    //build second wall in y-z-plane at x=4 -> x=3 should be forbidden
    Wall wall3;
    wall3.setBase(4,0,0);
    wall3.setNormal(1,0,0);
    ingredients.addWall(wall3);

    //check monomer 1 movement to backside of wall3 (wall1 was already checked above)
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check monomer 2 movement to backside of wall3
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 1 movements
    while((scmove.getDir().getX()==-1) || (scmove.getDir().getX()==1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));

    //check possible monomer 2 movements
    while((scmove.getDir().getX()==-1) || (scmove.getDir().getX()==1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_TRUE(ingredients.checkMove(ingredients,scmove));


    //************  addmove  ***************

    ingredients.clearAllWalls();

    //add former walls again
    ingredients.addWall(wall1);
    ingredients.addWall(wall2);
    ingredients.addWall(wall3);

    //create random system and check it by synchronize
    std::cout << "create random config by addmove and check by snychronize"<<std::endl;
    for(uint32_t i=0;i<2000;i++){
        bool accept_move(false);
        while(!accept_move){
            addmove.setPosition(std::rand()%32,std::rand()%32,std::rand()%32);
            if(addmove.check(ingredients)){
                addmove.apply(ingredients);
                accept_move=true;
            }
        }
    }
    EXPECT_NO_THROW(ingredients.synchronize());


    //check monomer added at wall1
    addmove.setPosition(1,0,0);
    EXPECT_FALSE(ingredients.checkMove(ingredients,addmove));

    //check possible monomer positions
    addmove.setPosition(2,0,0);
    EXPECT_TRUE(ingredients.checkMove(ingredients,addmove));

    //check monomer added at wall2
    addmove.setPosition(2,1,0);
    EXPECT_FALSE(ingredients.checkMove(ingredients,addmove));

    //check possible monomer positions
    addmove.setPosition(2,2,0);
    EXPECT_TRUE(ingredients.checkMove(ingredients,addmove));

    //check monomer added at wall3
    addmove.setPosition(3,2,0);
    EXPECT_FALSE(ingredients.checkMove(ingredients,addmove));

    //check possible monomer positions
    addmove.setPosition(4,2,0);
    EXPECT_TRUE(ingredients.checkMove(ingredients,addmove));

}


TEST(TestFeatureWall,ReadWriteRoutine)
{

    typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureWall) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients1; //to write bfm-file

    typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureWall) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients2; //to read bfm-file

    //prepare ingredients1 to write bfm-file
    ingredients1.setBoxX(32);
    ingredients1.setBoxY(32);
    ingredients1.setBoxZ(32);
    ingredients1.setPeriodicX(1);
    ingredients1.setPeriodicY(1);
    ingredients1.setPeriodicZ(1);

    Wall wall4;
    wall4.setBase(5,0,0);
    wall4.setNormal(1,0,0);
    ingredients1.addWall(wall4);

    Wall wall5;
    wall5.setBase(0,4,0);
    wall5.setNormal(0,1,0);
    ingredients1.addWall(wall5);

    EXPECT_NO_THROW(ingredients1.synchronize());

    TaskManager taskmanager;
    taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("test_wall.bfm",ingredients1,AnalyzerWriteBfmFile<IngredientsType>::NEWFILE));

    taskmanager.initialize();
    taskmanager.run();
    taskmanager.cleanup();

    TaskManager taskmanager1;
    taskmanager1.addUpdater(new UpdaterReadBfmFile<IngredientsType>("test_wall.bfm",ingredients2,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE));

    taskmanager1.initialize();
    taskmanager1.run();
    taskmanager1.cleanup();

    remove("test_wall.bfm");

    EXPECT_TRUE(ingredients1.getWalls()[0].getBase().getX()==ingredients2.getWalls()[0].getBase().getX());
    EXPECT_TRUE(ingredients1.getWalls()[1].getBase().getX()==ingredients2.getWalls()[1].getBase().getX());

    EXPECT_NO_THROW(ingredients1.synchronize());
    EXPECT_NO_THROW(ingredients2.synchronize());

}
