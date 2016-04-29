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

#include <LeMonADE/LeMonADE.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>

using namespace std;

TEST(TestFeatureExcludedVolumeSc,CheckMoves)
{
 
  typedef LOKI_TYPELIST_2(FeatureBondset< >,FeatureExcludedVolumeSc< >) Features;

  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  Ing ingredients;
  
  //prepare ingredients
    ingredients.setBoxX(32);
    ingredients.setBoxY(32);
    ingredients.setBoxZ(32);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    
    //one move of every type
    MoveBase basemove;
    MoveLocalSc scmove;
    MoveAddScMonomer addmove;
    
    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);
    
    ingredients.synchronize(ingredients);
    
    // **************   check base move   **************
    basemove.init(ingredients);
    EXPECT_TRUE(basemove.check(ingredients));
    //should change nothing
    basemove.apply(ingredients);
    EXPECT_EQ(ingredients.getMolecules()[0],VectorInt3(0,0,0));
    EXPECT_EQ(ingredients.getMolecules()[1],VectorInt3(2,0,0));
    EXPECT_EQ(ingredients.getMolecules()[2],VectorInt3(1,0,30));
    
    // **************   check sc move   **************
    //collossion to the right
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision 0 downwards via periodic boundaries
    while((scmove.getDir().getZ()!=-1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision to the left
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision 1 downwards via periodic boundaries
    while((scmove.getDir().getZ()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision upwards via periodic boundaries
    while((scmove.getDir().getZ()!=1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //some possible moves
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    while((scmove.getDir().getY()!=1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    while((scmove.getDir().getY()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    
    //shift monomer and try again
    ingredients.modifyMolecules()[1].setAllCoordinates(2,1,0);
    ingredients.synchronize(ingredients);
    
    //collossion to the right
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision 0 downwards via periodic boundaries
    while((scmove.getDir().getZ()!=-1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision to the left
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision 1 downwards via periodic boundaries
    while((scmove.getDir().getZ()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision upwards via periodic boundaries
    while((scmove.getDir().getZ()!=1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //use a localmove (without synchronize) and try again
    while((scmove.getDir().getY()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_TRUE(scmove.check(ingredients));
    scmove.apply(ingredients);
    
    //collossion to the right
    while((scmove.getDir().getX()!=1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision 0 downwards via periodic boundaries
    while((scmove.getDir().getZ()!=-1) || (scmove.getIndex()!=0)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision to the left
    while((scmove.getDir().getX()!=-1) || (scmove.getIndex()!=1)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    //collision upwards via periodic boundaries
    while((scmove.getDir().getZ()!=1) || (scmove.getIndex()!=2)) scmove.init(ingredients);
    EXPECT_FALSE(scmove.check(ingredients));
    
    // **************   check add move   **************
    
    addmove.init(ingredients);
    
    addmove.setPosition(0,0,0);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,0,0);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(2,0,0);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(3,0,0);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(4,0,0);
    EXPECT_TRUE(addmove.check(ingredients));
    
    addmove.setPosition(0,0,1);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(0,0,2);
    EXPECT_TRUE(addmove.check(ingredients));
    
    addmove.setPosition(1,0,-1);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,-1,-1);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,31,31);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,-1,0);
    EXPECT_TRUE(addmove.check(ingredients));
    
    //apply move and check position again
    addmove.apply(ingredients);
    addmove.init(ingredients);
    EXPECT_FALSE(addmove.check(ingredients));
    
}

TEST(TestFeatureExcludedVolumeSc,ApplyMoves)
{
 
  typedef LOKI_TYPELIST_2(FeatureBondset< >,FeatureExcludedVolumeSc< >) Features;

  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  Ing ingredients;
  
  //prepare ingredients
  ingredients.setBoxX(32);
  ingredients.setBoxY(32);
  ingredients.setBoxZ(32);
  ingredients.setPeriodicX(1);
  ingredients.setPeriodicY(1);
  ingredients.setPeriodicZ(1);
  
  ingredients.modifyMolecules().resize(1);
  ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
  
  ingredients.synchronize(ingredients);
  
  MoveBase basemove;
  basemove.init(ingredients);
  basemove.apply(ingredients);
  EXPECT_EQ(ingredients.getMolecules()[0],VectorInt3(0,0,0));
    
}
    
    
TEST(TestFeatureExcludedVolumeSc,CheckInterface)
{
	typedef LOKI_TYPELIST_1(FeatureExcludedVolumeSc< >) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing ingredients;

	 //prepare ingredients1
	    ingredients.setBoxX(32);
	    ingredients.setBoxY(32);
	    ingredients.setBoxZ(32);
	    ingredients.setPeriodicX(1);
	    ingredients.setPeriodicY(1);
	    ingredients.setPeriodicZ(1);
	    //set up ingredients
	    typename Ing::molecules_type& molecules=ingredients.modifyMolecules();

	    molecules.resize(3);
	    molecules[0].setAllCoordinates(9,10,10);
	    molecules[1].setAllCoordinates(13,10,10);

	    //initially is set to false
	    //before synchronize: latticeIsNotUpdated
	    EXPECT_FALSE(ingredients.isLatticeFilledUp());

	    ingredients.synchronize(ingredients);

	    //after synchronize: latticeIsUpdated
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

	    //manually setting it to false
	    ingredients.setLatticeFilledUp(false);
	    EXPECT_FALSE(ingredients.isLatticeFilledUp());

	    //manually setting it to true
	    ingredients.setLatticeFilledUp(true);
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

	    //before synchronize: latticeIsNotUpdated
	    ingredients.synchronize(ingredients);
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

}