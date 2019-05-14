/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2016 by
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
//#include <cstdlib>
#include "gtest/gtest.h"

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureExcludedVolumeBcc.h>

using namespace std;

class TestFeatureExcludedVolumeBcc: public ::testing::Test{
public:
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

TEST_F(TestFeatureExcludedVolumeBcc,Moves)
{

  typedef LOKI_TYPELIST_2(FeatureBondset< >,FeatureExcludedVolumeBcc< >) Features;

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
    UnknownMove basemove;
    MoveLocalBcc bccmove;
    MoveAddMonomerBcc<> addmove;

    ingredients.modifyMolecules().resize(4);
    ingredients.modifyMolecules()[0].setAllCoordinates(1,1,1);
    ingredients.modifyMolecules()[1].setAllCoordinates(1,5,1);
    ingredients.modifyMolecules()[2].setAllCoordinates(31,3,1);
    ingredients.modifyMolecules()[3].setAllCoordinates(1,1,3);

    EXPECT_NO_THROW(ingredients.synchronize(ingredients));

    // **************   check inconsistent moves   *****
    MoveLocalSc scmove;
    scmove.init(ingredients);
    EXPECT_ANY_THROW(scmove.check(ingredients));

    MoveAddMonomerSc<> addscmove;
    addscmove.init(ingredients);
    EXPECT_ANY_THROW(addscmove.check(ingredients));

    // **************   check base move   **************
    basemove.init(ingredients);
    EXPECT_TRUE(basemove.check(ingredients));
    //should change nothing
    basemove.apply(ingredients);
    EXPECT_EQ(ingredients.getMolecules()[0],VectorInt3(1,1,1));
    EXPECT_EQ(ingredients.getMolecules()[1],VectorInt3(1,5,1));
    EXPECT_EQ(ingredients.getMolecules()[2],VectorInt3(31,3,1));
    EXPECT_EQ(ingredients.getMolecules()[3],VectorInt3(1,1,3));

    // **************   check bcc move   **************
    //collossion monomer 2 x+1 to 0 and 1
    while((bccmove.getDir().getX()!=1) || (bccmove.getIndex()!=2)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    //collision monomer 0 y+1, x+1 to 2
    while((bccmove.getDir().getX()!=-1) || (bccmove.getDir().getY()!=1) || (bccmove.getIndex()!=0)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    //collision monomer 0 z+1 to 3
    while((bccmove.getDir().getZ()!=1) || (bccmove.getIndex()!=0)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    //collision monomer 1 x-1, y-1 to 2
    while((bccmove.getDir().getX()!=-1) || (bccmove.getDir().getY()!=-1) || (bccmove.getIndex()!=1)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    //collision monomer 3 z-1 to 0
    while((bccmove.getDir().getZ()!=-1) || (bccmove.getIndex()!=3)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    /* ################################ */
    //some possible moves
    while((bccmove.getDir().getZ()!=-1) || (bccmove.getDir().getX()!=1) || (bccmove.getIndex()!=0)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    while((bccmove.getDir().getX()!=-1) || (bccmove.getIndex()!=2)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    while((bccmove.getDir().getZ()!=1) || (bccmove.getIndex()!=3)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    while((bccmove.getDir().getX()!=1) || (bccmove.getIndex()!=1)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    //shift monomer and try again
    ingredients.modifyMolecules()[1].setAllCoordinates(4,2,2);
    EXPECT_NO_THROW(ingredients.synchronize(ingredients));

    //some possible moves
    while((bccmove.getDir().getX()!=-1) || (bccmove.getIndex()!=1)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    //some possible moves
    while((bccmove.getDir().getX()!=1) || (bccmove.getDir().getZ()!=1) || (bccmove.getIndex()!=3)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));

    //use a localmove (without synchronize) and try again
    while((bccmove.getDir().getX()!=-1) || (bccmove.getDir().getY()!=1) || (bccmove.getDir().getZ()!=-1) || (bccmove.getIndex()!=1)) bccmove.init(ingredients);
    EXPECT_TRUE(bccmove.check(ingredients));
    bccmove.apply(ingredients);

    //collision monomer 2 up to 3
    while((bccmove.getDir().getX()!=1) || (bccmove.getDir().getY()!=1) || (bccmove.getIndex()!=0)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    //collision monomer 2 down-left-z-down to 1
    while((bccmove.getDir().getX()!=-1) || (bccmove.getDir().getY()!=-1) ||  (bccmove.getIndex()!=1)) bccmove.init(ingredients);
    EXPECT_FALSE(bccmove.check(ingredients));

    // **************   check add move   **************
    ingredients.modifyMolecules().resize(2);
    ingredients.modifyMolecules()[0].setAllCoordinates(1,1,1);
    ingredients.modifyMolecules()[1].setAllCoordinates(5,1,1);
    EXPECT_NO_THROW(ingredients.synchronize(ingredients));

    addmove.init(ingredients);

    //position itsself occupied
    addmove.setPosition(1,1,1);
    EXPECT_FALSE(addmove.check(ingredients));

    //position itsself occupied
    addmove.setPosition(5,1,1);
    EXPECT_FALSE(addmove.check(ingredients));

    //position on "o" lattice occupied
    addmove.setPosition(0,0,0);
    EXPECT_FALSE(addmove.check(ingredients));

    //position on "o" lattice occupied
    addmove.setPosition(2,2,2);
    EXPECT_FALSE(addmove.check(ingredients));

    //position on "o" lattice occupied
    addmove.setPosition(4,0,0);
    EXPECT_FALSE(addmove.check(ingredients));

    //position on "o" lattice occupied
    addmove.setPosition(4,2,2);
    EXPECT_FALSE(addmove.check(ingredients));

    //position on "o" lattice occupied
    addmove.setPosition(6,2,2);
    EXPECT_FALSE(addmove.check(ingredients));

    //one position is even, the others odd
    addmove.setPosition(9,8,8);
    EXPECT_FALSE(addmove.check(ingredients));

    //one position is even, the others odd
    addmove.setPosition(8,8,9);
    EXPECT_FALSE(addmove.check(ingredients));

    //one position is even, the others odd
    addmove.setPosition(8,9,8);
    EXPECT_FALSE(addmove.check(ingredients));

    //creation of "lines" of monomer
    addmove.setPosition(3,1,1);
    EXPECT_FALSE(addmove.check(ingredients));

    //possible creations
    addmove.setPosition(3,-1,-1);
    EXPECT_TRUE(addmove.check(ingredients));  //false?

    //possible creations
    addmove.setPosition(3,31,31);
    EXPECT_TRUE(addmove.check(ingredients));

    //possible creations
    addmove.setPosition(7,1,1);
    EXPECT_TRUE(addmove.check(ingredients));

    //possible creations
    addmove.setPosition(-1,1,1);
    EXPECT_TRUE(addmove.check(ingredients));  //false?

    //possible creations
    addmove.setPosition(31,1,1);
    EXPECT_TRUE(addmove.check(ingredients));

    //apply move and check position again
    addmove.apply(ingredients);
    addmove.init(ingredients);
    EXPECT_FALSE(addmove.check(ingredients));

    //create random system and check it by synchronize
    for(uint32_t i=0;i<800;i++){
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

}

TEST_F(TestFeatureExcludedVolumeBcc,CheckInterface)
{
	typedef LOKI_TYPELIST_1(FeatureExcludedVolumeBcc< >) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> IngredientsType;
	IngredientsType ingredients;

	 //prepare ingredients1
	    ingredients.setBoxX(32);
	    ingredients.setBoxY(32);
	    ingredients.setBoxZ(32);
	    ingredients.setPeriodicX(1);
	    ingredients.setPeriodicY(1);
	    ingredients.setPeriodicZ(1);
	    //set up ingredients
	    typename IngredientsType::molecules_type& molecules=ingredients.modifyMolecules();

	    molecules.resize(3);
	    molecules[0].setAllCoordinates(10,10,10);
	    molecules[1].setAllCoordinates(14,10,10);

	    //initially is set to false
	    //before synchronize: latticeIsNotUpdated
	    EXPECT_FALSE(ingredients.isLatticeFilledUp());

	    EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	    //after synchronize: latticeIsUpdated
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

	    //manually setting it to false
	    ingredients.setLatticeFilledUp(false);
	    EXPECT_FALSE(ingredients.isLatticeFilledUp());

	    //manually setting it to true
	    ingredients.setLatticeFilledUp(true);
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

	    //before synchronize: latticeIsNotUpdated
	    EXPECT_NO_THROW(ingredients.synchronize(ingredients));
	    EXPECT_TRUE(ingredients.isLatticeFilledUp());

	    //set inconsistent monomer coordinates (not all even or all odd)
	    molecules[2].setAllCoordinates(24,1,1);
	    EXPECT_ANY_THROW(ingredients.synchronize(ingredients));

}
