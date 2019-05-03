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

#include "gtest/gtest.h"

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>


using namespace std;

class TestFeatureExcludedVolumeSc: public ::testing::Test{
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

TEST_F(TestFeatureExcludedVolumeSc,Moves)
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
    UnknownMove basemove;
    MoveLocalSc scmove;
    MoveAddMonomerSc<> addmove;

    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);

    ingredients.synchronize(ingredients);

    // **************   check inconsistent moves   *****
    MoveLocalBcc bccmove;
    bccmove.init(ingredients);
    EXPECT_ANY_THROW(bccmove.check(ingredients));

    MoveAddMonomerBcc<> addbccmove;
    addbccmove.init(ingredients);
    EXPECT_ANY_THROW(addbccmove.check(ingredients));


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

    //create random system and check it by synchronize
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

}

TEST_F(TestFeatureExcludedVolumeSc,DiagonalMoves)
{

  typedef LOKI_TYPELIST_2(FeatureBondset< >,FeatureExcludedVolumeSc< >) Features;

  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  Ing ingredients;

  //prepare ingredients
    ingredients.setBoxX(8);
    ingredients.setBoxY(8);
    ingredients.setBoxZ(8);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);

    //one move of every type
    MoveLocalScDiag DiagMove;

    ingredients.modifyMolecules().resize(3);
    VectorInt3 InitPos(0,0,0);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,2,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);

    ingredients.synchronize(ingredients);


    // **************   check diagonal move   **************
    EXPECT_NO_THROW(DiagMove.init(ingredients));
    EXPECT_NO_THROW(DiagMove.init(ingredients,0));
    EXPECT_NO_THROW(DiagMove.init(ingredients,0,VectorInt3(1,1,0)));
    EXPECT_ANY_THROW(DiagMove.init(ingredients,0,VectorInt3(1,1,1)));
    
    EXPECT_NO_THROW(DiagMove.init(ingredients,0,VectorInt3(1,1,0)));
    EXPECT_FALSE(DiagMove.check(ingredients));

    // setup all possible move directions (3*3*3)=18
    VectorInt3 steps[18];
    //ordinary moves
    steps[0]=VectorInt3(1,0,0);
	  steps[1]=VectorInt3(-1,0,0);
	  steps[2]=VectorInt3(0,1,0);
	  steps[3]=VectorInt3(0,-1,0);
	  steps[4]=VectorInt3(0,0,1);
	  steps[5]=VectorInt3(0,0,-1); //fail -> monomer 2
	  //diagonal moves
	  steps[6]=VectorInt3(0,1,1);
	  steps[7]=VectorInt3(0,-1,1);
	  steps[8]=VectorInt3(0,1,-1);//fail -> monomer 2
	  steps[9]=VectorInt3(0,-1,-1);//fail -> monomer 2
	  steps[10]=VectorInt3(1,0,1);
	  steps[11]=VectorInt3(1,0,-1);//fail -> monomer 2
	  steps[12]=VectorInt3(-1,0,1);
	  steps[13]=VectorInt3(-1,0,-1);
	  steps[14]=VectorInt3(1,1,0); //fail -> monomer 1
	  steps[15]=VectorInt3(1,-1,0);
	  steps[16]=VectorInt3(-1,1,0);
	  steps[17]=VectorInt3(-1,-1,0);

    for (size_t i =0  ;i < 18 ; i++ )
    {
      std::cout << "Use "<< steps[i] <<std::endl;
      if ( !(  i==5
	    || i==8
	    || i==9
	    || i==11
	    || i==14) )
      {
	      EXPECT_NO_THROW(DiagMove.init(ingredients,0,steps[i]));
	      EXPECT_TRUE(DiagMove.check(ingredients));
	      EXPECT_NO_THROW(DiagMove.apply(ingredients));
	      EXPECT_NO_THROW(ingredients.synchronize());
	      ingredients.modifyMolecules()[0].modifyVector3D()=InitPos;
	      EXPECT_NO_THROW(ingredients.synchronize());
	      EXPECT_EQ(InitPos, ingredients.getMolecules()[0].getVector3D() );
      }else{
	      EXPECT_NO_THROW(DiagMove.init(ingredients,0,steps[i]));
	      EXPECT_FALSE(DiagMove.check(ingredients));
	      EXPECT_NO_THROW(DiagMove.apply(ingredients));
	      EXPECT_ANY_THROW(ingredients.synchronize());
	      ingredients.modifyMolecules()[0].modifyVector3D()=InitPos;
	      EXPECT_NO_THROW(ingredients.synchronize());
	      EXPECT_EQ(InitPos, ingredients.getMolecules()[0].getVector3D() );
      }
    }
    // apply the move (0,1,1)
    EXPECT_NO_THROW(DiagMove.init(ingredients,0,steps[6]));
    EXPECT_TRUE(DiagMove.check(ingredients));
    EXPECT_NO_THROW(DiagMove.apply(ingredients));
    
    // now move step[9] (0,-1,-1) should be accepted
    EXPECT_NO_THROW(DiagMove.init(ingredients,0,steps[9]));
    EXPECT_TRUE(DiagMove.check(ingredients));
    EXPECT_NO_THROW(DiagMove.apply(ingredients));
    
}


TEST_F(TestFeatureExcludedVolumeSc,CheckInterface)
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
