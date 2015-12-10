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
#include <LeMonADE/feature/FeatureExcludedVolume.h>

using namespace std;

TEST(ExcludedVolumeTest,CheckScMove)
{
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureExcludedVolume<>) Features;

    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> Ing;
    Ing myIngredients;


    //prepare myIngredients1
    myIngredients.setBoxX(64);
    myIngredients.setBoxY(64);
    myIngredients.setBoxZ(64);
    myIngredients.setPeriodicX(1);
    myIngredients.setPeriodicY(1);
    myIngredients.setPeriodicZ(1);
    //set up ingrediens
    typename Ing::molecules_type& molecules=myIngredients.modifyMolecules();


    molecules.resize(3);
    molecules[0].setAllCoordinates(9,10,10);
    molecules[2].setAllCoordinates(13,10,10);

    //now check excluded volume for moves of monomer 0.
    //for this the testmove is used. at the same time monomer 1 is moved around
    //monomer 0 using update move, to test all possible situations
    MoveLocalSc testmove,updateMove;

    //initialize testmove with monomer index 1 and positive x-direction
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)) testmove.init(myIngredients);
    std::cout<<"testing move in positive x direction..."<<testmove.getDir()<<" monomer index "<<testmove.getIndex()<<std::endl;

    molecules[1].setAllCoordinates(11,9,9);
    myIngredients.synchronize(myIngredients);

    //now testmove should return false when checked
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //initialize update move to move monomer 1 in y direction
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);
    //check 11 10 9
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 11 9
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //initialize update move to move monomer 1 in z direction
    while(updateMove.getDir().getZ()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 11 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //move monomer 1 in negative y direction
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 10 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 9 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //move monomer 1 in z direction
    while(updateMove.getDir().getZ()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 9 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //move monomer 1 in y direction
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 10 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 11 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //now check moves in negative x direction. this is done in analogy to
    //the check in positive x direction. monomer 1 is moved around. this
    //time the testmove is initialized for monomer 2 (instead of monomer0)

    //initialize testmove with monomer 2 and negative x direction
    while((testmove.getDir().getX()!=-1) || (testmove.getIndex()!=2)) testmove.init(myIngredients);
    std::cout<<"testing move in negative x direction..."<<testmove.getDir()<<" monomer index "<<testmove.getIndex()<<std::endl;

    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);
    //check 11 10 9
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 11 9
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    while(updateMove.getDir().getZ()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 11 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 10 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 9 10
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    while(updateMove.getDir().getZ()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 9 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients);

    //check 11 10 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;

    //check 11 11 11
    updateMove.apply(myIngredients);
    EXPECT_FALSE(testmove.check(myIngredients));
    cout<<"pos "<<molecules[1]<<" "<<testmove.check(myIngredients)<<std::endl;


}

TEST(ExcludedVolumeTest,CheckInterface)
{
	typedef LOKI_TYPELIST_1(FeatureExcludedVolume<>) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	 //prepare myIngredients1
	    myIngredients.setBoxX(64);
	    myIngredients.setBoxY(64);
	    myIngredients.setBoxZ(64);
	    myIngredients.setPeriodicX(1);
	    myIngredients.setPeriodicY(1);
	    myIngredients.setPeriodicZ(1);
	    //set up ingrediens
	    typename Ing::molecules_type& molecules=myIngredients.modifyMolecules();


	    molecules.resize(3);
	    molecules[0].setAllCoordinates(9,10,10);
	    molecules[1].setAllCoordinates(13,10,10);

	    //initially is set to false
	    //before synchronize: latticeIsNotUpdated
	    EXPECT_FALSE(myIngredients.isLatticeFilledUp());

	    myIngredients.synchronize(myIngredients);

	    //after synchronize: latticeIsUpdated

	    EXPECT_TRUE(myIngredients.isLatticeFilledUp());

	    //manually setting it to false
	    myIngredients.setLatticeFilledUp(false);
	    EXPECT_FALSE(myIngredients.isLatticeFilledUp());

	    //manually setting it to true
	    myIngredients.setLatticeFilledUp(true);
	    EXPECT_TRUE(myIngredients.isLatticeFilledUp());

	    //before synchronize: latticeIsNotUpdated
	    myIngredients.synchronize(myIngredients);
	    EXPECT_TRUE(myIngredients.isLatticeFilledUp());

}
