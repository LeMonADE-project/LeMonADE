#include <iostream>

#include "gtest/gtest.h"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>

using namespace std;

class NNInteractionScTest: public ::testing::Test{
  /* suppress cout output for better readability -->un-/comment here:*/
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
  /* ** */
};


TEST_F(NNInteractionScTest,CheckApplyScMovePowerOfTwoLattice)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc<FeatureLatticePowerOfTwo>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients1;



    //prepare myIngredients1
    myIngredients1.setBoxX(64);
    myIngredients1.setBoxY(64);
    myIngredients1.setBoxZ(64);
    myIngredients1.setPeriodicX(1);
    myIngredients1.setPeriodicY(1);
    myIngredients1.setPeriodicZ(1);

    //set interaction between types 1,2
    myIngredients1.setNNInteraction(1,2,0.8);
    double epsilon0=0.8;

    //add two monomers. the test then proceeds like this: monomer 1 is moved
    //to all positions around monomer 0 that give contributions to the
    //acceptance probability. the probability is checked using a testmove
    //and compared with the expected value. in the testmove, monomer 0 is
    //simply moved one step in positive x-direction, and then back to its
    //original position.
    //the reverse move is checked using
    //a different testmove (backmove) and also compared with the expectation.
    //the monomer is then moved to a new position using updatemove.
    //all used moves are of type MoveLocalSc
    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(9,10,10);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);


    //now check EV for move in positive x-direction
    MoveLocalSc testmove,backmove,updateMove;

    molecules1[1].setAllCoordinates(12,9,10);
    myIngredients1.synchronize(myIngredients1);

    //cout<<"-------------> checking single Interaction same hight <--------------------"<<std::endl;

    testmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    //std::cout<<"testing move in positive x direction..."<<testmove.getDir()<<" monomer index "<<testmove.getIndex()<<std::endl;
    backmove.init(myIngredients1);
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    //std::cout<<"and move in negativ x direction..."<<backmove.getDir()<<" monomer index "<<backmove.getIndex()<<std::endl;

    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    EXPECT_EQ(molecules1[0].getX(),9);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);
    testmove.apply(myIngredients1);
    EXPECT_EQ(molecules1[0].getX(),10);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    EXPECT_EQ(molecules1[0].getX(),9);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);

    //move monomoer1 to 12 10 10
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-4.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(4.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 12 11 10
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(4.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-4.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to starting position
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube upwards in z-direction
    while((testmove.getDir().getZ()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);

    //cout<<"-------------> checking single Interaction above <--------------------"<<std::endl;
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    //check 12 10 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    //check 12 11 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to starting position
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);


    //now move the testcube 2*downwards in z-direction
    while((testmove.getDir().getZ()!=-1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);

    //cout<<"-------------> checking single Interaction below <--------------------"<<std::endl;
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    //check 12 10 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    //check 12 11 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //this time do not move back to the starting position, but one step
    //in y direction to test the area that was left out. the testcube
    //will for this be moved another step downwards.
    //cout<<"-------------> checking single Interaction twice below <--------------------"<<std::endl;
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube one step further down in z-direction to test the
    //layer of 12 cubes above the testcube
    while((testmove.getDir().getZ()!=-1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);

    //first move in negative x direction four times
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //now move in positive y once and then 4 times positive x


    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move in pos y and then four times negative x again
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to 11 9 10
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube up in z direction four times
    while((testmove.getDir().getZ()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);

    //now do the procedure again
    //cout<<"-------------> checking single Interaction twice above <--------------------"<<std::endl;


    //first move in negative x direction four times
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //now move in positive y once and then 4 times positive x


    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move in pos y and then four times negative x again
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);
}




//same test again with no power of two lattice

TEST_F(NNInteractionScTest,CheckApplyScMoveAnyLattice)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc<FeatureLattice>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients1;



    //prepare myIngredients1
    myIngredients1.setBoxX(65);
    myIngredients1.setBoxY(65);
    myIngredients1.setBoxZ(65);
    myIngredients1.setPeriodicX(1);
    myIngredients1.setPeriodicY(1);
    myIngredients1.setPeriodicZ(1);

    //set interaction between types 1,2
    myIngredients1.setNNInteraction(1,2,0.8);
    double epsilon0=0.8;

    //add two monomers. the test then proceeds like this: monomer 1 is moved
    //to all positions around monomer 0 that give contributions to the
    //acceptance probability. the probability is checked using a testmove
    //and compared with the expected value. in the testmove, monomer 0 is
    //simply moved one step in positive x-direction, and then back to its
    //original position.
    //the reverse move is checked using
    //a different testmove (backmove) and also compared with the expectation.
    //the monomer is then moved to a new position using updatemove.
    //all used moves are of type MoveLocalSc
    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(9,10,10);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);


    //now check EV for move in positive x-direction
    MoveLocalSc testmove,backmove,updateMove;

    molecules1[1].setAllCoordinates(12,9,10);
    myIngredients1.synchronize(myIngredients1);

    //cout<<"-------------> checking single Interaction same hight <--------------------"<<std::endl;

    testmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    //std::cout<<"testing move in positive x direction..."<<testmove.getDir()<<" monomer index "<<testmove.getIndex()<<std::endl;
    backmove.init(myIngredients1);
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    //std::cout<<"and move in negativ x direction..."<<backmove.getDir()<<" monomer index "<<backmove.getIndex()<<std::endl;

    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    EXPECT_EQ(molecules1[0].getX(),9);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);
    testmove.apply(myIngredients1);
    EXPECT_EQ(molecules1[0].getX(),10);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    EXPECT_EQ(molecules1[0].getX(),9);
    EXPECT_EQ(molecules1[0].getY(),10);
    EXPECT_EQ(molecules1[0].getZ(),10);

    //move monomoer1 to 12 10 10
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-4.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(4.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 12 11 10
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(4.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-4.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to starting position
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube upwards in z-direction
    while((testmove.getDir().getZ()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);

    //cout<<"-------------> checking single Interaction above <--------------------"<<std::endl;
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    //check 12 10 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    //check 12 11 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to starting position
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);


    //now move the testcube 2*downwards in z-direction
    while((testmove.getDir().getZ()!=-1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);

    //cout<<"-------------> checking single Interaction below <--------------------"<<std::endl;
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    //check 12 10 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);
    //check 12 11 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 11 12 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 7 11 10
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 10 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 9 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //check 8 12 10
    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    //this time do not move back to the starting position, but one step
    //in y direction to test the area that was left out. the testcube
    //will for this be moved another step downwards.
    //cout<<"-------------> checking single Interaction twice below <--------------------"<<std::endl;
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube one step further down in z-direction to test the
    //layer of 12 cubes above the testcube
    while((testmove.getDir().getZ()!=-1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);

    //first move in negative x direction four times
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //now move in positive y once and then 4 times positive x


    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move in pos y and then four times negative x again
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //move back to 11 9 10
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getY()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    updateMove.apply(myIngredients1);

    //now move the testcube up in z direction four times
    while((testmove.getDir().getZ()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);
    testmove.apply(myIngredients1);

    //now do the procedure again
    //cout<<"-------------> checking single Interaction twice above <--------------------"<<std::endl;


    //first move in negative x direction four times
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    //now move in positive y once and then 4 times positive x


    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-2.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(2.0*epsilon0));
    backmove.apply(myIngredients1);

    //move in pos y and then four times negative x again
    while(updateMove.getDir().getY()!=1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);
    updateMove.apply(myIngredients1);
    while(updateMove.getDir().getX()!=-1 || updateMove.getIndex()!=1) updateMove.init(myIngredients1);


    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(-1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);

    updateMove.apply(myIngredients1);

    testmove.init(myIngredients1);
    backmove.init(myIngredients1);
    while((testmove.getDir().getX()!=1) || (testmove.getIndex()!=0)){ testmove.init(myIngredients1);}
    while((backmove.getDir().getX()!=-1) || (backmove.getIndex()!=0)){ backmove.init(myIngredients1);}
    testmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(testmove.getProbability(),exp(1.0*epsilon0));
    testmove.apply(myIngredients1);
    backmove.check(myIngredients1);
    EXPECT_DOUBLE_EQ(backmove.getProbability(),exp(-1.0*epsilon0));
    backmove.apply(myIngredients1);
}





TEST_F(NNInteractionScTest,getSetInteraction)
{
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureNNInteractionSc<FeatureLattice>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients;

    for(int32_t i=1;i<=255;i++)
    {
        for(int32_t j=1;j<=255;j++)
        {
            EXPECT_DOUBLE_EQ(myIngredients.getNNInteraction(i,j),0.0);
        }
    }
    EXPECT_THROW(myIngredients.getNNInteraction(-1,5),std::runtime_error);
    EXPECT_THROW(myIngredients.getNNInteraction(1,-5),std::runtime_error);
    EXPECT_THROW(myIngredients.getNNInteraction(0,1),std::runtime_error);
    EXPECT_THROW(myIngredients.getNNInteraction(1,0),std::runtime_error);
    EXPECT_THROW(myIngredients.getNNInteraction(0,256),std::runtime_error);
    EXPECT_THROW(myIngredients.getNNInteraction(256,0),std::runtime_error);

    //set interaction between types 1,2
    myIngredients.setNNInteraction(1,2,0.8);
    EXPECT_DOUBLE_EQ(myIngredients.getNNInteraction(1,2),0.8);

    myIngredients.setNNInteraction(1,2,-0.8);
    EXPECT_DOUBLE_EQ(myIngredients.getNNInteraction(1,2),-0.8);


    myIngredients.setNNInteraction(1,2,0.0);
    EXPECT_DOUBLE_EQ(myIngredients.getNNInteraction(1,2),0.0);
}

TEST_F(NNInteractionScTest,Synchronize)
{
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureNNInteractionSc<FeatureLattice>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients1;
    //prepare myIngredients1
    myIngredients1.setBoxX(32);
    myIngredients1.setBoxY(32);
    myIngredients1.setBoxZ(32);
    myIngredients1.setPeriodicX(1);
    myIngredients1.setPeriodicY(1);
    myIngredients1.setPeriodicZ(1);

    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(9,10,10);
    molecules1[1].setAllCoordinates(1,1,1);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);

    myIngredients1.synchronize(myIngredients1);

    VectorInt3 pos;
    //go through lattice and count contacts at each position
    for(int32_t x=0;x<32;x++)
    {
        for(int32_t y=0;y<32;y++)
        {
            for(int32_t z=0;z<32;z++)
            {
                //get the lattice entry
                pos.setAllCoordinates(x,y,z);
                if(pos==VectorInt3(9,10,10))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(9,11,10))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(9,10,11))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(9,11,11))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(10,10,10))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(10,11,10))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(10,10,11))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(10,11,11))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(1,1,1))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(1,2,1))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(1,1,2))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(1,2,2))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(2,1,1))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(2,2,1))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(2,1,2))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(2,2,2))
                    EXPECT_EQ(2,myIngredients1.getLatticeEntry(pos));
                else
                    EXPECT_EQ(0,myIngredients1.getLatticeEntry(pos));


            }
        }
    }

    molecules1.resize(3);
    molecules1[2].setAllCoordinates(1,1,1);
    molecules1[2].setAttributeTag(2);
    EXPECT_THROW(myIngredients1.synchronize(myIngredients1),std::runtime_error);


}

TEST_F(NNInteractionScTest,ReadWrite)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc<FeatureLattice>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients;


    //set interaction between types 1,2
    myIngredients.setNNInteraction(1,2,0.8);
    myIngredients.setNNInteraction(1,3,-0.8);
    myIngredients.setNNInteraction(2,3,0.0);
    myIngredients.setNNInteraction(9,4,10.0);

    //prepare myIngredients
    myIngredients.setBoxX(32);
    myIngredients.setBoxY(32);
    myIngredients.setBoxZ(32);
    myIngredients.setPeriodicX(1);
    myIngredients.setPeriodicY(1);
    myIngredients.setPeriodicZ(1);

    typename Ing1::molecules_type& molecules1=myIngredients.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(9,10,10);
    molecules1[1].setAllCoordinates(1,1,1);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);

    myIngredients.synchronize(myIngredients);

    AnalyzerWriteBfmFile<Ing1> outfile("./interactionRW.test",myIngredients,AnalyzerWriteBfmFile<Ing1>::NEWFILE);
    outfile.initialize();
    outfile.execute();

    Ing1 myIngredients2;
    UpdaterReadBfmFile<Ing1> infile("./interactionRW.test",myIngredients2,UpdaterReadBfmFile<Ing1>::READ_LAST_CONFIG_FAST);
    infile.initialize();

    for(int32_t i=1;i<=255;i++)
    {
        for(int32_t j=1;j<=255;j++)
        {
            EXPECT_DOUBLE_EQ(myIngredients.getNNInteraction(i,j),myIngredients2.getNNInteraction(i,j));
        }
    }

    outfile.closeFile();
    infile.closeFile();
    //remove the temporary file
  EXPECT_EQ(0,remove("./interactionRW.test"));


}


TEST_F(NNInteractionScTest,ApplyMoveAddMonomerSc)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc<FeatureLattice>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 myIngredients1;

    //prepare myIngredients1
    myIngredients1.setBoxX(65);
    myIngredients1.setBoxY(65);
    myIngredients1.setBoxZ(65);
    myIngredients1.setPeriodicX(1);
    myIngredients1.setPeriodicY(1);
    myIngredients1.setPeriodicZ(1);

    myIngredients1.synchronize();

    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();

    MoveAddMonomerSc<> addMonomer;

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(256);
    EXPECT_THROW(addMonomer.apply(myIngredients1),std::runtime_error);

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(0);
    EXPECT_THROW(addMonomer.apply(myIngredients1),std::runtime_error);

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(-1);
    EXPECT_THROW(addMonomer.apply(myIngredients1),std::runtime_error);

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(2);
    addMonomer.apply(myIngredients1);

    EXPECT_EQ(2,myIngredients1.getLatticeEntry(5,6,7));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(6,6,7));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(5,7,7));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(6,7,7));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(5,6,8));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(6,6,8));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(5,7,8));
    EXPECT_EQ(2,myIngredients1.getLatticeEntry(6,7,8));

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(255);
    addMonomer.apply(myIngredients1);

    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(5,6,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(6,6,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(5,7,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(6,7,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(5,6,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(6,6,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(5,7,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(6,7,8)));


}
