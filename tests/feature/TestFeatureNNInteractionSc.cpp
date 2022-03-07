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
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc) Features1;
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
    molecules1[0].setInteractionTag(1);
    molecules1[1].setInteractionTag(2);


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
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc) Features1;
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
    molecules1[0].setInteractionTag(1);
    molecules1[1].setInteractionTag(2);


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
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureNNInteractionSc) Features1;
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
    typedef LOKI_TYPELIST_3(FeatureBondset<>,FeatureNNInteractionSc,FeatureExcludedVolumeSc<>) Features1;
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
    molecules1[0].setInteractionTag(1);
    molecules1[1].setInteractionTag(2);

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
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(9,11,10))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(9,10,11))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(9,11,11))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(10,10,10))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(10,11,10))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(10,10,11))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(10,11,11))
                    EXPECT_EQ(1,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(1,1,1))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(1,2,1))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(1,1,2))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(1,2,2))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(2,1,1))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(2,2,1))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(2,1,2))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else if(pos==VectorInt3(2,2,2))
                    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(pos));
                else
                    EXPECT_EQ(0,myIngredients1.getInteractionLatticeEntry(pos));


            }
        }
    }

    molecules1.resize(3);
    molecules1[2].setAllCoordinates(1,1,1);
    molecules1[2].setInteractionTag(2);
    EXPECT_THROW(myIngredients1.synchronize(myIngredients1),std::runtime_error);


}

TEST_F(NNInteractionScTest,ReadWrite)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc) Features1;
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
    molecules1[0].setInteractionTag(1);
    molecules1[1].setInteractionTag(2);

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
    typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionSc, FeatureExcludedVolumeSc<>) Features1;
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
    addMonomer.setInteractionTag(256);
    EXPECT_THROW(addMonomer.apply(myIngredients1),std::runtime_error);

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setInteractionTag(0);
    EXPECT_NO_THROW(addMonomer.apply(myIngredients1));

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setInteractionTag(-1);
    EXPECT_THROW(addMonomer.apply(myIngredients1),std::runtime_error);

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setInteractionTag(2);
    addMonomer.apply(myIngredients1);

    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(5,6,7));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(6,6,7));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(5,7,7));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(6,7,7));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(5,6,8));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(6,6,8));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(5,7,8));
    EXPECT_EQ(2,myIngredients1.getInteractionLatticeEntry(6,7,8));

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setInteractionTag(255);
    addMonomer.apply(myIngredients1);

    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(5,6,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(6,6,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(5,7,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(6,7,7)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(5,6,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(6,6,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(5,7,8)));
    EXPECT_EQ(255,int32_t(myIngredients1.getInteractionLatticeEntry(6,7,8)));


}

class InteractionTagTest: public ::testing::Test{
  
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureNNInteractionSc) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> MyIngredients;
  
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

TEST_F(InteractionTagTest,MonomerInteractionTagGetSet)
{
  MonomerInteractionTag tag;
  EXPECT_EQ(0,static_cast<int32_t>(tag.getInteractionTag()));
  tag.setInteractionTag(-1);
  //castet to uint8_t
  EXPECT_EQ(255,static_cast<int32_t>(tag.getInteractionTag()));
  tag.setInteractionTag(15);
  EXPECT_EQ(15,static_cast<int32_t>(tag.getInteractionTag()));
}

TEST_F(InteractionTagTest,exportRead)
{
    MyIngredients ingredients;
    MyIngredients::molecules_type const& molecules= ingredients.getMolecules();
    FileImport<MyIngredients> file ("tests/interactionTagTest.test",ingredients);

    //scan file for !mcs and read-in first frame
    file.initialize();

    EXPECT_EQ(molecules[0].getInteractionTag(),0); /*this one has the default value*/
    EXPECT_EQ(molecules[1].getInteractionTag(),5);
    EXPECT_EQ(molecules[2].getInteractionTag(),5);
    EXPECT_EQ(molecules[3].getInteractionTag(),5);
    EXPECT_EQ(molecules[4].getInteractionTag(),2);
    EXPECT_EQ(molecules[6].getInteractionTag(),5);
    EXPECT_EQ(molecules[14].getInteractionTag(),5);
    EXPECT_EQ(molecules[15].getInteractionTag(),2);

}

TEST_F(InteractionTagTest,exportWrite)
{
    MyIngredients setupIngrediens;
    MyIngredients ingredients;

    MyIngredients::molecules_type const& setupMolecules= setupIngrediens.getMolecules();
    MyIngredients::molecules_type const& molecules= ingredients.getMolecules();

    setupIngrediens.setBoxX(100);
    setupIngrediens.setPeriodicX(true);
    setupIngrediens.setBoxY(100);
    setupIngrediens.setPeriodicY(true);
    setupIngrediens.setBoxZ(100);
    setupIngrediens.setPeriodicZ(true);
    setupIngrediens.modifyMolecules().resize(5);
    setupIngrediens.modifyMolecules()[0].setInteractionTag(1);
    setupIngrediens.modifyMolecules()[1].setInteractionTag(1);
    setupIngrediens.modifyMolecules()[4].setInteractionTag(2);
    //write to file and read back in
    AnalyzerWriteBfmFile<MyIngredients> outfile("tests/interactionTagTestOut.test",setupIngrediens,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
    outfile.initialize();

    FileImport<MyIngredients> infile ("tests/interactionTagTestOut.test",ingredients);

    //scan file for !mcs and read-in first frame
    infile.initialize();

    EXPECT_EQ(molecules[0].getInteractionTag(),1);
    EXPECT_EQ(molecules[1].getInteractionTag(),1);
    EXPECT_EQ(molecules[2].getInteractionTag(),0); /*has default value*/
    EXPECT_EQ(molecules[3].getInteractionTag(),0); /*has default value*/
    EXPECT_EQ(molecules[4].getInteractionTag(),2);

    //remove the temporary file
    EXPECT_EQ(0,remove("tests/interactionTagTestOut.test"));
}

TEST_F(InteractionTagTest,ReadInteractionTagsClass)
{
    //this test is for testing the reaction of the read class to incorrectly
    //formatted input
    MyIngredients ingredients;
    ReadInteractionTags<MyIngredients> read(ingredients);

    ingredients.modifyMolecules().resize(10);

    std::stringstream stream1;
    read.setInputStream(&stream1);
    stream1<<"\ni am a parrot\n";
    EXPECT_THROW(read.execute(),std::runtime_error);


    std::stringstream stream2;
    read.setInputStream(&stream2);
    stream2<<"\n1:2:2\n";
    EXPECT_THROW(read.execute(),std::runtime_error);


    std::stringstream stream3;
    read.setInputStream(&stream3);
    stream3<<"\n1-2-3\n";
    EXPECT_THROW(read.execute(),std::runtime_error);


    std::stringstream stream4;
    read.setInputStream(&stream4);
    stream4<<"\n1-a:4\n";
    EXPECT_THROW(read.execute(),std::runtime_error);

    std::stringstream stream5;
    read.setInputStream(&stream5);
    stream5<<"\n1-5:c\n";
    EXPECT_THROW(read.execute(),std::runtime_error);
}

TEST_F(InteractionTagTest,CopyConstructor)
{
	typedef MyIngredients::molecules_type MyMolecules;
	MyMolecules molecules1;

	//check if attributes are copied corriectly by copy constructor
	molecules1.resize(5);
	molecules1[0].setInteractionTag(1);
	molecules1[1].setInteractionTag(2);
	molecules1[2].setInteractionTag(3);
	molecules1[3].setInteractionTag(4);
	molecules1[4].setInteractionTag(5);
	//create new objects from molecules1

	typedef ConfigureSystem<VectorInt3,Features,8> Config8;
	typedef Ingredients<Config8> MyIngredients8;
	typedef MyIngredients8::molecules_type MyMolecules8;
	MyMolecules molecules2(molecules1);
	MyMolecules8 molecules3(molecules2);

	EXPECT_EQ(molecules1[0].getInteractionTag(),molecules2[0].getInteractionTag());
	EXPECT_EQ(molecules1[0].getInteractionTag(),molecules3[0].getInteractionTag());

	EXPECT_EQ(molecules1[1].getInteractionTag(),molecules2[1].getInteractionTag());
	EXPECT_EQ(molecules1[1].getInteractionTag(),molecules3[1].getInteractionTag());

	EXPECT_EQ(molecules1[2].getInteractionTag(),molecules2[2].getInteractionTag());
	EXPECT_EQ(molecules1[2].getInteractionTag(),molecules3[2].getInteractionTag());

	EXPECT_EQ(molecules1[3].getInteractionTag(),molecules2[3].getInteractionTag());
	EXPECT_EQ(molecules1[3].getInteractionTag(),molecules3[3].getInteractionTag());

	EXPECT_EQ(molecules1[4].getInteractionTag(),molecules2[4].getInteractionTag());
	EXPECT_EQ(molecules1[4].getInteractionTag(),molecules3[4].getInteractionTag());

}