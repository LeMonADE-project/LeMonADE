#include <iostream>

#include "gtest/gtest.h"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/feature/FeatureNNInteractionBcc.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>

using namespace std;

class NNInteractionBccTest: public ::testing::Test{
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


TEST_F(NNInteractionBccTest,CheckApplyBccMovePowerOfTwoLattice)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionBcc<FeatureLatticePowerOfTwo>) Features1;
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
    //all used moves are of type MoveLocalBcc
    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(10,10,10);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);

    molecules1[1].setAllCoordinates(13,7,7);
    myIngredients1.synchronize(myIngredients1);

    //define moves in all possible directions for later use
    MoveLocalBcc DXDYDZ,DXDYmDZ,DXmDYDZ,DXmDYmDZ,mDXDYDZ,mDXDYmDZ,mDXmDYDZ,mDXmDYmDZ;

    DXDYDZ.init(myIngredients1);
    while((DXDYDZ.getDir()[0]!=1 ) ||(DXDYDZ.getDir()[1]!=1 ) ||(DXDYDZ.getDir()[2]!=1 ) || (DXDYDZ.getIndex()!=1)){ DXDYDZ.init(myIngredients1);}

    DXDYmDZ.init(myIngredients1);
    while((DXDYmDZ.getDir()[0]!=1 ) ||(DXDYmDZ.getDir()[1]!=1 ) ||(DXDYmDZ.getDir()[2]!=-1 ) || (DXDYmDZ.getIndex()!=1)){ DXDYmDZ.init(myIngredients1);}

    DXmDYDZ.init(myIngredients1);
    while((DXmDYDZ.getDir()[0]!=1 ) ||(DXmDYDZ.getDir()[1]!=-1 ) ||(DXmDYDZ.getDir()[2]!=1 ) || (DXmDYDZ.getIndex()!=1)){ DXmDYDZ.init(myIngredients1);}

    DXmDYmDZ.init(myIngredients1);
    while((DXmDYmDZ.getDir()[0]!=1 ) ||(DXmDYmDZ.getDir()[1]!=-1 ) ||(DXmDYmDZ.getDir()[2]!=-1 ) || (DXmDYmDZ.getIndex()!=1)){ DXmDYmDZ.init(myIngredients1);}

    mDXDYDZ.init(myIngredients1);
    while((mDXDYDZ.getDir()[0]!=-1 ) ||(mDXDYDZ.getDir()[1]!=1 ) ||(mDXDYDZ.getDir()[2]!=1 ) || (mDXDYDZ.getIndex()!=1)){ mDXDYDZ.init(myIngredients1);}

    mDXDYmDZ.init(myIngredients1);
    while((mDXDYmDZ.getDir()[0]!=-1 ) ||(mDXDYmDZ.getDir()[1]!=1 ) ||(mDXDYmDZ.getDir()[2]!=-1 ) || (mDXDYmDZ.getIndex()!=1)){ mDXDYmDZ.init(myIngredients1);}

    mDXmDYDZ.init(myIngredients1);
    while((mDXmDYDZ.getDir()[0]!=-1 ) ||(mDXmDYDZ.getDir()[1]!=-1 ) ||(mDXmDYDZ.getDir()[2]!=1 ) || (mDXmDYDZ.getIndex()!=1)){ mDXmDYDZ.init(myIngredients1);}

    mDXmDYmDZ.init(myIngredients1);
    while((mDXmDYmDZ.getDir()[0]!=-1 ) ||(mDXmDYmDZ.getDir()[1]!=-1 ) ||(mDXmDYmDZ.getDir()[2]!=-1 ) || (mDXmDYmDZ.getIndex()!=1)){ mDXmDYmDZ.init(myIngredients1);}

    //now start checking moves of monomer 1 around the position of monomer 0

    ///--- first move upwards in a zig-zag line from (3,-3,-3) to (3,-1,3) relative to monomer 0
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //last step to the next row
    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //----- now move downwards in zig-zag from (3,-1,3) to (3,1,-3)
    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(-epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(-epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(-epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //------ now move upwards in zig-zag to (1,1,3)

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //last step to the next row
    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //--- now again zig-za downwards to (-1,3,-3)
    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(-epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //---move up again to (-3,1,3) relative to monomer 0
    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(-epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(-epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(-epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //---move down to (-2,0,-2), which has a contact, then move to (-1,1,-3), which doesn't
    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);


    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //--  now check the position below monomer 0, i.e. (0,0,-2)
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(-epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--now check the position (0,-2,-2), go there in a first step
    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),7);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(-epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--now check (-2,-2,-2), go there in a first step
    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),7);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--check -2,-2,0, go there in a first step
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(-epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    //-- check (0,-2,0), go there first
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    //--now check (0,-2,2), go there in one step first
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //-- now the last unchecked position with a contact is (0,0,2)
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),13);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);
}

//same test but with any lattice

TEST_F(NNInteractionBccTest,CheckApplyBccMoveAnyLattice)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionBcc<FeatureLattice>) Features1;
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
    //all used moves are of type MoveLocalBcc
    typename Ing1::molecules_type& molecules1=myIngredients1.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(10,10,10);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);

    molecules1[1].setAllCoordinates(13,7,7);
    myIngredients1.synchronize(myIngredients1);

    //define moves in all possible directions for later use
    MoveLocalBcc DXDYDZ,DXDYmDZ,DXmDYDZ,DXmDYmDZ,mDXDYDZ,mDXDYmDZ,mDXmDYDZ,mDXmDYmDZ;

    DXDYDZ.init(myIngredients1);
    while((DXDYDZ.getDir()[0]!=1 ) ||(DXDYDZ.getDir()[1]!=1 ) ||(DXDYDZ.getDir()[2]!=1 ) || (DXDYDZ.getIndex()!=1)){ DXDYDZ.init(myIngredients1);}

    DXDYmDZ.init(myIngredients1);
    while((DXDYmDZ.getDir()[0]!=1 ) ||(DXDYmDZ.getDir()[1]!=1 ) ||(DXDYmDZ.getDir()[2]!=-1 ) || (DXDYmDZ.getIndex()!=1)){ DXDYmDZ.init(myIngredients1);}

    DXmDYDZ.init(myIngredients1);
    while((DXmDYDZ.getDir()[0]!=1 ) ||(DXmDYDZ.getDir()[1]!=-1 ) ||(DXmDYDZ.getDir()[2]!=1 ) || (DXmDYDZ.getIndex()!=1)){ DXmDYDZ.init(myIngredients1);}

    DXmDYmDZ.init(myIngredients1);
    while((DXmDYmDZ.getDir()[0]!=1 ) ||(DXmDYmDZ.getDir()[1]!=-1 ) ||(DXmDYmDZ.getDir()[2]!=-1 ) || (DXmDYmDZ.getIndex()!=1)){ DXmDYmDZ.init(myIngredients1);}

    mDXDYDZ.init(myIngredients1);
    while((mDXDYDZ.getDir()[0]!=-1 ) ||(mDXDYDZ.getDir()[1]!=1 ) ||(mDXDYDZ.getDir()[2]!=1 ) || (mDXDYDZ.getIndex()!=1)){ mDXDYDZ.init(myIngredients1);}

    mDXDYmDZ.init(myIngredients1);
    while((mDXDYmDZ.getDir()[0]!=-1 ) ||(mDXDYmDZ.getDir()[1]!=1 ) ||(mDXDYmDZ.getDir()[2]!=-1 ) || (mDXDYmDZ.getIndex()!=1)){ mDXDYmDZ.init(myIngredients1);}

    mDXmDYDZ.init(myIngredients1);
    while((mDXmDYDZ.getDir()[0]!=-1 ) ||(mDXmDYDZ.getDir()[1]!=-1 ) ||(mDXmDYDZ.getDir()[2]!=1 ) || (mDXmDYDZ.getIndex()!=1)){ mDXmDYDZ.init(myIngredients1);}

    mDXmDYmDZ.init(myIngredients1);
    while((mDXmDYmDZ.getDir()[0]!=-1 ) ||(mDXmDYmDZ.getDir()[1]!=-1 ) ||(mDXmDYmDZ.getDir()[2]!=-1 ) || (mDXmDYmDZ.getIndex()!=1)){ mDXmDYmDZ.init(myIngredients1);}

    //now start checking moves of monomer 1 around the position of monomer 0

    ///--- first move upwards in a zig-zag line from (3,-3,-3) to (3,-1,3) relative to monomer 0
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //last step to the next row
    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //----- now move downwards in zig-zag from (3,-1,3) to (3,1,-3)
    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(-epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(-epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(-epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //------ now move upwards in zig-zag to (1,1,3)

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),13);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),12);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //last step to the next row
    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //--- now again zig-za downwards to (-1,3,-3)
    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(-epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    mDXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYmDZ.getProbability(),exp(epsilon0));
    mDXDYmDZ.apply(myIngredients1);
    mDXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),13);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //---move up again to (-3,1,3) relative to monomer 0
    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(-epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),8);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(-epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYDZ.getProbability(),exp(-epsilon0));
    DXDYDZ.apply(myIngredients1);
    DXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),12);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),13);

    //---move down to (-2,0,-2), which has a contact, then move to (-1,1,-3), which doesn't
    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(-epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),10);

    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);


    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),11);
    EXPECT_EQ(molecules1[1].getZ(),7);

    //--  now check the position below monomer 0, i.e. (0,0,-2)
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(-epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--now check the position (0,-2,-2), go there in a first step
    DXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYmDZ.getProbability(),exp(epsilon0));
    DXmDYmDZ.apply(myIngredients1);
    DXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),7);

    mDXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYDZ.getProbability(),exp(-epsilon0));
    mDXmDYDZ.apply(myIngredients1);
    mDXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--now check (-2,-2,-2), go there in a first step
    mDXmDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXmDYmDZ.getProbability(),exp(epsilon0));
    mDXmDYmDZ.apply(myIngredients1);
    mDXmDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),7);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),8);

    //--check -2,-2,0, go there in a first step
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),7);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),9);

    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(-epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),8);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    //-- check (0,-2,0), go there first
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),10);

    //--now check (0,-2,2), go there in one step first
    DXmDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXmDYDZ.getProbability(),exp(epsilon0));
    DXmDYDZ.apply(myIngredients1);
    DXmDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),11);
    EXPECT_EQ(molecules1[1].getY(),7);
    EXPECT_EQ(molecules1[1].getZ(),11);

    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(-epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),8);
    EXPECT_EQ(molecules1[1].getZ(),12);

    //-- now the last unchecked position with a contact is (0,0,2)
    mDXDYDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(mDXDYDZ.getProbability(),exp(epsilon0));
    mDXDYDZ.apply(myIngredients1);
    mDXDYDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),9);
    EXPECT_EQ(molecules1[1].getY(),9);
    EXPECT_EQ(molecules1[1].getZ(),13);

    DXDYmDZ.check(myIngredients1);
    EXPECT_DOUBLE_EQ(DXDYmDZ.getProbability(),exp(-epsilon0));
    DXDYmDZ.apply(myIngredients1);
    DXDYmDZ.resetProbability();
    EXPECT_EQ(molecules1[1].getX(),10);
    EXPECT_EQ(molecules1[1].getY(),10);
    EXPECT_EQ(molecules1[1].getZ(),12);
}





TEST_F(NNInteractionBccTest,getSetInteraction)
{
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureNNInteractionBcc<FeatureLattice>) Features1;
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

TEST_F(NNInteractionBccTest,Synchronize)
{
    typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureNNInteractionBcc<FeatureLattice>) Features1;
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
    molecules1[0].setAllCoordinates(9,9,9);
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
                if(pos==VectorInt3(9,9,9))
                    EXPECT_EQ(1,myIngredients1.getLatticeEntry(pos));
                else if(pos==VectorInt3(1,1,1))
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

TEST_F(NNInteractionBccTest,ReadWrite)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionBcc<FeatureLattice>) Features1;
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
    molecules1[0].setAllCoordinates(9,9,9);
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


TEST_F(NNInteractionBccTest,ApplyMoveAddMonomerBcc)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureNNInteractionBcc<FeatureLattice>) Features1;
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

    MoveAddMonomerBcc<> addMonomer;

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

    addMonomer.init(myIngredients1);
    addMonomer.setPosition(5,6,7);
    addMonomer.setTag(255);
    addMonomer.apply(myIngredients1);

    EXPECT_EQ(255,int32_t(myIngredients1.getLatticeEntry(5,6,7)));

}
