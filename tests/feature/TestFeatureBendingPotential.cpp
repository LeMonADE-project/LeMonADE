#include <iostream>

#include "gtest/gtest.h"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBendingPotential.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>

using namespace std;
/*****************************************************************************/
/**
 * @file
 * @brief Tests for the class FeatureBendingPotential
 *
 * @author Ankush Checkervarty
 * @date 01.08.2019
 * */
/*****************************************************************************/


class BendingPotentialTest: public ::testing::Test{

public:    
    /** This function setups a linear chain which has only one
     * bondvector=(2,0,0) between monomers.
     */ 
    template < class IngredientsType>
	void setSimpleLinearChain(IngredientsType& ingredients, int32_t length)
	{
            ingredients.modifyMolecules().resize(length);
            
            for(int32_t i=0;i<length;i++){
                ingredients.modifyMolecules()[i].setAllCoordinates(2*i,5,5);
                ingredients.modifyMolecules()[i].setAttributeTag(1);}
                
            for(int32_t i=1;i<length;i++)
                ingredients.modifyMolecules().connect(i,i-1);
                
                ingredients.synchronize(ingredients);
	}
	
    /** This function setups a linear chain with two monomer types 
     * which has only one bondvector=(2,0,0) between monomers.
     */ 
    template < class IngredientsType>
	void setMixedLinearChain(IngredientsType& ingredients, int32_t length)
	{
            ingredients.modifyMolecules().resize(length);
            
            for(int32_t i=0;i<length/2;i++){
                ingredients.modifyMolecules()[i].setAllCoordinates(2*i,5,5);
                ingredients.modifyMolecules()[i].setAttributeTag(1);}
            
            for(int32_t i=length/2;i<length;i++){
                ingredients.modifyMolecules()[i].setAllCoordinates(2*i,5,5);
                ingredients.modifyMolecules()[i].setAttributeTag(2);}
                
            for(int32_t i=1;i<length;i++)
                ingredients.modifyMolecules().connect(i,i-1);
                
                ingredients.synchronize(ingredients);
	}

    
    
  /* suppress cout output for better readability -->un-/comment here:*/
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


TEST_F(BendingPotentialTest, CheckLocalmoveProbability)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureBendingPotential) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> Ing;
    Ing myIngredients;



    //prepare myIngredients.
    myIngredients.setBoxX(64);
    myIngredients.setBoxY(64);
    myIngredients.setBoxZ(64);
    myIngredients.setPeriodicX(1);
    myIngredients.setPeriodicY(1);
    myIngredients.setPeriodicZ(1);
    myIngredients.modifyBondset().addBFMclassicBondset();
    
    //setup a simple Linear chain.
    setSimpleLinearChain(myIngredients,5);
    
    //set lipid chainsEnds
    myIngredients.setChainEnds(myIngredients);
    
    //set interaction for type 1
    myIngredients.setBendingPotential(myIngredients,1,0.5);
    double bpStrength=0.5;

    //setup a move type and direction to move.
    MoveLocalSc move;
    VectorInt3 dir(0,1,0);

    //first let's choose the middle most monomer    
    move.init(myIngredients,2,dir);
    float probMiddleMost=1;
    
    VectorInt3 bondVector1,bondVector2;
    
    //Initial bondvectors at monomer position.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);
    
    float probPresentMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength);
 
    //Future bondvectors at monomer position after the move.
    bondVector1.setAllCoordinates(2,1,0);
    bondVector2.setAllCoordinates(2,-1,0);
 
    float probFutureMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength);
    
    probMiddleMost *= probFutureMonPos/probPresentMonPos;

 
    //Initial bondvectors at monomer position one direction positive to the moved monomer.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);

    float probPresentOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength);
    
    bondVector1.setAllCoordinates(2,-1,0);
    bondVector2.setAllCoordinates(2,0,0);

    float probFutureOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength);

    probMiddleMost *= probFutureOnePositive/probPresentOnePositive;
    
    //Since the changes in bond vectors in positive and negative direction is symmetric 
    //we can multiply same probability(calculated for positive) for negative direction.  

    probMiddleMost *= probFutureOnePositive/probPresentOnePositive;

    move.check(myIngredients);

    //use less than instead of equal to get rid of slight mistmatch problem. 
    EXPECT_LE(move.getProbability()-probMiddleMost,1e-8);
    
    
    //choose the monomer second last to end of the chain.
    move.init(myIngredients,3,dir);
    
    float probSecondLastMon = 1;
    
    //It will have contribution from the monomer position change.
    probSecondLastMon *= probFutureMonPos/probPresentMonPos;

    //and contribution from position one positive to monomer.
    probSecondLastMon *= probFutureOnePositive/probPresentOnePositive;
    
    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probSecondLastMon,1e-8);
    
    //Start of chain monomer is symmetric to second last monomer to end of chain.
    move.init(myIngredients,1,dir);

    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probSecondLastMon,1e-8);
    
    //Monomer at the start of the chain. 
    move.init(myIngredients,0,dir);
    
    float probEndMon=1;
    
    //It has only contribution from monomer position one positive to monomer position.
    probEndMon *= probFutureOnePositive/probPresentOnePositive;
    
    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon,1e-8);
    
    //Similarly symmetric for the other end monomer
    move.init(myIngredients,4,dir);

    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon,1e-8);
}


TEST_F(BendingPotentialTest, CheckforMixedChain)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureBendingPotential) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> Ing;
    Ing myIngredients;



    //prepare myIngredients.
    myIngredients.setBoxX(64);
    myIngredients.setBoxY(64);
    myIngredients.setBoxZ(64);
    myIngredients.setPeriodicX(1);
    myIngredients.setPeriodicY(1);
    myIngredients.setPeriodicZ(1);
    myIngredients.modifyBondset().addBFMclassicBondset();
    
    //setup a mixed Linear chain.
    setMixedLinearChain(myIngredients,6);
    
    //set lipid chainsEnds
    myIngredients.setChainEnds(myIngredients);
    
    //set interaction for type 1
    myIngredients.setBendingPotential(myIngredients,1,0.5);
    myIngredients.setBendingPotential(myIngredients,2,0.2);
    double bpStrength1=0.5;
    double bpStrength2=0.2;

    //setup a move type and direction to move.
    MoveLocalSc move;
    VectorInt3 dir(0,1,0);

    //first middle monomer in type 1 part of the chain
    //This is where bit shift entry chainsEnds is useful.    
    move.init(myIngredients,1,dir);
    float probMiddleMost1=1;
    
    VectorInt3 bondVector1,bondVector2;
    
    //Initial bondvectors at monomer position.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);
    
    float probPresentMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength1);
 
    //Future bondvectors at monomer position after the move.
    bondVector1.setAllCoordinates(2,1,0);
    bondVector2.setAllCoordinates(2,-1,0);
 
    float probFutureMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength1);
    
    probMiddleMost1 *= probFutureMonPos/probPresentMonPos;
    
    //There can't be any other contributions as the first chain 
    //only three monomers with type1 
    move.check(myIngredients);
    
    //use less than instead of equal to get rid of slight mistmatch problem. 
    EXPECT_LE(move.getProbability()-probMiddleMost1,1e-8);
    
    //choose the monomer end of the chain for attribute 1.
    move.init(myIngredients,0,dir);
    
    //Initial bondvectors at monomer position one direction positive to the moved monomer.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);

    float probPresentOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength1);
    
    bondVector1.setAllCoordinates(2,-1,0);
    bondVector2.setAllCoordinates(2,0,0);

    float probFutureOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength1);
   
    
    float probEndMon1=1;
    
    //It has only contribution from monomer position one positive to monomer position.
    probEndMon1 *= probFutureOnePositive/probPresentOnePositive;
    
    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon1,1e-8);
    
    //Similarly symmetric for the other end monomer for type 1.
    move.init(myIngredients,2,dir);

    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon1,1e-8);
    
    /********Type2********/
    
    // middle monomer in Type 2
    move.init(myIngredients,4,dir);
    EXPECT_EQ(myIngredients.getMolecules()[4].getAttributeTag(),2);
    float probMiddleMost2=1;

    //Initial bondvectors at monomer position.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);
    
    //multiply with the bpStrength2 this time.
    probPresentMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength2);
 
    //Future bondvectors at monomer position after the move.
    bondVector1.setAllCoordinates(2,1,0);
    bondVector2.setAllCoordinates(2,-1,0);
 
    probFutureMonPos=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength2);
    
    probMiddleMost2 *= probFutureMonPos/probPresentMonPos;
    
    move.check(myIngredients);
    
    //use less than instead of equal to get rid of slight mistmatch problem. 
    EXPECT_LE(move.getProbability()-probMiddleMost2,1e-8);
    
    
    //choose the monomer end of the chain for attribute 2.
    move.init(myIngredients,3,dir);
    
    //Initial bondvectors at monomer position one direction positive to the moved monomer.
    bondVector1.setAllCoordinates(2,0,0);
    bondVector2.setAllCoordinates(2,0,0);
    
    //bpStrength2
    probPresentOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength2);
    
    bondVector1.setAllCoordinates(2,-1,0);
    bondVector2.setAllCoordinates(2,0,0);

    probFutureOnePositive=exp(-BendingPotentials::simpleHarmonic(bondVector1,bondVector2)*bpStrength2);
   
    
    float probEndMon2=1;
    
    //It has only contribution from monomer position one positive to monomer position.
    probEndMon2 *= probFutureOnePositive/probPresentOnePositive;
    
    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon2,1e-8);
    
    //Similarly symmetric for the other end monomer for type 2.
    move.init(myIngredients,5,dir);

    move.check(myIngredients);
    EXPECT_LE(move.getProbability()-probEndMon2,1e-8);
    
}

TEST_F(BendingPotentialTest, CheckReadWriteAndbondVectorToIndex)
{
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBondset<>,FeatureBendingPotential) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> Ing;
    Ing myIngredients;
    myIngredients.modifyBondset().addBFMclassicBondset();


    //set interaction between type
    myIngredients.setBendingPotential(myIngredients,1,0.8);
    myIngredients.setBendingPotential(myIngredients,2,-0.8);
    myIngredients.setBendingPotential(myIngredients,3,0.0);
    myIngredients.setBendingPotential(myIngredients,9,10.0);

    //prepare myIngredients
    myIngredients.setBoxX(32);
    myIngredients.setBoxY(32);
    myIngredients.setBoxZ(32);
    myIngredients.setPeriodicX(1);
    myIngredients.setPeriodicY(1);
    myIngredients.setPeriodicZ(1);
    
    typename Ing::molecules_type& molecules1=myIngredients.modifyMolecules();
    molecules1.resize(2);
    molecules1[0].setAllCoordinates(9,10,10);
    molecules1[1].setAllCoordinates(1,1,1);
    molecules1[0].setAttributeTag(1);
    molecules1[1].setAttributeTag(2);

    myIngredients.synchronize(myIngredients);

    AnalyzerWriteBfmFile<Ing> outfile("./bpStrength.test",myIngredients,AnalyzerWriteBfmFile<Ing>::NEWFILE);
    outfile.initialize();
    outfile.execute();

    Ing myIngredients2;
    UpdaterReadBfmFile<Ing> infile("./bpStrength.test",myIngredients2,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_FAST);
    infile.initialize();
    
    for(int32_t i=1;i<=10;i++)
        EXPECT_LE(myIngredients.getBendingPotential(i)-myIngredients2.getBendingPotential(i),1e-7);
    

    outfile.closeFile();
    infile.closeFile();
    //remove the temporary file
   EXPECT_EQ(0,remove("./bpStrength.test"));
   
   /***testing the new bondVectorToIndex function**/
   std::map <int32_t,VectorInt3>::const_iterator it;

   std::vector<int32_t> store(180,0);
   
   for(it=myIngredients.getBondset().begin();it!=myIngredients.getBondset().end();it++){
       int32_t iD=myIngredients.bondVectorToIndex(it->second);
       store[iD]++;}
   
   int32_t count1=0;
   for(int32_t i=0;i<store.size();i++)
       if(store[i]==1) count1++;
    
    EXPECT_EQ(108,count1);

}
