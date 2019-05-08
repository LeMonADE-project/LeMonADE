/*****************************************************************************/
/**
 * @file
 * @brief test for Feature Spring Potential
 * @author Martin
 * @date 08.05.2019
 * */
/*****************************************************************************/
#include "gtest/gtest.h"


#include <iostream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

#include <LeMonADE/feature/FeatureSpringPotentialTwoGroups.h>


class TestFeatureSpringPotentialTwoGroups: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureSpringPotentialTwoGroups) Features;
  typedef ConfigureSystem<VectorInt3,Features> ConfigType;
  typedef Ingredients < ConfigType> IngredientsType;
  IngredientsType ingredients;
  
  /* suppress cout output for better readability -->un-/comment here: */
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
  /* ** */
};


TEST_F(TestFeatureSpringPotentialTwoGroups,Basics){
  // test default constructor
  FeatureSpringPotentialTwoGroups feature;
  EXPECT_DOUBLE_EQ(0.0,feature.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.0,feature.getSpringConstant());
  
  //test constructor in ingredients
  EXPECT_DOUBLE_EQ(0.0,ingredients.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.0,ingredients.getSpringConstant());

  //setup a small system
  ingredients.setSpringConstant(0.25);
  ingredients.setEquilibriumLength(4);
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(0,0,0);

  // check the spring parameters and the monomer extension
  EXPECT_DOUBLE_EQ(4.00,ingredients.getEquilibriumLength());
  EXPECT_DOUBLE_EQ(0.25,ingredients.getSpringConstant());
  EXPECT_EQ(0,ingredients.getMolecules()[0].getMonomerGroupTag());

  // check the monomer extension setter and getter
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(1));
  EXPECT_EQ(1,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(2));
  EXPECT_EQ(2,ingredients.getMolecules()[0].getMonomerGroupTag());
  EXPECT_NO_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(0));
  EXPECT_EQ(0,ingredients.getMolecules()[0].getMonomerGroupTag());

  // other values than 0,1,2 are not allowed
  EXPECT_ANY_THROW(ingredients.modifyMolecules()[0].setMonomerGroupTag(3));
  
}
//   
//   TEST_F(TestFeatureSpringPotentialTwoGroups,InitializeGroupSorting){
//     //setup 
//     ingredients.setAffectedMonomerType(1);
//     ingredients.setSpringConstant(0.25);
//     ingredients.setEquilibriumLength(4);
//     ingredients.setBoxX(16);
//     ingredients.setBoxY(16);
//     ingredients.setBoxZ(16);
//     ingredients.setPeriodicX(true);
//     ingredients.setPeriodicY(true);
//     ingredients.setPeriodicZ(true);
//     ingredients.modifyBondset().addBFMclassicBondset();
//   
//      //molecule one
//     ingredients.modifyMolecules().addMonomer(0,0,6);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(1);
//     ingredients.modifyMolecules().addMonomer(2,0,6);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(1);
//     ingredients.modifyMolecules().addMonomer(4,0,6);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(1);
//     ingredients.modifyMolecules().connect(0,1);
//     ingredients.modifyMolecules().connect(1,2);
//     /*   o-o-o   */
//     //molecule two
//     ingredients.modifyMolecules().addMonomer(0,0,2);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(4);
//     ingredients.modifyMolecules().addMonomer(2,0,2);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(4);
//     ingredients.modifyMolecules().addMonomer(4,0,2);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(4);
//     ingredients.modifyMolecules().addMonomer(6,0,2);
//     ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(4);
//     ingredients.modifyMolecules().connect(3,4);
//     ingredients.modifyMolecules().connect(4,5);
//     ingredients.modifyMolecules().connect(5,6);
//     /*   o-o-o-o  */
//     
//     EXPECT_NO_THROW(ingredients.synchronize());
//     EXPECT_EQ(7,ingredients.getMolecules().size());
//     EXPECT_EQ(3,ingredients.getAffectedMonomerGroup().size());
//     EXPECT_EQ(0,ingredients.getAffectedMonomerGroup().at(0));
//     EXPECT_EQ(1,ingredients.getAffectedMonomerGroup().at(1));
//     EXPECT_EQ(2,ingredients.getAffectedMonomerGroup().at(2));
//     
//   }
