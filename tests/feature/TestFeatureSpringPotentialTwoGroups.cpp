/*****************************************************************************/
/**
 * @file
 * @brief test for Feature Virtual Spring towards a wall in z-direction
 * @author Martin
 * @date 31.08.2017
 * */
/*****************************************************************************/
#include "gtest/gtest.h"


#include <iostream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

#include "FeatureSpringPotentialTwoGroups.h"


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


//   TEST_F(TestFeatureSpringPotentialTwoGroups,Constructor){
//     // test default constructor
//     FeatureVirtualSpringWall feature;
//     EXPECT_DOUBLE_EQ(0.0,feature.getEquilibriumLength());
//     EXPECT_DOUBLE_EQ(0.0,feature.getSpringConstant());
//     EXPECT_EQ(0,feature.getAffectedMonomerType());
//     
//     //test constructor in ingredients
//     EXPECT_DOUBLE_EQ(0.0,ingredients.getEquilibriumLength());
//     EXPECT_DOUBLE_EQ(0.0,ingredients.getSpringConstant());
//     EXPECT_EQ(0,ingredients.getAffectedMonomerType());
//     
//   }
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
