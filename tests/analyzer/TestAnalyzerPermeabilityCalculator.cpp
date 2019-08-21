/*--------------------------------------------------------------------------------
 *    ooo      L   attice-based  |
 *  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 * o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
 * oo---0---oo  A   lgorithm and  |
 * o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
 *  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
 *    ooo                        |
 * ----------------------------------------------------------------------------------
 *
 * This file is part of LeMonADE.
 *
 * LeMonADE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LeMonADE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * --------------------------------------------------------------------------------*/

/*****************************************************************************/
/**
 * @file
 * @brief Tests for the class AnalyzerPermeabilityCalculator
 *
 * @author Ankush Checkervarty
 * @date 09.08.2019
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/analyzer/AnalyzerPermeabilityCalculator.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/DepthIterator.h>


using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class AnalyzerPermeabilityCalculatorTest
 * @brief prepare system and check functions for AnalyzerPermeabilityCalculatorTest
 * */
/*****************************************************************************/
class AnalyzerPermeabilityCalculatorTest: public ::testing::Test{
protected:
    void setInitConfig()
    {
            ingredients.setBoxX(64);
            ingredients.setBoxY(64);
            ingredients.setBoxZ(64);
            ingredients.setPeriodicX(true);
            ingredients.setPeriodicY(true);
            ingredients.setPeriodicZ(true);
            ingredients.synchronize(ingredients);
    }

    /** This function setups the system with two lipid monomers, two 
     *  disconnected object monomers and a solvent monomer.
     *  First two particles are objects and 3rd one is solvent. 
     */ 
	void setConfig(int32_t midplane, int32_t width)
	{

            ingredients.modifyMolecules().resize(5);
            
            ingredients.modifyMolecules()[0].setAllCoordinates(2,3,midplane+width);
            ingredients.modifyMolecules()[0].setAttributeTag(1);
            
            ingredients.modifyMolecules()[1].setAllCoordinates(2,3,midplane-width);
            ingredients.modifyMolecules()[1].setAttributeTag(2);
            
            ingredients.modifyMolecules()[2].setAllCoordinates(2,5,midplane+width);
            ingredients.modifyMolecules()[2].setAttributeTag(7);
            
            ingredients.modifyMolecules()[3].setAllCoordinates(2,5,midplane-width);
            ingredients.modifyMolecules()[3].setAttributeTag(8);
            
            ingredients.modifyMolecules()[4].setAllCoordinates(2,5,midplane+factor2sigma*width+1);
            ingredients.modifyMolecules()[4].setAttributeTag(5);
            ingredients.synchronize(ingredients);
	}
	
        
	//! set the system into one of two model configs, designed to give easily calculated Rg^2
	std::string passComments(std::ifstream& file)
	{
        std::string linecontent;
        
        for(int32_t i=0;;i++){
            std::getline(file,linecontent);
            if(linecontent[0]!='#'&&linecontent.size()>0)
                break;}
        
        return linecontent;
	}

	//define system
	typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureAttributes<>) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients < Config> MyIngredients;
	MyIngredients ingredients;
    int32_t factor2sigma;

	//! nested derived class from AnalyzerPermeabilityCalculator to test protected functions
	class PmAnalyzerDerived:public AnalyzerPermeabilityCalculator<MyIngredients>
	{
	public:
        using AnalyzerPermeabilityCalculator::AnalyzerPermeabilityCalculator;
        using AnalyzerPermeabilityCalculator::execute;
        using AnalyzerPermeabilityCalculator::particlesInsideBoundaries;
        using AnalyzerPermeabilityCalculator::particlesInsideBoundariesOld;
        using AnalyzerPermeabilityCalculator::counterSolvent;
        using AnalyzerPermeabilityCalculator::counterObjects;
        using AnalyzerPermeabilityCalculator::midplane;
        using AnalyzerPermeabilityCalculator::sigma;
        
	};
	/* suppress cout output for better readability -->un/comment here:*/
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


TEST_F(AnalyzerPermeabilityCalculatorTest, CheckConstructorandPermeabilityupdate)
{
    setInitConfig();
    
    int32_t midplane=32;
    int32_t width=4;
    factor2sigma=3;
    
    setConfig(midplane,width);
    /**setConfig() produces
     * Particle 0:object1 inside the boundaries from positive side.
     * Particle 1:object2 inside the boundaries from negative side.
     * Particle 2:solvent outside the boundaries.
     */
    std::vector<MonomerGroup< MyIngredients::molecules_type> > objects;
    hasOneOfTheseTwoTypes<7,8> hastype;
    fill_connected_groups(ingredients.getMolecules(),objects,MonomerGroup<MyIngredients::molecules_type>((ingredients.getMolecules())),hastype);	
    
    
    //making a bin size same as box, so there is only one bin.
    PmAnalyzerDerived pmAnalyzer(ingredients,objects,64);
    
    //checking the constructor
    EXPECT_EQ(midplane,pmAnalyzer.midplane[0][0]);
    EXPECT_EQ(width,pmAnalyzer.sigma);
    
    EXPECT_EQ(3,pmAnalyzer.particlesInsideBoundariesOld.size());
    EXPECT_EQ(1,pmAnalyzer.particlesInsideBoundariesOld[0]);
    EXPECT_EQ(-1,pmAnalyzer.particlesInsideBoundariesOld[1]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[2]);

    EXPECT_EQ(3,pmAnalyzer.particlesInsideBoundaries.size());
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundaries[0]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundaries[1]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundaries[2]);
    
    EXPECT_EQ(0,pmAnalyzer.counterObjects);
    EXPECT_EQ(0,pmAnalyzer.counterSolvent);
    
    //moving objects out of boundaries to opposite sides.
    //and moving solvent positive side from boundary.
    ingredients.modifyMolecules()[2].setZ(midplane-factor2sigma*width-1);
    ingredients.modifyMolecules()[3].setZ(midplane+factor2sigma*width+1);
    ingredients.modifyMolecules()[4].setZ(midplane+width);

    pmAnalyzer.execute();

    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[0]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[1]);
    EXPECT_EQ(1,pmAnalyzer.particlesInsideBoundariesOld[2]);
    
    //Two object translocated
    EXPECT_EQ(2,pmAnalyzer.counterObjects);
    EXPECT_EQ(0,pmAnalyzer.counterSolvent);
    
    //Objects move back, solvent is opposite side of midplane
    ingredients.modifyMolecules()[2].setZ(midplane-(factor2sigma-1)*width);
    ingredients.modifyMolecules()[3].setZ(midplane+(factor2sigma-1)*width);
    ingredients.modifyMolecules()[4].setZ(midplane-(factor2sigma-1)*width);
    
    pmAnalyzer.execute();

    EXPECT_EQ(-1,pmAnalyzer.particlesInsideBoundariesOld[0]);
    EXPECT_EQ(1,pmAnalyzer.particlesInsideBoundariesOld[1]);
    
    //Analyzer remember where solvent entered into boundaries
    EXPECT_EQ(1,pmAnalyzer.particlesInsideBoundariesOld[2]);
    
    EXPECT_EQ(2,pmAnalyzer.counterObjects);
    EXPECT_EQ(0,pmAnalyzer.counterSolvent);
    
    //Objects move to other side of midplane, solvent 
    //moved out opposite side from where entered in.
    ingredients.modifyMolecules()[2].setZ(midplane+(factor2sigma-1)*width);
    ingredients.modifyMolecules()[3].setZ(midplane-(factor2sigma-1)*width);
    ingredients.modifyMolecules()[4].setZ(midplane-factor2sigma*width-1);

    pmAnalyzer.execute();
    //analyzer remembers where it entered through
    EXPECT_EQ(-1,pmAnalyzer.particlesInsideBoundariesOld[0]);
    EXPECT_EQ(1,pmAnalyzer.particlesInsideBoundariesOld[1]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[2]);
    
    EXPECT_EQ(2,pmAnalyzer.counterObjects);
    EXPECT_EQ(1,pmAnalyzer.counterSolvent);
    
    //Objects move out of boundaries however one of them from
    //same side that it came in. Thus only one event is counted. 
    ingredients.modifyMolecules()[2].setZ(midplane-factor2sigma*width-1);
    ingredients.modifyMolecules()[3].setZ(midplane-factor2sigma*width-3);
    ingredients.modifyMolecules()[4].setZ(midplane-factor2sigma*width-5);

    pmAnalyzer.execute();

    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[0]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[1]);
    EXPECT_EQ(0,pmAnalyzer.particlesInsideBoundariesOld[2]);
    
    EXPECT_EQ(3,pmAnalyzer.counterObjects);
    EXPECT_EQ(1,pmAnalyzer.counterSolvent);

}

TEST_F(AnalyzerPermeabilityCalculatorTest, CheckFileDumpandMidplaneUpdate)
{
    setInitConfig();
    int32_t midplane=32;
    int32_t width=4;
    
    setConfig(midplane,width);

    std::vector<MonomerGroup< MyIngredients::molecules_type> > objects;
    hasOneOfTheseTwoTypes<7,8> hastype;
    fill_connected_groups(ingredients.getMolecules(),objects,MonomerGroup<MyIngredients::molecules_type>((ingredients.getMolecules())),hastype);	
    
    //making a bin size same as box, so there is only one bin.
    std::string outputFilename="TestDump.dat";
    PmAnalyzerDerived pmAnalyzer1(ingredients,objects,64,3,outputFilename,1,150,1);

    pmAnalyzer1.counterObjects=1;
    pmAnalyzer1.counterSolvent=1;
    
    pmAnalyzer1.execute();
    
    std::ifstream file(outputFilename.c_str());
    
    EXPECT_TRUE(file.is_open());
    
    std::string linecontent = passComments(file);
    EXPECT_EQ(linecontent[0],'0');
    EXPECT_EQ(linecontent[2],'1');
    EXPECT_EQ(linecontent[4],'1');
    
    pmAnalyzer1.counterObjects=2;
    pmAnalyzer1.counterSolvent=3;
    

    ingredients.modifyMolecules().setAge(1);
    pmAnalyzer1.execute();
    
    std::getline(file,linecontent);
    EXPECT_EQ(linecontent[0],'1');
    EXPECT_EQ(linecontent[2],'2');
    EXPECT_EQ(linecontent[4],'3');

    remove("TestDump.dat");    
    //midplane update Check
    PmAnalyzerDerived pmAnalyzer2(ingredients,objects,64,3,outputFilename,100,1,1);
    
    EXPECT_EQ(midplane,pmAnalyzer2.midplane[0][0]);
    EXPECT_EQ(width,pmAnalyzer2.sigma);
    
    midplane=35;
    width=2;
    setConfig(midplane,width);
    pmAnalyzer2.execute();
    
    EXPECT_EQ(midplane,pmAnalyzer2.midplane[0][0]);
    EXPECT_EQ(width,pmAnalyzer2.sigma);
    
    midplane=28;
    width=6;
    setConfig(midplane,width);
    pmAnalyzer2.execute();

    EXPECT_EQ(midplane,pmAnalyzer2.midplane[0][0]);
    EXPECT_EQ(width,pmAnalyzer2.sigma);

}

TEST_F(AnalyzerPermeabilityCalculatorTest, CheckPoreUpdate)
{
    setInitConfig();
    int32_t midplane=32;
    int32_t width=4;
    setConfig(midplane,width);
    
    std::vector<MonomerGroup< MyIngredients::molecules_type> > objects;
    hasOneOfTheseTwoTypes<7,8> hastype;
    fill_connected_groups(ingredients.getMolecules(),objects,MonomerGroup<MyIngredients::molecules_type>((ingredients.getMolecules())),hastype);
    
    PmAnalyzerDerived pmAnalyzer3(ingredients,objects,32);
    
    EXPECT_EQ(midplane, pmAnalyzer3.midplane[0][0]);
    EXPECT_EQ(midplane, pmAnalyzer3.midplane[0][1]);
    EXPECT_EQ(midplane, pmAnalyzer3.midplane[1][0]);
    EXPECT_EQ(midplane, pmAnalyzer3.midplane[1][1]);
    
    midplane=35;
    width=4;
    setConfig(midplane,width);
    
    std::string outputFilename="TestDump.dat";
    PmAnalyzerDerived pmAnalyzer4(ingredients,objects,16,3,outputFilename,100,1,1);

    EXPECT_EQ(midplane, pmAnalyzer4.midplane[0][0]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[0][1]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[1][0]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[1][1]);

    EXPECT_NE(midplane, pmAnalyzer4.midplane[2][0]);
    EXPECT_NE(midplane, pmAnalyzer4.midplane[0][2]);
    EXPECT_NE(midplane, pmAnalyzer4.midplane[2][2]);
    
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[3][0]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[3][1]);
    
    ingredients.modifyMolecules().addMonomer(2*16+1,0,midplane+width);
    ingredients.modifyMolecules()[5].setAttributeTag(1);
    ingredients.modifyMolecules().addMonomer(2*16+1,0,midplane-width);
    ingredients.modifyMolecules()[6].setAttributeTag(2);
    
    pmAnalyzer4.execute();
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[2][0]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[2][1]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[3][0]);
    EXPECT_EQ(midplane, pmAnalyzer4.midplane[3][1]);
    remove("TestDump.dat");
        
}
