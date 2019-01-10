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
 * @brief Tests for the classes AnalyzerRadiusOfGyration
 *
 * @author Hauke Rabbel
 * @date 05.07.2016
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>
#include <LeMonADE/utility/Vector3D.h>


using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class AnalyzerRadiusOfGyrationTest
 * @brief prepare system and check functions for AnalyzerRadiusOfGyrationTest
 * */
/*****************************************************************************/
class AnalyzerRadiusOfGyrationTest: public ::testing::Test{
protected:
	//! check if file with given name exists
	bool fileExists(string _filename){
		std::ifstream file(_filename.c_str());
		if(file.is_open()){
			file.close();
			return true;
		}
		else return false;
	}
	//! count the number of lines containing data (i.e. no comments beginning with #)
	uint32_t fileNLines(string _filename){
		std::ifstream file(_filename.c_str());
		if(file.is_open()){
			uint32_t linecount=0;
			while(file.eof()==false){
				std::string comment;
				std::getline(file,comment);
				if(comment.size()>0 && comment[0]!='#')
					linecount++;
			}
			file.close();
			return linecount;
		}
		else return 0;

	}

	//! get the value from lineNo,columnNo in file given
	double getValueFromFile(string _filename,size_t lineNo,size_t columnNo){
		std::ifstream file(_filename.c_str());
		if(file.is_open()){
			std::string line;
			uint32_t linecount=0;
			while(file.eof()==false){

				std::getline(file,line);
				if(line.size()>0 && line[0]!='#'){
					if(linecount==lineNo) break;
					linecount++;
				}
			}
			std::stringstream lineStream(line);
			uint32_t columncount=0;
			double value;
			while(columncount<=columnNo && lineStream.good()){
				lineStream>>value;
				columncount++;
			}
			if(!lineStream.good())
				throw std::runtime_error("AnalyzerRadiusOfGyrationTest::getValue(): could not find column");

			file.close();
			return value;
		}
		else throw std::runtime_error("AnalyzerRadiusOfGyrationTest::getValue(): could not open file");

	}

	//! set the system into one of two model configs, designed to give easily calculated Rg^2
	void setConfig1()
	{
		ingredients.modifyMolecules().resize(4);
		ingredients.modifyMolecules()[0].setAllCoordinates(-2,0,0);
		ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
		ingredients.modifyMolecules()[2].setAllCoordinates(0,-2,0);
		ingredients.modifyMolecules()[3].setAllCoordinates(0,2,0);
	}

	//! set the system into one of two model configs, designed to give easily calculated Rg^2
	void setConfig2()
	{
		ingredients.modifyMolecules().resize(4);
		ingredients.modifyMolecules()[0].setAllCoordinates(-3,0,0);
		ingredients.modifyMolecules()[1].setAllCoordinates(3,0,0);
		ingredients.modifyMolecules()[2].setAllCoordinates(0,-3,0);
		ingredients.modifyMolecules()[3].setAllCoordinates(0,3,0);
	}
	//define system
	typedef LOKI_TYPELIST_1(FeatureAttributes<>) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients < Config> MyIngredients;
	MyIngredients ingredients;

	//! nested derived class from AnalyzerRadiusOfGyration to test protected function
	class RgAnalyzerDerived:public AnalyzerRadiusOfGyration<MyIngredients>
	{
	public:
		RgAnalyzerDerived(const MyIngredients& ing)
		:AnalyzerRadiusOfGyration<MyIngredients >(ing){}

		void setMonomerGroups(std::vector<MonomerGroup<MyIngredients::molecules_type> > groupVector)
		{
			AnalyzerRadiusOfGyration<MyIngredients>::setMonomerGroups(groupVector);
		}

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

/*****************************************************************************/
/**
 * @fn TEST_F(AnalyzerRadiusOfGyrationTest, CheckConstructor)
 * @brief Check if constructor uses the correct default file name
 * */
/*****************************************************************************/
TEST_F(AnalyzerRadiusOfGyrationTest, CheckConstructor)
{
	AnalyzerRadiusOfGyration<MyIngredients> analyzer(ingredients);
	analyzer.initialize();
	analyzer.cleanup();

	//expect the standard output file to exist
	EXPECT_TRUE(fileExists("Rg2TimeSeries.dat"));
	remove("Rg2TimeSeries.dat");
}

/*****************************************************************************/
/**
 * @fn TEST_F(AnalyzerRadiusOfGyrationTest, FileName)
 * @brief Test function for setting output file name
 * */
/*****************************************************************************/
TEST_F(AnalyzerRadiusOfGyrationTest, FileName)
{
	AnalyzerRadiusOfGyration<MyIngredients> analyzer(ingredients,"outfile.dat");
	analyzer.initialize();
	analyzer.cleanup();

	//expect the standard output file to exist
	EXPECT_TRUE(fileExists("outfile.dat"));

	analyzer.setOutputFile("Rgoutfile.dat");
	analyzer.cleanup();
	EXPECT_TRUE(fileExists("Rgoutfile.dat"));

	EXPECT_FALSE(fileExists("Rg2TimeSeries.dat"));
	remove("outfile.dat");
	remove("Rgoutfile.dat");

}

/*****************************************************************************/
/**
 * @fn TEST_F(AnalyzerRadiusOfGyrationTest, BufferSize)
 * @brief Test if file dump after works after buffer is filled
 * */
/*****************************************************************************/
TEST_F(AnalyzerRadiusOfGyrationTest, BufferSize)
{
	setConfig1();
	AnalyzerRadiusOfGyration<MyIngredients> analyzer(ingredients);
	analyzer.setBufferSize(5);
	analyzer.initialize();

	for(int i=0;i<4;i++){
		analyzer.execute();
		EXPECT_FALSE(fileExists("Rg2TimeSeries.dat"));
	}
	analyzer.execute();
	EXPECT_TRUE(fileExists("Rg2TimeSeries.dat"));
	uint32_t nlines=fileNLines("Rg2TimeSeries.dat");
	EXPECT_EQ(nlines,5);

	for(int i=0;i<5;i++){
		analyzer.execute();
	}
	nlines=fileNLines("Rg2TimeSeries.dat");
	EXPECT_EQ(nlines,10);

	analyzer.setBufferSize(10);

	for(int i=0;i<10;i++){
		analyzer.execute();
	}
	nlines=fileNLines("Rg2TimeSeries.dat");
	EXPECT_EQ(nlines,20);
	analyzer.cleanup();
	remove("Rg2TimeSeries.dat");
}

/*****************************************************************************/
/**
 * @fn TEST_F(AnalyzerRadiusOfGyrationTest, CheckOutputValues)
 * @brief Test if output is calculated correctly based on two model configurations
 * */
/*****************************************************************************/
TEST_F(AnalyzerRadiusOfGyrationTest, CheckOutputValues)
{
	setConfig1();
	AnalyzerRadiusOfGyration<MyIngredients> analyzer(ingredients);
	analyzer.initialize();
	analyzer.execute();
	setConfig2();
	analyzer.execute();
	analyzer.cleanup();
	uint32_t nlines=fileNLines("Rg2TimeSeries.dat");
	EXPECT_EQ(nlines,2);

	double value=getValueFromFile("Rg2TimeSeries.dat",0,1);
	EXPECT_EQ(value,2.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,2);
	EXPECT_EQ(value,2.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,3);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,4);
	EXPECT_EQ(value,4.0);

	value=getValueFromFile("Rg2TimeSeries.dat",1,1);
	EXPECT_EQ(value,4.5);
	value=getValueFromFile("Rg2TimeSeries.dat",1,2);
	EXPECT_EQ(value,4.5);
	value=getValueFromFile("Rg2TimeSeries.dat",1,3);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",1,4);
	EXPECT_EQ(value,9.0);

}

/*****************************************************************************/
/**
 * @fn TEST_F(AnalyzerRadiusOfGyrationTest, CheckGroups)
 * @brief Test if output is calculated correctly based on two different groups
 * */
/*****************************************************************************/
TEST_F(AnalyzerRadiusOfGyrationTest, CheckGroups)
{
	MonomerGroup<typename MyIngredients::molecules_type> group1(ingredients.getMolecules());
	MonomerGroup<typename MyIngredients::molecules_type> group2(ingredients.getMolecules());

	group1.push_back(0);
	group1.push_back(1);
	group2.push_back(2);
	group2.push_back(3);

	std::vector<MonomerGroup<typename MyIngredients::molecules_type> > groupVector1;
	std::vector<MonomerGroup<typename MyIngredients::molecules_type> > groupVector2;
	groupVector1.push_back(group1);
	groupVector2.push_back(group2);

	setConfig1();
	//using the nested derived analyzer class here, which provides a public
	//interface for setting groups
	RgAnalyzerDerived analyzer(ingredients);
	analyzer.setMonomerGroups(groupVector1);
	analyzer.initialize();
	analyzer.execute();
	analyzer.setMonomerGroups(groupVector2);
	analyzer.execute();
	analyzer.cleanup();


	uint32_t nlines=fileNLines("Rg2TimeSeries.dat");
	EXPECT_EQ(nlines,2);

	double value=getValueFromFile("Rg2TimeSeries.dat",0,1);
	EXPECT_EQ(value,4.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,2);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,3);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",0,4);
	EXPECT_EQ(value,4.0);

	value=getValueFromFile("Rg2TimeSeries.dat",1,1);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",1,2);
	EXPECT_EQ(value,4.0);
	value=getValueFromFile("Rg2TimeSeries.dat",1,3);
	EXPECT_EQ(value,0.0);
	value=getValueFromFile("Rg2TimeSeries.dat",1,4);
	EXPECT_EQ(value,4.0);

}
