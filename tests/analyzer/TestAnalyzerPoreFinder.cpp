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
 * @brief Tests for the class AnalyzerPoreFinder
 *
 * @author Ankush Checkervarty
 * @date 09.08.2019
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <string> 

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/analyzer/AnalyzerPoreFinder.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/DepthIterator.h>


using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class AnalyzerPoreFinderTest
 * @brief prepare system and check functions for AnalyzerPoreFinderTest
 * */
/*****************************************************************************/
class AnalyzerPoreFinderTest: public ::testing::Test{
protected:
    
    void setInitConfig()
    {
            ingredients.setBoxX(32);
            ingredients.setBoxY(32);
            ingredients.setBoxZ(32);
            ingredients.setPeriodicX(true);
            ingredients.setPeriodicY(true);
            ingredients.setPeriodicZ(true);
            ingredients.synchronize(ingredients);
    }
    
    int32_t reduceDistanceInPeriodicSpace(int32_t distance, int32_t period)
    {
        if(distance > period/2){
            while(-(distance-period)<period/2) 
                distance -= (period);
            
            return distance;}
            
        else if(-distance > period/2){ 
            while((distance+period)<period/2) 
                distance += (period);
            
            return distance;}
            
        else 
            return distance;
    }

    /** This function sets up a pore with a given radius and pore center. 
     */ 
	void setPore(int32_t radius,VectorInt2 center)
	{
            int32_t boxX=ingredients.getBoxX();
            int32_t boxY=ingredients.getBoxY();
            int32_t boxZ=ingredients.getBoxZ();
            int32_t boxXm_1=ingredients.getBoxX()-1;
            int32_t boxYm_1=ingredients.getBoxY()-1;
            int32_t size=0;
            int32_t poreCenterX=center.getX();
            int32_t poreCenterY=center.getY();
            
            for(int32_t x=0;x<boxX;x++)
                for(int32_t y=0;y<boxY;y++){
                    
                    int32_t shift_x=reduceDistanceInPeriodicSpace(x-poreCenterX,boxX);
                    int32_t shift_y=reduceDistanceInPeriodicSpace(y-poreCenterY,boxY);
                    int32_t shift_xp1=shift_x+1;
                    int32_t shift_yp1=shift_y+1;

                    if(shift_x*shift_x+shift_y*shift_y>(radius)*(radius)&&
                       shift_xp1*shift_xp1+shift_y*shift_y>(radius)*(radius)&& 
                       shift_x*shift_x+shift_yp1*shift_yp1>(radius)*(radius)&&                      
                       shift_xp1*shift_xp1+shift_yp1*shift_yp1>(radius)*(radius) 
                    ){
                        ingredients.modifyMolecules().addMonomer(x,y,boxZ/2);
                        ingredients.modifyMolecules()[size].setAttributeTag(2);
                        size++;
                    }
                }
	}
        /** This function sets up two pores with given radii and pore centers.*/ 
	void setTwoPores(int32_t radius1,int32_t radius2, VectorInt2 center1, VectorInt2 center2)
	{
            int32_t boxX=ingredients.getBoxX();
            int32_t boxY=ingredients.getBoxY();
            int32_t boxZ=ingredients.getBoxZ();
            int32_t boxXm_1=ingredients.getBoxX()-1;
            int32_t boxYm_1=ingredients.getBoxY()-1;
            int32_t size=0;
            int32_t poreCenterX1=center1.getX();
            int32_t poreCenterY1=center1.getY();
            int32_t poreCenterX2=center2.getX();
            int32_t poreCenterY2=center2.getY();
            
            for(int32_t x=0;x<boxX;x++)
                for(int32_t y=0;y<boxY;y++){
                    
                    int32_t shift_x1=reduceDistanceInPeriodicSpace(x-poreCenterX1,boxX);
                    int32_t shift_y1=reduceDistanceInPeriodicSpace(y-poreCenterY1,boxY);
                    int32_t shift_x1p1=shift_x1+1;
                    int32_t shift_y1p1=shift_y1+1;
                    
                    int32_t shift_x2=reduceDistanceInPeriodicSpace(x-poreCenterX2,boxX);
                    int32_t shift_y2=reduceDistanceInPeriodicSpace(y-poreCenterY2,boxY);
                    int32_t shift_x2p1=shift_x2+1;
                    int32_t shift_y2p1=shift_y2+1;

                    if(shift_x1*shift_x1+shift_y1*shift_y1>(radius1)*(radius1)&&
                       shift_x1p1*shift_x1p1+shift_y1*shift_y1>(radius1)*(radius1)&& 
                       shift_x1*shift_x1+shift_y1p1*shift_y1p1>(radius1)*(radius1)&&                      
                       shift_x1p1*shift_x1p1+shift_y1p1*shift_y1p1>(radius1)*(radius1)&&
                       
                       shift_x2*shift_x2+shift_y2*shift_y2>(radius2)*(radius2)&&
                       shift_x2p1*shift_x2p1+shift_y2*shift_y2>(radius2)*(radius2)&& 
                       shift_x2*shift_x2+shift_y2p1*shift_y2p1>(radius2)*(radius2)&&                      
                       shift_x2p1*shift_x2p1+shift_y2p1*shift_y2p1>(radius2)*(radius2)
                    ){
                        ingredients.modifyMolecules().addMonomer(x,y,boxZ/2);
                        ingredients.modifyMolecules()[size].setAttributeTag(2);
                        size++;
                    }
                }
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
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO,FeatureLatticePowerOfTwo<int32_t>,FeatureAttributes<>) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients <Config> MyIngredients;
	MyIngredients ingredients;
        
        //!nested derived class from AnalyzerPoreFinder to test protected functions
	class PmAnalyzerDerived: public AnalyzerPoreFinder<MyIngredients>
	{
	public:
        using AnalyzerPoreFinder::AnalyzerPoreFinder;
        using AnalyzerPoreFinder::execute;
        using AnalyzerPoreFinder::cleanup;
        using AnalyzerPoreFinder::coordinatesOfPore;
        using AnalyzerPoreFinder::centriod;
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

TEST_F(AnalyzerPoreFinderTest, ClusterAnalysisOnePore)
{
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    setInitConfig();
    /**create a pore in the middle the box of radius
     */
    int32_t radius=3;
    int32_t poreCenterX=16;
    int32_t poreCenterY=16;
    VectorInt2 center(poreCenterX,poreCenterY);
    
    setPore(radius,center);
    /**setPore(3,16,16) produces
     * monomers with tag 2 all over the middle of box in z direction 
     * in xy plane leaving a radius of 3(pore).
     */
    
    std::vector<MonomerGroup< MyIngredients::molecules_type> > objects;
    hasOneOfTheseTwoTypes<7,8> hastype;
    fill_connected_groups(ingredients.getMolecules(),objects,MonomerGroup<MyIngredients::molecules_type>((ingredients.getMolecules())),hastype);
    
    PmAnalyzerDerived analPore(ingredients,objects);
    
    analPore.execute();

    float pi=3.142;
    int32_t area=pi*radius*radius;
 
    //Since this is comparing discrete lattice points and real space,
    //there is some uncertainty of 3-4 lattice points 
    int32_t diff= area-analPore.coordinatesOfPore.size();
    EXPECT_LE(abs(diff),4);
   
    //Only Centriod should be there 
    EXPECT_EQ(1,analPore.centriod.size());

    //This is again due to difference between discrete and real space 
    EXPECT_LE(abs(analPore.centriod[0].getX()-poreCenterX),1);

    EXPECT_LE(abs(analPore.centriod[0].getY()-poreCenterY),1);
        
    //Now let's test it for random radius and random center position
    //radius can be maximum half of box size
    radius=rand()%5+3;

    poreCenterX=rand()%32;
    poreCenterY=rand()%32;
    VectorInt2 center2(poreCenterX,poreCenterY);
    
    //set Molecules size to zero again.
    ingredients.modifyMolecules().resize(0);
    
    setPore(radius,center2);
    analPore.execute();
    
    area=pi*radius*radius;
    diff= area-analPore.coordinatesOfPore.size();
    EXPECT_LE(abs(diff),4);
   
    //only one Centriod should be there 
    EXPECT_EQ(1,analPore.centriod.size());

    //this is again due to difference between discrete and real space 
    EXPECT_LE(abs(analPore.centriod[0].getX()-poreCenterX),1);

    EXPECT_LE(abs(analPore.centriod[0].getY()-poreCenterY),1);
    
    radius=3;
    poreCenterX=16;
    poreCenterY=16;
    
    
    analPore.execute();
    //check the cleanup and files.
    analPore.cleanup();
    
    std::ifstream file("PoreCoordinates.dat");
    
    EXPECT_TRUE(file.is_open());
    
    std::string linecontent = passComments(file);
    EXPECT_EQ(linecontent[0],'0');
    EXPECT_EQ(linecontent[2],'0');
    EXPECT_EQ(linecontent[4],'0');
     
    for(int32_t x=0;x<16;x++)
        for(int32_t y=0;y<32;y++)
            std::getline(file,linecontent);

    for(int32_t y=0;y<16;y++)
        std::getline(file,linecontent);
     
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"16");
    EXPECT_EQ(linecontent.substr(6,1),"1");
    
    std::getline(file,linecontent);

    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"17");
    EXPECT_EQ(linecontent.substr(6,1),"1");
    
    file.close();
    remove("PoreCoordinates.dat");
    
 }
TEST_F(AnalyzerPoreFinderTest, ClusterAnalysisTwoPoreandCopolymers)
{
    setInitConfig();

    int32_t radius1=3;
    int32_t poreCenterX1=16;
    int32_t poreCenterY1=16;
    VectorInt2 center1(poreCenterX1,poreCenterY1);

    int32_t radius2=4;
    int32_t poreCenterX2=26;
    int32_t poreCenterY2=16;
    VectorInt2 center2(poreCenterX2,poreCenterY2);
    
    setTwoPores(radius1,radius2,center1,center2);

    //let's create some copolymers
    ingredients.modifyMolecules().addMonomer(poreCenterX1,poreCenterY1-9,16+5);
    int32_t last_Index=ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules()[last_Index].setAttributeTag(8);

    ingredients.modifyMolecules().addMonomer(poreCenterX1,poreCenterY1-3,16+5);
    last_Index=ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules()[last_Index].setAttributeTag(8);

    ingredients.modifyMolecules().addMonomer(poreCenterX1,poreCenterY1+3,16+5);
    last_Index=ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules()[last_Index].setAttributeTag(7);
    
    
    std::vector<MonomerGroup< MyIngredients::molecules_type> > objects;
    hasOneOfTheseTwoTypes<7,8> hastype;
    fill_connected_groups(ingredients.getMolecules(),objects,MonomerGroup<MyIngredients::molecules_type>((ingredients.getMolecules())),hastype);
    
    PmAnalyzerDerived analPore(ingredients,objects);
    
    analPore.execute();

    float pi=3.142;
    int32_t area=pi*radius1*radius1 + pi*radius2*radius2;
    int32_t diff= area-analPore.coordinatesOfPore.size();

    EXPECT_LE(abs(diff),4);
    
    //Two Centriods should be there.
    EXPECT_EQ(2,analPore.centriod.size());

    //According to the search pore with less value of x and y.
    EXPECT_LE(abs(analPore.centriod[0].getX()-poreCenterX1),1);

    EXPECT_LE(abs(analPore.centriod[0].getY()-poreCenterY1),1);
        
    EXPECT_LE(abs(analPore.centriod[1].getX()-poreCenterX2),1);

    EXPECT_LE(abs(analPore.centriod[1].getY()-poreCenterY2),1);

    //check the cleanup and files.
    analPore.cleanup();
    
    std::ifstream file("PoreCoordinates.dat");
    
    EXPECT_TRUE(file.is_open());
    
    std::string linecontent = passComments(file);
    
    //check for the values at the center.
    for(int32_t x=0;x<16;x++)
        for(int32_t y=0;y<32;y++)
            std::getline(file,linecontent);
        
    for(int32_t y=0;y<16;y++)
        std::getline(file,linecontent);
     
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"16");
    EXPECT_EQ(linecontent.substr(6,1),"1");
    
    //If there are two pores, the analyzer
    //treats them like two pore events.
    //Thus, value at 4 radius points should 0.5!
    for(int32_t y=16;y<16+4;y++)
        std::getline(file,linecontent);
    
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"20");
    EXPECT_EQ(linecontent.substr(6,3),"0.5");
    
    file.close();
    
    //check for copolymers.
    std::ifstream file1("PolymerCoordinates.dat");
    
    EXPECT_TRUE(file1.is_open());
    
    linecontent = passComments(file1);
    
    EXPECT_EQ(linecontent[0],'0');
    EXPECT_EQ(linecontent[2],'0');
    EXPECT_EQ(linecontent[4],'0');

    //For the first copolymer.
    for(int32_t x=0;x<16;x++)
        for(int32_t y=0;y<32;y++)
            std::getline(file1,linecontent);
        
    for(int32_t y=0;y<16-9;y++)
        std::getline(file1,linecontent);
        
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,1),"7");
    EXPECT_EQ(linecontent.substr(5,3),"0.5");
    
    //second copolymer
    for(int32_t y=16-9;y<16-3;y++)
        std::getline(file1,linecontent);
        
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"13");
    EXPECT_EQ(linecontent.substr(6,3),"0.5");
    
    //third copolymer
    for(int32_t y=16-3;y<16+3;y++)
        std::getline(file1,linecontent);
    
    EXPECT_EQ(linecontent.substr(0,2),"16");
    EXPECT_EQ(linecontent.substr(3,2),"19");
    EXPECT_EQ(linecontent.substr(6,3),"0.5");
    
    remove("PoreCoordinates.dat");
    remove("PolymerCoordinates.dat");

}    
