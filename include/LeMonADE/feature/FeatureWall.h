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

#ifndef LEMONADE_FEATURE_WALL_H
#define LEMONADE_FEATURE_WALL_H

#include <iostream>
#include <sstream>
#include <vector>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddScMonomer.h>
//#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>



class Wall
{
    
public:
    
    Wall(): base(), normal() {}
    
    const VectorInt3 getBase() const {
        return base;
    }
    
    VectorInt3 setBase(uint32_t basX_, uint32_t basY_, uint32_t basZ_) {
        base.setAllCoordinates(basX_,basY_,basZ_);
    }
    
    const VectorInt3 getNormal() const {
        return normal;
    }
    
    VectorInt3 setNormal(uint32_t norX_, uint32_t norY_, uint32_t norZ_) {
        normal.setAllCoordinates(norX_,norY_,norZ_);
    }
    
private:
    
    VectorInt3 base;
    
    VectorInt3 normal;
    
};



class FeatureWall: public Feature
{
public:
    
    //constructor
    FeatureWall() {}
    
    //destructor
    virtual ~FeatureWall(){}
    
    //read and write wall in/from bfm-file
    template<class IngredientsType> 
    void exportRead(FileImport<IngredientsType>& fileReader);
    
    template<class IngredientsType> 
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
    
    //checks if move is allowed by the Metropolis-criterion.
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move);
    
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveAddScMonomer& addmove); 
    
    template<class IngredientsType>
    void synchronize(const IngredientsType& ingredients);
    
    std::vector<Wall> getWalls() const
    {
        return walls;
    }
    
    void addWall(Wall wall)
    {
        walls.push_back(wall);
    }
    
    void clearAllWalls()
    {
        walls.clear();
    }
    
    
private:
    
    std::vector<Wall> walls;
        
};



template <class IngredientsType>
class WriteWall: public AbstractWrite<IngredientsType>
{
    public:
        
        WriteWall(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}
  
        void writeStream(std::ostream& strm)
        {
            const IngredientsType& ingredients=(this->getSource());
            
            for (size_t i=0; i<ingredients.getWalls().size(); i++) {
                strm << "#!wall: " << ingredients.getWalls()[i].getBase().getX() << " " << ingredients.getWalls()[i].getBase().getY() << " " << ingredients.getWalls()[i].getBase().getZ() << "\t" << ingredients.getWalls()[i].getNormal().getX() << " " << ingredients.getWalls()[i].getNormal().getY() << " " << ingredients.getWalls()[i].getNormal().getZ() << "\n";
            }
            strm << "\n\n";
        }
};



template < class IngredientsType >
class ReadWall : public ReadToDestination < IngredientsType >
{
    public:
        
        ReadWall(IngredientsType& destination):ReadToDestination< IngredientsType > (destination){}

        void execute();
};







/* checkMove */

template<class IngredientsType>
bool FeatureWall::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) 
{	
    uint32_t counter(0);
    
	for (size_t i = 0; i < ingredients.getWalls().size(); i++) { //check all walls in the system
        
		uint8_t direction;
        VectorInt3 unitX(1,0,0);
        VectorInt3 unitY(0,1,0);
        VectorInt3 unitZ(0,0,1);
		
		if (ingredients.getWalls()[i].getNormal() == unitX) { direction = 0; } //points direction of normal and sets variable direction to look only on this part of the new position vector later
		if (ingredients.getWalls()[i].getNormal() == unitY) { direction = 1; }
		if (ingredients.getWalls()[i].getNormal() == unitZ) { direction = 2; }

		if ((ingredients.getMolecules()[move.getIndex()] + move.getDir()).getCoordinate(direction) != ingredients.getWalls()[i].getBase().getCoordinate(direction)-1) { //if new position has NOT the same value in "direction" as the wall -1, permit the move
			counter++;
		}
	}
	
	if (counter == ingredients.getWalls().size()) {return true;}
	
	return false;
}



/* checkMove for addmove */

template<class IngredientsType>
bool FeatureWall::checkMove(const IngredientsType& ingredients, MoveAddScMonomer& addmove)
{	
    uint32_t counter(0);
    
	for (size_t i = 0; i < ingredients.getWalls().size(); i++) { //check all walls in the system
        
		uint8_t direction;
        VectorInt3 unitX(1,0,0);
        VectorInt3 unitY(0,1,0);
        VectorInt3 unitZ(0,0,1);
		
		if (ingredients.getWalls()[i].getNormal() == unitX) { direction = 0; } //points direction of normal and sets variable direction to look only on this part of the new position vector later
		if (ingredients.getWalls()[i].getNormal() == unitY) { direction = 1; }
		if (ingredients.getWalls()[i].getNormal() == unitZ) { direction = 2; }

		if (addmove.getPosition().getCoordinate(direction) != ingredients.getWalls()[i].getBase().getCoordinate(direction)-1) { //if new molecule has NOT the same value in "direction" as the wall -1, permit the move
			counter++;
		}
	}
	
	if (counter == ingredients.getWalls().size()) {return true;}
	
	return false;
}





/* synchronize */

template<class IngredientsType>
void FeatureWall::synchronize(const IngredientsType& ingredients)
{
    uint8_t direction;
    VectorInt3 unitX(1,0,0);
    VectorInt3 unitY(0,1,0);
    VectorInt3 unitZ(0,0,1);
        
    for (size_t i=0; i<ingredients.getMolecules().size(); i++) {
		
        for (size_t w = 0; w < ingredients.getWalls().size(); w++) {
                
            if (ingredients.getWalls()[w].getNormal() == unitX) { direction = 0; }
            if (ingredients.getWalls()[w].getNormal() == unitY) { direction = 1; }
            if (ingredients.getWalls()[w].getNormal() == unitZ) { direction = 2; }            
            
            if (ingredients.getMolecules()[i].getCoordinate(direction) == ingredients.getWalls()[w].getBase().getCoordinate(direction)-1) {
                std::ostringstream errorMessage;
                errorMessage << "FeatureWall::synchronize(const IngredientsType& ingredients): Invalid monomer position of monomer " << i << " at " << ingredients.getMolecules()[i] << " in wall: normal: " << ingredients.getWalls()[w].getNormal().getX() << ingredients.getWalls()[w].getNormal().getY() << ingredients.getWalls()[w].getNormal().getZ() << ", base: " << ingredients.getWalls()[w].getBase().getCoordinate(direction)-1 << ".\n";
                throw std::runtime_error(errorMessage.str());
            }
        }
    }
}







/* ReadWall execute() */

template < class IngredientsType >
void ReadWall<IngredientsType>::execute()
{
    std::cout << "reading #!wall...";
    
  
    uint32_t baseX, baseY, baseZ, normalX, normalY, normalZ;
    
    this->getInputStream() >> baseX >> baseY >> baseZ >> normalX >> normalY >> normalZ;
    
    if(!this->getInputStream().fail())
    { 
        Wall wall;
        wall.setBase(baseX,baseY,baseZ);
        wall.setNormal(normalX,normalY,normalZ);
        this->getDestination().addWall(wall);
        
        std::cout << "Base: " << baseX << " " << baseY << " " << baseZ << ", Normal: " << normalX << " " << normalY << " " << normalZ << std::endl;
    }
    else
        throw std::runtime_error("ReadWall<IngredientsType>::execute()\n Could not read wall");
}



/* exportRead wall */

template < class IngredientsType >
void FeatureWall::exportRead(FileImport<IngredientsType>& fileReader)
{
    IngredientsType& destination=fileReader.getDestination();
    fileReader.registerRead("#!wall:", new  ReadWall< IngredientsType > (destination));
}



/* exportWrite wall */

template < class IngredientsType >
void FeatureWall::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
    const IngredientsType& source=fileWriter.getIngredients_();
    fileWriter.registerWrite("#!wall:", new WriteWall <IngredientsType> (source));
}



#endif //LEMONADE_FEATURE_WALL_H
