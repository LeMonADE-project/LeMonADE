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
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>


/**
 * @file
 * @brief Enable arbitrary walls in the simulation box by the FeatureWall
 * @details For the \b sc-BFM this feature can hold an arbitrary
 * number of walls with normal vectors (1,0,0), (0,1,0) and (0,0,1) being the
 * unit vectors along the principal lattice directions.
 *
 * @todo Enable this feature for the \b bcc-BFM.
 * */

/**
 * @class Wall
 * @brief class providing a single wall with definitions, setter and getter functions
 * */
class Wall
{

public:
    //! standard constructor creating an empty wall
    Wall(): base(), normal() {}

    //! getter function for the base vector of the wall
    const VectorInt3 getBase() const {
        return base;
    }

    //! setter function for the base vector of the wall
    VectorInt3 setBase(uint32_t baseX_, uint32_t baseY_, uint32_t baseZ_) {
        base.setAllCoordinates(baseX_,baseY_,baseZ_);
    }

    //! getter function for the normal vector of the wall
    const VectorInt3 getNormal() const {
        return normal;
    }

    //! setter function for the normal vector of the wall
    VectorInt3 setNormal(uint32_t norX_, uint32_t norY_, uint32_t norZ_) {
      VectorInt3 test(norX_, norY_, norZ_);
      if(test==VectorInt3(1,0,0) || test==VectorInt3(0,1,0) || test==VectorInt3(0,0,1) ){
        normal.setAllCoordinates(norX_,norY_,norZ_);
      }else{
	throw std::runtime_error("normal vector should be P(1,0,0)");
      }
    }

private:

    //! base vector of the wall: arbitrary position somewhere on the wall to provide a well defined wall
    VectorInt3 base;

    //! normal vector of the wall: vector perpendiculat to the walls surface
    VectorInt3 normal;

};


/**
 * @class FeatureWall
 * @brief Feature holding a vector of walls.
 * */
class FeatureWall: public Feature
{
public:

    //! standard constructor
    FeatureWall() {}

    //! standard destructor
    virtual ~FeatureWall(){}

    //! read walls from bfm file
    template<class IngredientsType>
    void exportRead(FileImport<IngredientsType>& fileReader);

    //! write walls to bfm file
    template<class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

    //! check move function for local sc move
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move);

    //! check move function for add sc move
    template<class IngredientsType, class TagType>
    bool checkMove(const IngredientsType& ingredients,MoveAddMonomerSc<TagType>& addmove);

    //! implemantation of synchronize
    template<class IngredientsType>
    void synchronize(const IngredientsType& ingredients);

    //! getter function for the walls container
    std::vector<Wall> getWalls() const{
        return walls;
    }

    /**
     * @brief add function to add a wall to the walls container
     * @throw <std::runtime_error> monomer occupies a position on the walls
     **/
    void addWall(Wall wall){
      if(wall.getNormal().getLength()!=0.0)
        walls.push_back(wall);
      else
	throw std::runtime_error("wall is not well defined: normal vector has length 0");
    }

    //! empty walls container
    void clearAllWalls(){
        walls.clear();
    }


private:
    //! walls container
    std::vector<Wall> walls;

};


/*****************************************************************/
/**
 * @class WriteWall
 *
 * @brief Handles BFM-File-Write \b #!wall
 * @tparam [in] IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteWall: public AbstractWrite<IngredientsType>
{
    public:

        WriteWall(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

        void writeStream(std::ostream& strm){
            const IngredientsType& ingredients=(this->getSource());

            for (size_t i=0; i<ingredients.getWalls().size(); i++) {
                strm << "#!wall: " << ingredients.getWalls()[i].getBase().getX() << " " << ingredients.getWalls()[i].getBase().getY() << " " << ingredients.getWalls()[i].getBase().getZ() << "\t" << ingredients.getWalls()[i].getNormal().getX() << " " << ingredients.getWalls()[i].getNormal().getY() << " " << ingredients.getWalls()[i].getNormal().getZ() << "\n";
            }
            strm << "\n\n";
        }
};


/*****************************************************************/
/**
 * @class ReadWall
 *
 * @brief Handles BFM-File-Read \b #!wall
 * @tparam [in] IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadWall : public ReadToDestination < IngredientsType >
{
    public:

        ReadWall(IngredientsType& destination):ReadToDestination< IngredientsType > (destination){}

        void execute();
};


/**
 * @details checking if move is going to touch one of the walls
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move local sc move
 */
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

/**
 * @details checking if added monomer is going to touch one of the walls
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] addmove move to add sc monomer
 */
template<class IngredientsType, class TagType>
bool FeatureWall::checkMove(const IngredientsType& ingredients, MoveAddMonomerSc<TagType>& addmove)
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

/**
 * @brief Synchronize this feature with the system given as argument
 *
 * @details checking all monomer positions to be not in conflict with one of the walls
 *
 * @throw <std::runtime_error> monomer occupies a position on the walls
 * @param [in] ingredients a reference to the IngredientsType - mainly the system
 **/
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

/**
 * @brief Executes the reading routine to extract \b #!wall.
 *
 * @throw <std::runtime_error> walls could not be read.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
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


/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!wall: bx by bz	nx ny nz
 *
 * @param fileReader File importer for the bfm-file
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template < class IngredientsType >
void FeatureWall::exportRead(FileImport<IngredientsType>& fileReader)
{
    IngredientsType& destination=fileReader.getDestination();
    fileReader.registerRead("#!wall:", new  ReadWall< IngredientsType > (destination));
}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * #!wall: bx by bz	nx ny nz
 *
 * @param fileWriter File writer for the bfm-file.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
void FeatureWall::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
    const IngredientsType& source=fileWriter.getIngredients_();
    fileWriter.registerWrite("#!wall:", new WriteWall <IngredientsType> (source));
}

#endif //LEMONADE_FEATURE_WALL_H
