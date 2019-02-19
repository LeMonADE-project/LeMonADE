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

#ifndef LEMONADE_FEATURE_FEATUREBOX
#define LEMONADE_FEATURE_FEATUREBOX

#include <iostream>
#include <exception>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBoxRead.h>
#include <LeMonADE/feature/FeatureBoxWrite.h>

#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/FileImport.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

/*****************************************************************/
/**
 *
 * @file
 * @class FeatureBox
 *
 * @brief Implementation of a Cartesian simulation box with boundary conditions.
 *
 * @details This feature provides the information about a Cartesian simulation box (Rectangular cuboid).
 * It holds the Length, Width, and Height of the box and also the boundary condition along each axis.
 * This Feature alone never provides a direct lattice etc. This FeatureBox holds general information
 * for the simulation and direct usage are implemented elsewhere.
 * For this see FeatureLatticeBase, FeatureLattice or FeatureLatticePowerOfTwo.
 * It also provides read/write functionality for the commands \b !box_x ,  \b !box_y,  \b !box_z ,
 * \b !periodic_x , \b !periodic_y, \b !periodic_z.
 **/
/*****************************************************************/
class FeatureBox : public Feature
{

public:

	//! Default constructor. Set Length=Width=Height=0 and P.B.C. as false (hard walls)
	FeatureBox();

	//! Default destructor (empty)
	virtual ~FeatureBox() {};

	//! Returns the number of all possible lattice entries - the volume of the box.
	uint64_t getNumberOfLatticeSites() const ;

	/**
	 * @brief Synchronize this feature with the system given as argument
	 *
	 * @details Synchronize this FeatureBox with the system given as argument.
	 * It check if the size of the simulationbox and periodic boundary is set.
	 *
	 * @param ing a reference to the IngredientsType - mainly the system
	 **/
	template < class IngredientsType > void synchronize(IngredientsType& ing) {
		assertBoxSizeSet();
		assertPeriodicitySet();
		assertParticlesInBox(ing.getMolecules());
	};

	/**
	 * @brief Set the Cartesian x-component of the box (Length).
	 *
	 * @param x The new length of the box.
	 */
	void setBoxX(int32_t x);

	/**
	 * @brief Set the Cartesian y-component of the box (Width).
	 *
	 * @param y The new width of the box.
	 */
	void setBoxY(int32_t y);

	/**
	 * @brief Set the Cartesian z-component of the box (Height).
	 *
	 * @param z The height of the box.
	 */
	void setBoxZ(int32_t z);

	/**
	 * @brief Get the Cartesian x-component of the box (Length).
	 *
	 * @return The length of the box.
	 */
	int32_t  getBoxX() const ;

	/**
	 * @brief Get the Cartesian y-component of the box (Width).
	 *
	 * @return The width of the box.
	 */
	int32_t  getBoxY() const ;

	/**
	 * @brief Get the Cartesian z-component of the box (Height).
	 *
	 * @return The height of the box.
	 */
	int32_t  getBoxZ() const ;

	/**
	 * @brief Set the periodic boundary condition in x-direction of the box (Length).
	 *
	 * @param x Periodic boundary condition on (true). Otherwise false.
	 */
	void setPeriodicX(bool x);

	/**
	 * @brief Set the periodic boundary condition in y-direction of the box (Width).
	 *
	 * @param y Periodic boundary condition on (true). Otherwise false.
	 */
	void setPeriodicY(bool y);

	/**
	 * @brief Set the periodic boundary condition in z-direction of the box (Height).
	 *
	 * @param z Periodic boundary condition on (true). Otherwise false.
	 */
	void setPeriodicZ(bool z);

	/**
	 * @brief Get the periodic boundary condition in x-direction of the box (Length).
	 *
	 * @return True if Periodic boundary condition is set on. Otherwise false.
	 */
	bool  isPeriodicX() const ;

	/**
	 * @brief Get the periodic boundary condition in y-direction of the box (Width).
	 *
	 * @return True if Periodic boundary condition is set on. Otherwise false.
	 */
	bool  isPeriodicY() const ;

	/**
	 * @brief Get the periodic boundary condition in z-direction of the box (Height).
	 *
	 * @return True if Periodic boundary condition is set on. Otherwise false.
	 */
	bool  isPeriodicZ() const ;

	/**
	 * @brief Returns if the rectangular cuboid is a simple cube with  all length, height, width are all the same.
	 *
	 * @return True if rectangular cuboid is a cube. Otherwise false.
	 */
	bool isCubic() const;


	/**
	 * @brief Export the relevant functionality for reading bfm-files to the responsible reader object
	 *
	 * @details The function is called by the Ingredients class when an object of type Ingredients
	 * is associated with an object of type FileImport. The export of the Reads is thus
	 * taken care automatically when it becomes necessary.\n
	 * Registered Read-In Commands:
	 * * !box_x
	 * * !box_y
	 * * !box_z
	 * * !periodic_x
	 * * !periodic_y
	 * * !periodic_z
	 *
	 * @param fileReader File importer for the bfm-file
	 * @tparam IngredientsType Features used in the system. See Ingredients.
	 **/
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader)
	{
	//  ModelFeature provides an empty exportRead(), which is overwritten here.

	fileReader.registerRead("!box_x",new ReadBoxX <FeatureBox> (*this));
	fileReader.registerRead("!box_y",new ReadBoxY <FeatureBox> (*this));
	fileReader.registerRead("!box_z",new ReadBoxZ <FeatureBox> (*this));
	fileReader.registerRead("!periodic_x",new ReadPeriodicX <FeatureBox> (*this));
	fileReader.registerRead("!periodic_y",new ReadPeriodicY <FeatureBox> (*this));
	fileReader.registerRead("!periodic_z",new ReadPeriodicZ <FeatureBox> (*this));

	// EmptyModelFeature is used as the "dead end" of a the genareted class using Loki-typelist
	// See  GenerateContextType in GenericMonteCarlo.h.
	// The call of other independent Feature::exportReads can be handled by using an actual Feature "Holder"
	// such as ModelFeatureHolder GenerateContextType ... , which triggers function calls
	// of the Features and their precedors.

	};

	/**
	 * @brief Export the relevant functionality for writing bfm-files to the responsible writer object
	 *
	 * @details The function is called by the Ingredients class when an object of type Ingredients
	 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
	 * taken care automatically when it becomes necessary.\n
	 * Registered Write-Out Commands:
	 * * !box_x
	 * * !box_y
	 * * !box_z
	 * * !periodic_x
	 * * !periodic_y
	 * * !periodic_z
	 *
	 * @param fileWriter File writer for the bfm-file.
	 **/
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter)const {
	fileWriter.registerWrite("!box_x",new WriteBoxX<FeatureBox>(*this));
	fileWriter.registerWrite("!box_y",new WriteBoxY<FeatureBox>(*this));
	fileWriter.registerWrite("!box_z",new WriteBoxZ<FeatureBox>(*this));
	fileWriter.registerWrite("!periodic_x",new WritePeriodicX<FeatureBox>(*this));
	fileWriter.registerWrite("!periodic_y",new WritePeriodicY<FeatureBox>(*this));
	fileWriter.registerWrite("!periodic_z",new WritePeriodicZ<FeatureBox>(*this));
	}

	/**
	 * @brief For all unknown moves: this does nothing
	 *
	 * @details Returns true for all moves other than the ones that have specialized versions of this function.
	 * This dummy function is implemented for generality.
	 *
	 * @param [in] ingredients A reference to the IngredientsType - mainly the system
	 * @param [in] move General move other than MoveLocalBase.
	 * @return true Always!
	 */
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const
	{
		//accept unknown moves
		return true;
	}


	//for moves of type MoveLocalBase check periodicity

	/**
	 * @brief Overloaded for MoveLocalBase (MoveLocalSc or BcccLocalMove)
	 *
	 * @details Returns true if the moves doesn´t violate the p.b.c. The periodicity is set to \a false ,
	 * the faces of rectangular cuboid behave like hard walls. In fact, the new position \a pos of the monomer can not exceed
	 * the limit in Length [pos.getX()<(getBoxX()-1)], Width [pos.getY()<(getBoxY()-1)] or Height [pos.getZ()<(getBoxZ()-1)]
	 * or be smaller than 0 in the corresponding directions. \n
	 * If the periodicity is set to \a true , the Move is not limited to the (virtual) simulation box.
	 *
	 * @param [in] ingredients A reference to the IngredientsType - mainly the system
	 * @param [in] move General move other than MoveLocalBase.
	 * @return True if move is allowed with p.b.c. or rejected (false).
	 *
	 * @todo This implementation of "walls" maybe slow down the algorithm. Discuss a new FeatureWall
	 */
	template<class IngredientsType,class LocalMoveType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBase<LocalMoveType>& move) const
	{
		typename IngredientsType::molecules_type::vertex_type pos=ingredients.getMolecules()[move.getIndex()];
		pos+=move.getDir();
		if(	(periodicX || (pos.getX()<(getBoxX()-1) && pos.getX()>=0) ) &&
			(periodicY || (pos.getY()<(getBoxY()-1) && pos.getY()>=0) ) &&
			(periodicZ || (pos.getZ()<(getBoxZ()-1) && pos.getZ()>=0) )
		)
			return true;
		else
			return false;

	}

	/**
	 * @brief Overloaded for MoveConnectSc
	 *
	 * @details Returns true if the moves doesn´t violate the p.b.c. esp. the new link is not created of walls.
	 * The periodicity is set to \a false ,
	 * the faces of rectangular cuboid behave like hard walls. In fact, the new position \a pos of the monomer can not exceed
	 * the limit in Length [pos.getX()<(getBoxX()-1)], Width [pos.getY()<(getBoxY()-1)] or Height [pos.getZ()<(getBoxZ()-1)]
	 * or be smaller than 0 in the corresponding directions. \n
	 * If the periodicity is set to \a true , the Move is not limited to the (virtual) simulation box.
	 *
	 * @param [in] ingredients A reference to the IngredientsType - mainly the system
	 * @param [in] move General move other than MoveConnectSc.
	 * @return True if move is allowed with p.b.c. or rejected (false).
	 *
	 * @todo This implementation of "walls" maybe slow down the algorithm. Discuss a new FeatureWall
	 */
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectSc& move) const
	{
		typename IngredientsType::molecules_type::vertex_type targetedPos=ingredients.getMolecules()[move.getIndex()];
		targetedPos+=move.getDir(); // // pos += P+-(2,0,0)
		if(	(periodicX || (targetedPos.getX()<(getBoxX()-1) && targetedPos.getX()>=0) ) &&
			(periodicY || (targetedPos.getY()<(getBoxY()-1) && targetedPos.getY()>=0) ) &&
			(periodicZ || (targetedPos.getZ()<(getBoxZ()-1) && targetedPos.getZ()>=0) )
		)
			return true;
		else
			return false;

	}

	/**
	 * @brief Overloaded for MoveAddMonomerBase(MoveAddMonomerSc and MoveAddMonomerBcc)
	 *
	 * @details Returns true if the moves doesn´t violate the p.b.c. The periodicity is set to \a false ,
	 * the faces of rectangular cuboid behave like hard walls. In fact, the new position \a pos of the monomer can not exceed
	 * the limit in Length [pos.getX()<(getBoxX()-1)], Width [pos.getY()<(getBoxY()-1)] or Height [pos.getZ()<(getBoxZ()-1)]
	 * or be smaller than 0 in the corresponding directions. \n
	 * If the periodicity is set to \a true , the Move is not limited to the (virtual) simulation box.
	 *
	 * @param [in] ingredients A reference to the IngredientsType - mainly the system
	 * @param [in] addmove Move of type MoveAddMonomerSc or MoveAddMonomerBcc
	 * @return True if move is allowed with p.b.c. or rejected (false).
	 *
	 */
	template<class IngredientsType,class AddMoveType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBase<AddMoveType, TagType>& addmove) const
	{
		VectorInt3 pos=addmove.getPosition();

		if(	(periodicX || (pos.getX()<(getBoxX()-1) && pos.getX()>=0) ) &&
			(periodicY || (pos.getY()<(getBoxY()-1) && pos.getY()>=0) ) &&
			(periodicZ || (pos.getZ()<(getBoxZ()-1) && pos.getZ()>=0) )
		)
			return true;
		else
			return false;

	}


private:

  //! Function for asserting then variables \a boxX (Length) has not been set.
  void assertBoxSizeSetX() const ;

  //! Function for asserting then variables \a boxY (Width) has not been set.
  void assertBoxSizeSetY() const ;

  //! Function for asserting then variables \a boxZ (Height) has not been set.
  void assertBoxSizeSetZ() const ;

  //! Function for asserting then variables \a boxX (Length), \a boxY (Width), \a boxZ (Height) has not been set.
  void assertBoxSizeSet() const ;

  //! Function for asserting then variables \a periodicX (p.b.c in x-direction) has not been set.
  void assertPeriodicitySetX() const ;

  //! Function for asserting then variables \a periodicY (p.b.c in y-direction) has not been set.
  void assertPeriodicitySetY() const ;

  //! Function for asserting then variables \a periodicZ (p.b.c in z-direction) has not been set.
  void assertPeriodicitySetZ() const ;

  //! Function for asserting then variables \a periodicX, \a periodicY, and \a periodicZ has not been set.
  void assertPeriodicitySet() const ;

  template<class MoleculesType>
  void assertParticlesInBox(const MoleculesType& molecules);

  //! Holding the information about the Length of simulation box (rectangular cuboid)
  int32_t boxX;

  //! Holding the information about the Width of simulation box (rectangular cuboid)
  int32_t boxY;

  //! Holding the information about the Height of simulation box (rectangular cuboid)
  int32_t boxZ;

  //! Holding the information about the p.b.c in x-direction
  bool periodicX;

  //! Holding the information about the p.b.c in y-direction
  bool periodicY;

  //! Holding the information about the p.b.c in z-direction
  bool periodicZ;

  //! Holding the information if periodicity in x-direction was initialized (false-not set; true-set)
  bool periodicInitX;

  //! Holding the information if periodicity in y-direction was initialized (false-not set; true-set)
  bool periodicInitY;

  //! Holding the information if periodicity in z-direction was initialized (false-not set; true-set)
  bool periodicInitZ;


};

/**
 * @brief This function checks if monomers are inside the simulation box if p.b.c are not present
 *
 * @throw <std::runtime_error> if p.b.c are false and  monomer position are "outside" of the box
 *
 * @param molecules Reference to all Molecules in the system.
 */
template<class MoleculesType>
void FeatureBox::assertParticlesInBox(const MoleculesType& molecules)
{

	if( (!isPeriodicX()) || (!isPeriodicY()) || (!isPeriodicZ()) )
	{
		for (size_t n=0;n<molecules.size();++n)
		{
			if ( ( (!periodicX) && (molecules[n].getX()>=(boxX-1) || molecules[n].getX()<0) )||
				 ( (!periodicY) && (molecules[n].getY()>=(boxY-1) || molecules[n].getY()<0) )||
				 ( (!periodicZ) && (molecules[n].getZ()>=(boxZ-1) || molecules[n].getZ()<0) )
			)
			{
				std::stringstream errormessage;
				errormessage<<"**** FeatureBox::synchronize(): particle position of monomer "<<n<< " ("<< molecules[n].getX() <<"," << molecules[n].getY() << "," << molecules[n].getZ() << ") " <<" outside box ***\n";
				throw std::runtime_error(errormessage.str());
			}
		}
	}

}
#endif
