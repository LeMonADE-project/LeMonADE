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

#ifndef LEMONADE_FEATURE_FEATURELATTICE_H
#define LEMONADE_FEATURE_FEATURELATTICE_H

#include <LeMonADE/feature/FeatureLatticeBase.h>

/*****************************************************************************/
/**
 * @file
 * @date 2014/07/14
 * @author Ron
 *
 * @class FeatureLattice
 * @brief Provides lattice with variable value type for general lattice sizes.
 *        If you have 2^n lattice in all direction use the optimized version FeatureLatticePowerOfTwo .
 *
 * @tparam ValueType type of the lattice value, default is \a bool.
 * */
/*****************************************************************************/
template< class ValueType=bool>
class FeatureLattice: public FeatureLatticeBase< FeatureLattice<ValueType> > {
public:

	FeatureLattice(){};
	virtual ~FeatureLattice(){};

	//! Move the value on the lattice to a new position. Delete the Value on the old position
	void moveOnLattice(const VectorInt3& oldPos, const VectorInt3& newPos);

	//! Move the value on the lattice to a new position. Delete the Value on the old position
	void moveOnLattice(const int xOldPos, const int yOldPos, const int zOldPos, const int xNewPos, const int yNewPos, const int zNewPos);

	//! Get the lattice value at a certain point
	ValueType getLatticeEntry(const VectorInt3& pos) const;

	//! Get the lattice value at a certain point
	ValueType getLatticeEntry(const int x, const int y, const int z) const;

	//! Set the value on a lattice point
	void setLatticeEntry(const VectorInt3& pos, ValueType val);

	//! Set the value on a lattice point
	void setLatticeEntry(const int x, const int y, const int z, ValueType val);

	//! Synchronize this feature with the system given as argument
	template<class IngredientsType> void synchronize(IngredientsType& ing);

private:
	//! Functions for folding absolute coordinates into the lattice in X
	uint32_t foldBackX(int value) const;

	//! Functions for folding absolute coordinates into the lattice in Y
	uint32_t foldBackY(int value) const;

	//! Functions for folding absolute coordinates into the lattice in Z
	uint32_t foldBackZ(int value) const;
};

/******************************************************************************/
/***************************definition of members******************************/

/******************************************************************************/
/**
 * @details Move a lattice point to a new position. The data at the new position
 * is overwritten and the old position is set to ValueType() (most Zero).
 *
 * @param oldPos old position
 * @param newPos new position
 */
/******************************************************************************/
template<class ValueType>
inline void FeatureLattice<ValueType>::moveOnLattice(const VectorInt3& oldPos, const VectorInt3& newPos)
{
	uint32_t xOld,yOld,zOld;
	xOld=foldBackX(oldPos[0]);
	yOld=foldBackY(oldPos[1]);
	zOld=foldBackZ(oldPos[2]);

	this->lattice[foldBackX(newPos[0])+(foldBackY(newPos[1]) * this->xPro)+(foldBackZ(newPos[2])*this->proXY)] = this->lattice[xOld+(yOld * this->xPro)+(zOld*this->proXY)];
	this->lattice[xOld+(yOld * this->xPro)+(zOld*this->proXY)] = ValueType();

}

/******************************************************************************/
/**
 * @details Move a lattice point to a new position. The data at the new position
 * is overwritten and the old position is set to ValueType() (most Zero).
 *
 * @param[in] xOldPos x-coordinate of old position in absolute coordinates
 * @param[in] yOldPos y-coordinate of old position in absolute coordinates
 * @param[in] zOldPos z-coordinate of old position in absolute coordinates
 * @param[in] xNewPos x-coordinate of new position in absolute coordinates
 * @param[in] yNewPos y-coordinate of new position in absolute coordinates
 * @param[in] zNewPos z-coordinate of new position in absolute coordinates
 */
/******************************************************************************/
template<class ValueType>
inline void FeatureLattice<ValueType>::moveOnLattice(const int xOldPos, const int yOldPos, const int zOldPos, const int xNewPos, const int yNewPos, const int zNewPos)
{
	uint32_t xOld,yOld,zOld;
	xOld=foldBackX(xOldPos);
	yOld=foldBackY(yOldPos);
	zOld=foldBackZ(zOldPos);

	this->lattice[foldBackX(xNewPos)+(foldBackY(yNewPos) * this->xPro)+(foldBackZ(zNewPos)*this->proXY)] = this->lattice[xOld+(yOld * this->xPro)+(zOld*this->proXY)];
	this->lattice[xOld+(yOld * this->xPro)+(zOld*this->proXY)] = ValueType();

}

/**
 * Get the value stored on the lattice at coordinates given by VectorInt3 \a pos.
 * @param[in] pos specified position
 * @return \p ValueType value on the specified position \a pos
 */
template<class ValueType>
inline ValueType FeatureLattice<ValueType>::getLatticeEntry(const VectorInt3& pos) const
{
	return (this->lattice[foldBackX(pos[0])+(foldBackY(pos[1])* this->xPro)+(foldBackZ(pos[2])*this->proXY)]);
}


/**
 * Get the value stored on the lattice at coordinates given by Cartesian x y z coordinates.
 * @param[in] x x-coordinate on the Cartesian lattice
 * @param[in] y y-coordinate on the Cartesian lattice
 * @param[in] z z-coordinate on the Cartesian lattice
 * @return \a ValueType value on the specified position at x y z
 */
template<class ValueType>
inline ValueType FeatureLattice<ValueType>::getLatticeEntry(const int x, const int y, const int z) const
{
	return(this->lattice[foldBackX(x)+(foldBackY(y)*this->xPro)+(foldBackZ(z)*this->proXY)]);
}


/**
 * Set the value \a val on the lattice at coordinates given by \a pos.
 * @param[in] pos specified position
 * @param[in] val \e ValueType to set on \a pos
 */
template<class ValueType>
inline void FeatureLattice<ValueType>::setLatticeEntry(const VectorInt3& pos, ValueType val)
{
	this->lattice[foldBackX(pos[0])+(foldBackY(pos[1])* this->xPro)+(foldBackZ(pos[2])*this->proXY)]=val;
}


/**
 * Set the value \a val on the lattice at coordinates given by \a pos.
 * @param[in] x x-coordinate on the Cartesian lattice
 * @param[in] y y-coordinate on the Cartesian lattice
 * @param[in] z z-coordinate on the Cartesian lattice
 * @param[in] val \e ValueType to set on \a pos
 */
template<class ValueType>
inline void FeatureLattice<ValueType>::setLatticeEntry(const int x, const int y, const int z, ValueType val)
{
	this->lattice[foldBackX(x)+(foldBackY(y)* this->xPro)+(foldBackZ(z)*this->proXY)]=val;
}


/**
 * Fold back the absolute coordinate into the relative coordinate in X by modulo operation.
 *
 * @param value absolute coordinate in X
 * @return \a uint32_t relative coordinate in X
 */
template<class ValueType>
inline uint32_t FeatureLattice<ValueType>::foldBackX(int value) const{
	return (((value%this->_boxX)+this->_boxX)%this->_boxX);
}

/**
 * Fold back the absolute coordinate into the relative coordinate in Y by modulo operation.
 *
 * @param value absolute coordinate in Y
 * @return \a uint32_t relative coordinate in Y
 */
template<class ValueType>
inline uint32_t FeatureLattice<ValueType>::foldBackY(int value) const{
	return (((value%this->_boxY)+this->_boxY)%this->_boxY);
}

/**
 * Fold back the absolute coordinate into the relative coordinate in Z by modulo operation.
 *
 * @param value absolute coordinate in Z
 * @return \a uint32_t relative coordinate in Z
 */
template<class ValueType>
inline uint32_t FeatureLattice<ValueType>::foldBackZ(int value) const{
	return (((value%this->_boxZ)+this->_boxZ)%this->_boxZ);
}


/**
 * Synchronize this feature with the system given as argument. Creates and recreates the lattice.
 * This method set all lattice entries as value-initialized \a ValueType() (most cases Zero). It does \a not populate the lattice.
 * The synchronization is valid for general lattice sizes in all directions.
 *
 * @param ing a reference to the IngredientsType - mainly the system
 **/
template<class ValueType>
template<class IngredientsType>
void FeatureLattice<ValueType>::synchronize(IngredientsType& ing) {

	//if the lattice is already initialized, free the memory first
		if(this->lattice)
			delete[] this->lattice;


		this->_boxX=ing.getBoxX();
		this->_boxY=ing.getBoxY();
		this->_boxZ=ing.getBoxZ();

		this->boxXm1=this->_boxX-1;
		this->boxYm1=this->_boxY-1;
		this->boxZm1=this->_boxZ-1;

		// determine the shift values for first multiplication
		this->xPro=this->_boxX;

		// determine the shift values for second multiplication
		this->proXY=this->_boxX*this->_boxY;

		std::cout << "use bit shift for boxX: ("<< this->xPro << " ) = " << (this->xPro) << " = " << (this->_boxX) << std::endl;
		std::cout << "use bit shift for boxX*boxY: ("<< this->proXY << " ) = " << (this->proXY) << " = " << (this->_boxX*this->_boxY) << std::endl;

		// check if shift is correct
		if ( (this->_boxX != (this->xPro)) || ((this->_boxX*this->_boxY) != (this->proXY)) )
		{
			throw  std::runtime_error("Could not determine value for lattice indexing.\n");
		}

		//allocate memory
		this->setupLattice();
		//set all values to 0
		this->clearLattice();
}

#endif /* LEMONADE_FEATURE_FEATURELATTICE_H */
