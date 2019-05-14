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

#ifndef LEMONADE_FEATURE_FEATURELATTICEBASE_H
#define LEMONADE_FEATURE_FEATURELATTICEBASE_H

#include <iostream>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/utility/Vector3D.h>

/**
 * @file FeatureLatticeBase.h
 * @date 2014/07/14
 * @author Ron
 *
 **/

/**
 * @brief forward declaration for FeatureLatticeBase
 */
template<typename SpecializedClass> class FeatureLatticeBase;


/**
 * @brief Provides lattice with variable value type. This is the primitive base class for the specialized
 *        versions FeatureLattice and FeatureLatticePowerOfTwo using Curiously Recurring Template Pattern (CRTP) to avoid virtual functions.
 *
 * @details All lattices must be derived from this base class in order to be processed
 * 			correctly by the Features. All FeatureLattice should implement the functions foldBackX, foldBackY and foldBackZ for folding to relative coordinates.
 * 			Since this class simply serves as a common base, the functions here don't do
 * 			anything particular.
 *
 * @tparam <ValueType> type of the lattice value, defaults to bool.
 * @tparam <SpecializedClass> name of the specialized class.
 **/
template<template<typename> class SpecializedClass, typename ValueType>
class FeatureLatticeBase< SpecializedClass<ValueType> >:public Feature
{
public:
	//! This Feature requires a box.
	typedef LOKI_TYPELIST_1(FeatureBox) required_features_front;

	FeatureLatticeBase();
	virtual ~FeatureLatticeBase();

	FeatureLatticeBase(const FeatureLatticeBase& copyFeatureLatticeBase);

 	FeatureLatticeBase& operator= (const FeatureLatticeBase &FeatureLatticeBaseSource);


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

	//! Synchronize with system
	template<class IngredientsType> void synchronize(IngredientsType& val);

	//! Allocate memory for the lattice
	void setupLattice();

	//! Allocate memory for the lattice
	void setupLattice(uint32_t boxX, uint32_t boxY,uint32_t boxZ );


	//! Set lattice entries to 0
	void clearLattice();

	//delocate memory
        void deleteLattice();

protected:

	//! Hold the value of lattice size in X
	uint32_t _boxX;

	//! Hold the value of lattice size in Y
	uint32_t _boxY;

	//! Hold the value of lattice size in Z
	uint32_t _boxZ;

	//! Hold the value of lattice size in X subtracting one for folding
	uint32_t boxXm1;

	//! Hold the value of lattice size in Y subtracting one for folding
	uint32_t boxYm1;

	//! Hold the value of lattice size in Z subtracting one for folding
	uint32_t boxZm1;

	//! Hold the value of for indexing the lattice (FeatureLattice: boxX; FeatureLatticePowerOfTwo: log(boxX,2))
	uint32_t xPro;

	//! Hold the value of for indexing the lattice (FeatureLattice: boxX*boxY; FeatureLatticePowerOfTwo: log(boxX,2)*log(boxY,2))
	uint32_t proXY;

	/**
	 * @brief Linearized 3D lattice of type ValueType to 1D.
	 *
	 * @details Linearized 3D lattice of type \a ValueType to 1D by formula in:
	 * FeatureLattice: lattice[idx]=lattice[x+y*xPro+z*proXY]
	 * FeatureLatticePowerOfTwo: lattice[idx]=lattice[x+(y<<xPro)+(z<<proXY)]
	 */
	ValueType* lattice;
};

/******************************************************************************/
/***************************definition of members******************************/

/******************************************************************************/
//!constructor
template<template<typename> class SpecializedClass, typename ValueType>
FeatureLatticeBase<SpecializedClass<ValueType> >::FeatureLatticeBase()
	:_boxX(0),_boxY(0),_boxZ(0),boxXm1(0),boxYm1(0),boxZm1(0),xPro(0),proXY(0),lattice(NULL)
{

}

/******************************************************************************/
//!destructor
template<template<typename> class SpecializedClass, typename ValueType>
FeatureLatticeBase<SpecializedClass<ValueType> >::~FeatureLatticeBase()
{
    this->deleteLattice();
	// free memory
// 	if(lattice)
// 		delete[] lattice;
}

//!destructor
template<template<typename> class SpecializedClass, typename ValueType>
void FeatureLatticeBase<SpecializedClass<ValueType> >::deleteLattice()
{
    // free memory
    if(lattice != NULL) delete[] lattice;
    lattice = NULL;
}


/**
 * @todo testing!!!
 */
template<template<typename> class SpecializedClass, typename ValueType>
FeatureLatticeBase<SpecializedClass<ValueType> >::FeatureLatticeBase(const FeatureLatticeBase<SpecializedClass<ValueType> >& copyFeatureLatticeBase){

	std::cout << "copyFeatureLatticeBase" << std::endl;

	_boxX = copyFeatureLatticeBase._boxX;
	_boxY = copyFeatureLatticeBase._boxY;
	_boxZ = copyFeatureLatticeBase._boxZ;

	boxXm1 = copyFeatureLatticeBase.boxXm1;
	boxYm1 = copyFeatureLatticeBase.boxYm1;
	boxZm1 = copyFeatureLatticeBase.boxZm1;

	xPro = copyFeatureLatticeBase.xPro;
	proXY = copyFeatureLatticeBase.proXY;

	lattice = new ValueType[_boxX*_boxY*_boxZ];

}


/**
 * @todo testing!!!
 */
template<template<typename> class SpecializedClass, typename ValueType>
FeatureLatticeBase<SpecializedClass<ValueType> >&  FeatureLatticeBase<SpecializedClass<ValueType> >::operator = (const FeatureLatticeBase<SpecializedClass<ValueType> > &FeatureLatticeBaseSource)
{

    // check for self-assignment by comparing the address of the
    // implicit object and the parameter
    if (this == &FeatureLatticeBaseSource)
        return *this;

    uint64_t oldSize = _boxX*_boxY*_boxZ;
    uint64_t newSize = FeatureLatticeBaseSource._boxX*FeatureLatticeBaseSource._boxY*FeatureLatticeBaseSource._boxZ;

    _boxX  = FeatureLatticeBaseSource._boxX;
    _boxY  = FeatureLatticeBaseSource._boxY;
    _boxZ  = FeatureLatticeBaseSource._boxZ;

    boxXm1 = FeatureLatticeBaseSource.boxXm1;
    boxYm1 = FeatureLatticeBaseSource.boxYm1;
    boxZm1 = FeatureLatticeBaseSource.boxZm1;

    xPro   = FeatureLatticeBaseSource.xPro;
    proXY  = FeatureLatticeBaseSource.proXY;

    if ( oldSize != newSize )
    {
        this->deleteLattice();
        lattice = new ValueType[_boxX*_boxY*_boxZ];
    }

    // do the copy
    // return the existing object
    return *this;
}

/******************************************************************************/
/**
 * @details Move a lattice point to a new position. The data at the new position
 * is overwritten and the old position is set to ValueType() (most Zero).
 *
 * @param oldPos old position
 * @param newPos new position
 */
/******************************************************************************/
template<template<typename> class SpecializedClass, typename ValueType>
inline void FeatureLatticeBase<SpecializedClass<ValueType> >::moveOnLattice(const VectorInt3& oldPos, const VectorInt3& newPos)
{
	static_cast<SpecializedClass<ValueType>* >(this)->moveOnLattice(oldPos, newPos);
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

template<template<typename> class SpecializedClass, typename ValueType>
inline void FeatureLatticeBase<SpecializedClass<ValueType> >::moveOnLattice(const int xOldPos, const int yOldPos, const int zOldPos, const int xNewPos, const int yNewPos, const int zNewPos)
{
	static_cast<SpecializedClass<ValueType>* >(this)->moveOnLattice(xOldPos, yOldPos, zOldPos, xNewPos, yNewPos, zNewPos);
}

/**
 * Get the value stored on the lattice at coordinates given by VectorInt3 \a pos.
 * @param[in] pos specified position
 * @return \p ValueType value on the specified position \a pos
 */
template<template<typename> class SpecializedClass, typename ValueType>
inline ValueType FeatureLatticeBase<SpecializedClass<ValueType> >::getLatticeEntry(const VectorInt3& pos) const
{
	return static_cast<const SpecializedClass<ValueType>* >(this)->getLatticeEntry(pos);
}


/**
 * Get the value stored on the lattice at coordinates given by Cartesian x y z coordinates.
 * @param[in] x x-coordinate on the Cartesian lattice
 * @param[in] y y-coordinate on the Cartesian lattice
 * @param[in] z z-coordinate on the Cartesian lattice
 * @return \a ValueType value on the specified position at x y z
 */
template<template<typename> class SpecializedClass, typename ValueType>
inline ValueType FeatureLatticeBase<SpecializedClass<ValueType> >::getLatticeEntry(const int x, const int y, const int z) const
{
	return static_cast<const SpecializedClass<ValueType>* >(this)->getLatticeEntry(x, y, z);
}


/**
 * Set the value \a val on the lattice at coordinates given by \a pos.
 * @param[in] pos specified position
 * @param[in] val \e ValueType to set on \a pos
 */
template<template<typename> class SpecializedClass, typename ValueType>
inline void FeatureLatticeBase<SpecializedClass<ValueType> >::setLatticeEntry(const VectorInt3& pos, ValueType val)
{
	static_cast<SpecializedClass<ValueType>* >(this)->setLatticeEntry(pos, val);
}


/**
 * Set the value \a val on the lattice at coordinates given by \a pos.
 * @param[in] x x-coordinate on the Cartesian lattice
 * @param[in] y y-coordinate on the Cartesian lattice
 * @param[in] z z-coordinate on the Cartesian lattice
 * @param[in] val \e ValueType to set on \a pos
 */
template<template<typename> class SpecializedClass, typename ValueType>
inline void FeatureLatticeBase<SpecializedClass<ValueType> >::setLatticeEntry(const int x, const int y, const int z, ValueType val)
{
	static_cast<SpecializedClass<ValueType>* >(this)->setLatticeEntry(x, y, z, val);
}


/**
 * Synchronize this feature with the system given as argument. Creates and recreates the lattice.
 * This method set all lattice entries as value-initialized \a ValueType() (most cases Zero).
 * It does \a not populate the lattice.
 * The synchronization is done by the derived Feature class.
 *
 * @param ing a reference to the IngredientsType - mainly the system
 **/
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
void FeatureLatticeBase<SpecializedClass<ValueType> >::synchronize(IngredientsType& val)
{
	static_cast<const SpecializedClass<ValueType>* >(this)->synchronize(val);
}



/**
 * This method allocates memory for the lattice and fill it with the native value of \a ValueType.
 */
template<template<typename> class SpecializedClass, typename ValueType>
void FeatureLatticeBase<SpecializedClass<ValueType> >::setupLattice()
{

	std::cout<<"setting up lattice...";

	// Allocate memory
	lattice = new ValueType[_boxX*_boxY*_boxZ];

	for(uint32_t i = 0; i < _boxX*_boxY*_boxZ; i++)
		lattice[i]=ValueType();

	std::cout<<"done with size " << (_boxX*_boxY*_boxZ*sizeof(ValueType)) << " bytes = " << (_boxX*_boxY*_boxZ*sizeof(ValueType)/(1024.0*1024.0)) << " MB for lattice" <<std::endl;

}

/**
 * This method allocates memory for the lattice and fill it with the native value of \a ValueType.
 */
template<template<typename> class SpecializedClass, typename ValueType>
void FeatureLatticeBase<SpecializedClass<ValueType> >::setupLattice(uint32_t boxX, uint32_t boxY,uint32_t boxZ )
{

	_boxX = boxX;
	_boxY = boxY;
	_boxZ = boxZ;

	this->boxXm1=this->_boxX-1;
	this->boxYm1=this->_boxY-1;
	this->boxZm1=this->_boxZ-1;

	// determine the shift values for first multiplication
	this->xPro=this->_boxX;

	// determine the shift values for second multiplication
	this->proXY=this->_boxX*this->_boxY;

	std::cout << "use bit shift for boxX: ("<< this->xPro << " ) = " << (this->xPro) << " = " << (this->_boxX) << std::endl;
	std::cout << "use bit shift for boxX*boxY: ("<< this->proXY << " ) = " << (this->proXY) << " = " << (this->_boxX*this->_boxY) << std::endl;


	std::cout<<"setting up lattice...";

	// Allocate memory
	lattice = new ValueType[_boxX*_boxY*_boxZ];

	for(uint32_t i = 0; i < _boxX*_boxY*_boxZ; i++)
		lattice[i]=ValueType();

	std::cout<<"done with size " << (_boxX*_boxY*_boxZ*sizeof(ValueType)) << " bytes = " << (_boxX*_boxY*_boxZ*sizeof(ValueType)/(1024.0*1024.0)) << " MB for lattice" <<std::endl;

}


/**
 * All lattice entries are set to the native value of \a ValueType (in most cases Zero).
 * This method only clears the lattice. It will not destroy nor recreate the array.
 */
template<template<typename> class SpecializedClass, typename ValueType>
void FeatureLatticeBase<SpecializedClass<ValueType> >::clearLattice()
{

	// initialize to native value (=0)
	for(uint32_t i = 0; i < _boxX*_boxY*_boxZ; i++)
		lattice[i]=ValueType();
}

#endif /* LEMONADE_FEATURE_FEATURELATTICEBASE_H */
