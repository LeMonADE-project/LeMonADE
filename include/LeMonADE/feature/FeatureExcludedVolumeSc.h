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

#ifndef LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMESC_H
#define LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMESC_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>

/*****************************************************************************/
/**
 * @file
 * @date   2016/02/18
 * @author Hauke/Toni
 *
 * @class FeatureExcludedVolumeSc
 * @brief This Feature adds excluded volume check to the system.
 * It is specialised to the simple cubic lattice.
 *
 * @details The excluded volume is checked by a lattice occupation algorithm
 * on a simple cubic lattice.
 * The feature is implemented as a class template, where the template parameters
 * \a LatticeClassType and \a LatticeValueTypes specify, what
 * type of lattice is used and which kind of value is saved on the lattice.
 * This feature requires FeatureLattice.
 * Notice: only implementations of Sc moves included to avoid a mix of
 * FeatureExcludedVolumeSc and MoveLocalBcc.
 *
 * @tparam <LatticeClassType<LatticeValueType>> name of the specialized class.
 * Default is FeatureLattice.
 *
 * @tparam <LatticeValueType> type of the lattice value. Default is bool.
 * */
/*****************************************************************************/

/**
 * @brief Forward declaration. Equivalent to FeatureExcludedVolumeSc< FeatureLattice <bool> >
 */
template<class LatticeClassType = FeatureLattice<> > class FeatureExcludedVolumeSc;

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

template<template<typename> class LatticeClassType, typename LatticeValueType>
class FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> > : public Feature {
public:
	//! This Feature requires a lattice.
	typedef LOKI_TYPELIST_1(LatticeClassType<LatticeValueType>) required_features_front;

	//constructor
	FeatureExcludedVolumeSc() :
			latticeFilledUp(false)
	{
	}

	/**
	 * Returns true if the underlying lattice is synchronized and all excluded volume condition
	 * (e.g. monomer/vertex occupies lattice edges) is applied.
	 * Returns false if this feature is out-of-sync.
	 *
	 * @return true if this feature is synchronized
	 * 		   false if this feature is out-of-sync.
	 **/
	bool isLatticeFilledUp() const {
		return latticeFilledUp;
	}

	/**
	 * Set's the need of synchronization of this feature e.g. escp. if the underlying lattice needs
	 * to refilled and if all excluded volume condition needs to be updated.
	 *
	 * @param[in] latticeFilledUp Specified if ExVol should be refilled (false) or everything is in-sync (true).
	 *
	 **/
	void setLatticeFilledUp(bool latticeFilledUp) {
		this->latticeFilledUp = latticeFilledUp;
	}

	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

	//! check move for sc local move
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalSc& move) const;

	//! check move for sc local diagonal move 
	template<class IngredientsType>	
	bool checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move ) const;
	//! check bcc move: Throw error if wrong lattice Type is used
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const;

	//! check move for adding an sc monomer
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const;

	//! check bcc addmove: Throw error if wrong lattice Type is used
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const;

	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move);

	//! apply move for local sc moves
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);
	
	//! apply move for local sc diagonal moves 
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalScDiag& move);
	
	//! apply move for adding an sc monomer
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);

	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

protected:

	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;

};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveBase& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveBase& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing,const MoveBase& move)
{

}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalBcc& move )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveLocalBcc, which is the wrong one in this case
 * @return false, throws an exception
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const
{
	throw std::runtime_error("*****FeatureExcludedVolumeSc::check MoveLocalBcc: wrong lattice type ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerBcc& addmove )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveAddMonomerBcc, which is the wrong one in this case
 * @return false, throws an exception
 */
/*****************************************************************************
*/

template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType, class TagType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& addmove) const
{
	throw std::runtime_error("*****FeatureExcludedVolumeSc::check MoveAddMonomerBcc: wrong lattice type ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalSc& move )const
 * @brief checks excluded volume for moves of type MoveLocalSc
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * @return if move is allowed (\a true) or rejected (\a false).
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template < class IngredientsType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalSc& move ) const
{
	if(!latticeFilledUp)
	  throw std::runtime_error("*****FeatureExcludedVolumeSc::checkMove....lattice is not populated. Run synchronize!\n");

	//get the position of the monomer to be moved (assume "lower left corner")
	VectorInt3 refPos=ingredients.getMolecules()[move.getIndex()];
	//get the direction of the move
	VectorInt3 direction=move.getDir();

	//shift refPos in the direction of the move. this defines then the
	//point around which the volume occupation must be checked
	if(direction.getX()>0 || direction.getY()>0 || direction.getZ()>0) refPos+=direction;
	refPos+=direction;

	//now the nine closest positions in the plane of refPos perpendicular to the
	//move direction must be checked

	/*get two directions perpendicular to vector directon of the move*/
	VectorInt3 perp1,perp2;
  	/* first perpendicular direction is either (0 1 0) or (1 0 0)*/
	int32_t x1=((direction.getX()==0) ? 1 : 0);
	int32_t y1=((direction.getX()!=0) ? 1 : 0);
	perp1.setX(x1);
	perp1.setY(y1);
	perp1.setZ(0);

	/* second perpendicular direction is either (0 0 1) or (0 1 0)*/
	int32_t y2=((direction.getZ()==0) ? 0 : 1);
	int32_t z2=((direction.getZ()!=0) ? 0 : 1);
	perp2.setX(0);
	perp2.setY(y2);
	perp2.setZ(z2);


    //check if the lattice sites are free
    if	(
            (ingredients.getLatticeEntry(refPos))||
            (ingredients.getLatticeEntry(refPos+perp1))||
            (ingredients.getLatticeEntry(refPos+perp2))||
            (ingredients.getLatticeEntry(refPos+perp1+perp2))
            ) return false;
    else return true;

}


/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move )const
 * @brief checks excluded volume for moves of type MoveLocalScDiag
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * @return if move is allowed (\a true) or rejected (\a false).
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template < class IngredientsType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move ) const
{
	if(!latticeFilledUp)
	    throw std::runtime_error("*****FeatureExcludedVolumeSc::checkMove....lattice is not populated. Run synchronize!\n");
	if ( move.getDir()==VectorInt3(1,0,0)||move.getDir()==VectorInt3(-1,0,0)||
	     move.getDir()==VectorInt3(0,1,0)||move.getDir()==VectorInt3(0,-1,0)||
	     move.getDir()==VectorInt3(0,0,1)||move.getDir()==VectorInt3(0,0,-1) )
	{
		  //get the position of the monomer to be moved (assume "lower left corner")
		  VectorInt3 refPos=ingredients.getMolecules()[move.getIndex()];
		  //get the direction of the move
		  VectorInt3 direction=move.getDir();

		  //shift refPos in the direction of the move. this defines then the
		  //point around which the volume occupation must be checked
		  if(direction.getX()>0 || direction.getY()>0 || direction.getZ()>0) refPos+=direction;
		  refPos+=direction;

		  //now the nine closest positions in the plane of refPos perpendicular to the
		  //move direction must be checked

		  /*get two directions perpendicular to vector directon of the move*/
		  VectorInt3 perp1,perp2;
		  /* first perpendicular direction is either (0 1 0) or (1 0 0)*/
		  int32_t x1=((direction.getX()==0) ? 1 : 0);
		  int32_t y1=((direction.getX()!=0) ? 1 : 0);
		  perp1.setX(x1);
		  perp1.setY(y1);
		  perp1.setZ(0);

		  /* second perpendicular direction is either (0 0 1) or (0 1 0)*/
		  int32_t y2=((direction.getZ()==0) ? 0 : 1);
		  int32_t z2=((direction.getZ()!=0) ? 0 : 1);
		  perp2.setX(0);
		  perp2.setY(y2);
		  perp2.setZ(z2);


	      //check if the lattice sites are free
	      if(     (ingredients.getLatticeEntry(refPos))||
		      (ingredients.getLatticeEntry(refPos+perp1))||
		      (ingredients.getLatticeEntry(refPos+perp2))||
		      (ingredients.getLatticeEntry(refPos+perp1+perp2))
		) return false;
	      else return true;
	}
	else 
	{
		//get the position of the monomer to be moved (assume "lower left corner")
		VectorInt3 refPos=ingredients.getMolecules()[move.getIndex()];
		//get the direction of the move
		VectorInt3 direction=move.getDir();
		VectorInt3 vec1, vec2,vec3;
		//vec1 and vec2 are the orthogonal projections of the direction on the
		//coordinate system. vec3 is perpendicular to both other vectors. 
		if(direction.getX()==0){
		  vec3=VectorInt3(1,0,0);
		  if(direction.getY()==1){vec1=VectorInt3(0,1,0);refPos+=vec1;}
		  else{vec1=VectorInt3(0,-1,0);}
		  if(direction.getZ()==1){vec2=VectorInt3(0,0,1);refPos+=vec2;}
		  else{vec2=VectorInt3(0,0,-1); }
		}else if(direction.getY()==0){
		  vec3=VectorInt3(0,1,0);
		  if(direction.getX()==1){vec1=VectorInt3(1,0,0);refPos+=vec1;}
		  else{vec1=VectorInt3(-1,0,0); }
		  if(direction.getZ()==1){vec2=VectorInt3(0,0,1);refPos+=vec2;}
		  else{vec2=VectorInt3(0,0,-1); }
		}else if(direction.getZ()==0){
		  vec3=VectorInt3(0,0,1);
		  if(direction.getY()==1){vec1=VectorInt3(0,1,0);refPos+=vec1;}
		  else{vec1=VectorInt3(0,-1,0); }
		  if(direction.getX()==1){vec2=VectorInt3(1,0,0);refPos+=vec2;}
		  else{vec2=VectorInt3(-1,0,0); }
		}
		//check if the lattice sites are free
		if(
		  (ingredients.getLatticeEntry(refPos+vec1))||
		  (ingredients.getLatticeEntry(refPos+vec2))||
		  (ingredients.getLatticeEntry(refPos+vec1+vec2))||
		  (ingredients.getLatticeEntry(refPos+vec1+vec2+vec3)||
		  (ingredients.getLatticeEntry(refPos+vec1+vec3))||
		  (ingredients.getLatticeEntry(refPos+vec2+vec3)))
		  ) return false;
		else return true;
	}
	
}
/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalSc<LocalMoveType>& move)
 * @brief Updates the lattice occupation according to the move for moves of type MoveLocalSc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalSc& move)
{
	//get old position and direction of the move
	VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
	VectorInt3 direction=move.getDir();

	/*get two directions perpendicular to vector directon of the move*/
	VectorInt3 perp1,perp2;
  	/* first perpendicular direction is either (0 1 0) or (1 0 0)*/
	int32_t x1=((direction.getX()==0) ? 1 : 0);
	int32_t y1=((direction.getX()!=0) ? 1 : 0);
	perp1.setX(x1);
	perp1.setY(y1);
	perp1.setZ(0);

	/* second perpendicular direction is either (0 0 1) or (0 1 0)*/
	int32_t y2=((direction.getZ()==0) ? 0 : 1);
	int32_t z2=((direction.getZ()!=0) ? 0 : 1);
	perp2.setX(0);
	perp2.setY(y2);
	perp2.setZ(z2);


	if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) oldPos-=direction;
	direction*=2;

    VectorInt3 oldPlusDir=oldPos+direction;
	//change lattice occupation accordingly
    ing.moveOnLattice(oldPos,oldPlusDir);
    ing.moveOnLattice(oldPos+perp1,oldPlusDir+perp1);
    ing.moveOnLattice(oldPos+perp2,oldPlusDir+perp2);
    ing.moveOnLattice(oldPos+perp1+perp2,oldPlusDir+perp1+perp2);

}
/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
 * @brief updates the lattice ocupation according to the move for types MoveLocalScDiag 
 * 
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeSc<LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
{       
	if(move.getDir()==VectorInt3(1,0,0)||move.getDir()==VectorInt3(-1,0,0)||
	   move.getDir()==VectorInt3(0,1,0)||move.getDir()==VectorInt3(0,-1,0)||
	   move.getDir()==VectorInt3(0,0,1)||move.getDir()==VectorInt3(0,0,-1))
	{
	    //get old position and direction of the move
	    VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
	    VectorInt3 direction=move.getDir();	
	    
	    /*get two directions perpendicular to vector directon of the move*/
	    VectorInt3 perp1,perp2;
	    /* first perpendicular direction is either (0 1 0) or (1 0 0)*/
	    int32_t x1=((direction.getX()==0) ? 1 : 0);
	    int32_t y1=((direction.getX()!=0) ? 1 : 0);
	    perp1.setX(x1);
	    perp1.setY(y1);
	    perp1.setZ(0);
	    
	    /* second perpendicular direction is either (0 0 1) or (0 1 0)*/
	    int32_t y2=((direction.getZ()==0) ? 0 : 1);
	    int32_t z2=((direction.getZ()!=0) ? 0 : 1);
	    perp2.setX(0);
	    perp2.setY(y2);
	    perp2.setZ(z2);
	    
	    
	    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) oldPos-=direction;
	    direction*=2;
	    
	    VectorInt3 oldPlusDir=oldPos+direction;
	    
	    //change lattice occupation accordingly
	    ing.moveOnLattice(oldPos,oldPlusDir);
	    ing.moveOnLattice(oldPos+perp1,oldPlusDir+perp1);
	    ing.moveOnLattice(oldPos+perp2,oldPlusDir+perp2);
	    ing.moveOnLattice(oldPos+perp1+perp2,oldPlusDir+perp1+perp2);
	}
	else
	{
	    //get the position of the monomer to be moved (assume "lower left corner")
	    VectorInt3 refPos=ing.getMolecules()[move.getIndex()];
	    //get the direction of the move
	    VectorInt3 direction=move.getDir();
	    VectorInt3 vec1, vec2,vec3;
	    //vec1 and vec2 are the orthogonal projections of the direction on the
	    //coordinate system. vec3 is perpendicular to both other vectors. 
	    if(direction.getX()==0){
	      vec3=VectorInt3(1,0,0);
	      if(direction.getY()==1){vec1=VectorInt3(0,1,0);refPos+=vec1;}
	      else{vec1=VectorInt3(0,-1,0);}
	      if(direction.getZ()==1){vec2=VectorInt3(0,0,1);refPos+=vec2;}
	      else{vec2=VectorInt3(0,0,-1); }
	    }else if(direction.getY()==0){
	      vec3=VectorInt3(0,1,0);
	      if(direction.getX()==1){vec1=VectorInt3(1,0,0);refPos+=vec1;}
	      else{vec1=VectorInt3(-1,0,0); }
	      if(direction.getZ()==1){vec2=VectorInt3(0,0,1);refPos+=vec2;}
	      else{vec2=VectorInt3(0,0,-1); }
	    }else if(direction.getZ()==0){
	      vec3=VectorInt3(0,0,1);
	      if(direction.getY()==1){vec1=VectorInt3(0,1,0);refPos+=vec1;}
	      else{vec1=VectorInt3(0,-1,0); }
	      if(direction.getX()==1){vec2=VectorInt3(1,0,0);refPos+=vec2;}
	      else{vec2=VectorInt3(-1,0,0); }
	    }
	    //change lattice occupation accordingly
	    //it behaves like a mirror transformation
	    ing.moveOnLattice(refPos-vec1,refPos+vec1);
	    ing.moveOnLattice(refPos-vec2,refPos+vec2);
	    ing.moveOnLattice(refPos-vec1-vec2,refPos+vec1+vec2);
	    ing.moveOnLattice(refPos-vec1+vec3,refPos+vec1+vec3);
	    ing.moveOnLattice(refPos-vec2+vec3,refPos+vec2+vec3);
	    ing.moveOnLattice(refPos-vec1-vec2+vec3,refPos+vec1+vec2+vec3);
	}
}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief check excluded volume for insertion of a monomer on a sc lattice
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerSc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template < class IngredientsType, class TagType>
bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move ) const
{
  if (!latticeFilledUp)
	throw std::runtime_error("*****FeatureExcludedVolumeSc::checkMove....lattice is not populated. Run synchronize!\n");

  //check if the lattice sites are free
  VectorInt3 pos=move.getPosition();
  VectorInt3 dx(1,0,0);
  VectorInt3 dy(0,1,0);
  VectorInt3 dz(0,0,1);

  if	(
	 (ingredients.getLatticeEntry(pos))||
	 (ingredients.getLatticeEntry(pos+dx))||
	 (ingredients.getLatticeEntry(pos+dy))||
	 (ingredients.getLatticeEntry(pos+dx+dy))||
	 (ingredients.getLatticeEntry(pos+dz))||
	 (ingredients.getLatticeEntry(pos+dz+dx))||
	 (ingredients.getLatticeEntry(pos+dz+dy))||
	 (ingredients.getLatticeEntry(pos+dz+dx+dy))
	 ) return false;
  else
    return true;
}


/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief apply excluded volume for MoveAddMonomerSc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerSc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType, class TagType>
void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move)
{
  VectorInt3 pos=move.getPosition();
  VectorInt3 dx(1,0,0);
  VectorInt3 dy(0,1,0);
  VectorInt3 dz(0,0,1);

  ing.setLatticeEntry(pos,1);
  ing.setLatticeEntry(pos+dx,1);
  ing.setLatticeEntry(pos+dy,1);
  ing.setLatticeEntry(pos+dx+dy,1);
  ing.setLatticeEntry(pos+dz,1);
  ing.setLatticeEntry(pos+dz+dx,1);
  ing.setLatticeEntry(pos+dz+dy,1);
  ing.setLatticeEntry(pos+dz+dx+dy,1);
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::synchronize(IngredientsType& ingredients)
{
	//note: the lattice entries are set to 0 before by the
	//synchronize function of FeatureLattice
	std::cout << "FeatureExcludedVolumeSc::synchronizing lattice occupation...\n";
	fillLattice(ingredients);
	std::cout << "done\n";
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::fillLattice(IngredientsType& ingredients)
 * @brief This function populates the lattice directly with positions from molecules.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::fillLattice(IngredientsType& ingredients)
{
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	//copy the lattice occupation from the monomer coordinates
	for(size_t n=0;n<molecules.size();n++)
	{
		VectorInt3 pos=ingredients.getMolecules()[n];

		if( ingredients.getLatticeEntry(pos)!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,0,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,1,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,0,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,1,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,0,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,1,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,1,1))!=0
		)
		{
			throw std::runtime_error("********** FeatureExcludedVolume::fillLattice: multiple lattice occupation ******************");
		}
		else
		{
			//here we simply set a one on every occupied lattice
			//site. this assumes that 1 can be cast to  LatticeClassType<LatticeValueType> ,
			//even though in principle  LatticeClassType<LatticeValueType>  could be anything.
			//note that this may be just a preliminiary initialization,
			//as other features may assign more specific values to the
			//lattice site.
			VectorInt3 pos=molecules[n];

			ingredients.setLatticeEntry(pos,1);
			ingredients.setLatticeEntry(pos+VectorInt3(1,0,0),1);
			ingredients.setLatticeEntry(pos+VectorInt3(0,1,0),1);
			ingredients.setLatticeEntry(pos+VectorInt3(1,1,0),1);
			ingredients.setLatticeEntry(pos+VectorInt3(0,0,1),1);
			ingredients.setLatticeEntry(pos+VectorInt3(1,0,1),1);
			ingredients.setLatticeEntry(pos+VectorInt3(0,1,1),1);
			ingredients.setLatticeEntry(pos+VectorInt3(1,1,1),1);
		}

	}
	latticeFilledUp=true;
}

#endif
