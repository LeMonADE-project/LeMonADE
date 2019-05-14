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

#ifndef LEMONADE_CORE_MOLECULESWRITE_H
#define LEMONADE_CORE_MOLECULESWRITE_H

#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/utility/Vector3D.h>


/***********************************************************************/
/**
 * @file
 * @brief Writing routines for bfm-Reads !number_of_monomers !bonds, !mcs,
 * !add_bonds and !remove_bonds
 **/
/***********************************************************************/

/***********************************************************************/
/**
 * @class WriteNrOfMonomers
 * @brief Handles BFM-File-Write \b !number_of_monomers.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteNrOfMonomers: public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !number_of_monomers into the header of the bfm-file.
  WriteNrOfMonomers(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !number_of_monomers.
  void writeStream(std::ostream& strm){
    strm<<"!number_of_monomers="<<this->getSource().getMolecules().size()<<"\n\n";
  }
};

/***********************************************************************/
/**
 * @class WriteMcs
 * @brief Handles BFM-File-Write \b !mcs.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteMcs: public AbstractWrite<IngredientsType>
{
public:
	//! Writes \b !mcs in every output of the bfm-file.
	WriteMcs(const IngredientsType& src):AbstractWrite<IngredientsType>(src){};
	void writeStream(std::ostream& strm);

private:
	std::string writeSolventBlock(std::pair<size_t,size_t>) const;
	std::string compressNumber(int32_t) const;
	int32_t foldBack(int32_t,uint32_t) const;
	//!maximum length of a line in compressed solvent format
	static const int32_t maxLineLength=200;
};

/*******************Implementation of members  ******************************/

//! Executes the routine to write \b !mcs.
template <class IngredientsType>
void WriteMcs<IngredientsType>::writeStream(std::ostream& strm)
{
	//get references to molecules
	const typename IngredientsType::molecules_type& molecules(this->getSource().getMolecules());

	//get reference to map containing the indices of particles which
	//are written in a compressed fashion (solvent)
	const std::map<size_t,size_t>& compressedIndices=this->getSource().getCompressedOutputIndices();
	std::map<size_t,size_t>::const_iterator itCompressedIndices;
	//iterator always points to the next pair of compressed indices
	if(compressedIndices.size()!=0) itCompressedIndices=compressedIndices.begin();
	else itCompressedIndices=compressedIndices.end();

	//first the mcs is written into this stringstream
	std::stringstream contents;

	//write command and age
	contents<<"\n!mcs="<< molecules.getAge();

	//write monomers, bonds,solvents
	for(size_t n=0;n< molecules.size();n++)
	{
		//if next monomer is solvent, write solvent block
		if(itCompressedIndices!=compressedIndices.end())
		{
			if(itCompressedIndices->first==n)
			{
				contents<<writeSolventBlock(*itCompressedIndices);
				n=itCompressedIndices->second+1;
				itCompressedIndices++;
				if(n>=molecules.size()) break;
			}
		}

		//write position of next non-solvent
		if(molecules.getNumLinks(n)==0)
		{
			contents<<"\n";
			contents << molecules[n].getVector3D();
		}
		else if(n==0)
		{
			contents<<"\n";
			contents << molecules[n].getVector3D();
			contents<<" ";
		}
		//if the actual monomer is a chainstart, start a new subchain in !mcs Read
		else if (!molecules.areConnected(n,n-1)){
			contents<<"\n";
			contents << molecules[n].getVector3D();
			contents << " ";
		}
		//otherwise write the bond
		else
		{
			contents<<char(this->getSource().getBondset().getBondIdentifier(
				(molecules[n].getX())-(molecules[n-1].getX()),
				(molecules[n].getY())-(molecules[n-1].getY()),
				(molecules[n].getZ())-(molecules[n-1].getZ()) ) );

			contents.flush();
		}
	}

	contents<<"\n\n";
	contents.flush();
	strm << contents.str();
	strm.flush();
}

/**
 * @brief Returns a string containing compressed solvent coordinates
 * @details Compresses the monomers in the index range given by \a solventIndices
 * into the solvent format and returns the resulting string. The compression works
 * as follows: The coordinates in the box are indexed in a linear fashion as
 * index=foldedZ*boxX*boxY+foldedY*boxX+foldedX
 * where folded(XYZ) are the coordinates of the solvent folded into the box. Then
 * what is written into the solvent string is always the distance between adjacent
 * indices. To achieve further compression, the distance is not written in
 * decimal numbers, but as the corresponding ASCII (example below). Since control
 * characters are not used, only characters ASCII=33-128 are used. This means that
 * ASCII=33 (!) corresponds to distance 0, ASCII=34(") to distance 1, and so on.
 * If a distance between two indices is larger than 94 (128-34), the distance
 * is instead written as a decimal number, where a space (ASCII=32) is used as
 * a separator before and after this decimal number.
 * Additionally, the length of one line is restricted to 200. If the line becomes
 * longer, the string continues in the next line, which then begins with "sc ",
 * standing for "solvent continue".
 * Example:
 * Monomers 0-9 are solvents at positions (10 0 0),(11 0 0)...(19 0 0)
 * Then, a range of non solvent comes, and Monomers (200 - 1000) are again solvent
 * The output would look like this
 *
 * !mcs 0
 * solvent +"""""""""
 * 2 45 67 sgsoi540jsogjw59trjgjs4
 * 40 56 77 adifja23954w4865wjgf84672sigkjxjhg
 * (...some more chains...)
 * solvent  !3JH$)"/$?"/$="jhws (until 200 characters...)
 * sc ?38()/%3w2nc 1000 gewoi8 (until 200 characters...)
 * sc (...and so on)
 *
 * In the first line starting with "solvent" there is a "+", corresponding to
 * ASCII=43, i.e. subtracting the offset of 33 leads to the index 10, which
 * corresponds to the position of monomer 0. Monomer 1 has index 11, so the
 * distance is 11-10=1, and the character corresponding to 1+33 (offset) is ".
 *
 * The next lines contain polymer chains in the regular format
 * Afterwards, there are again compressed solvent monomers, spread over several
 * lines. In the second of these lines there is a decimal number 1000 with a
 * space before and after, to indicate this is to be read as a decimal distance.
 *
 * @throw <std::runtime_error> if indices are out of range
 **/
template<class IngredientsType>
std::string WriteMcs<IngredientsType>::writeSolventBlock(std::pair<size_t,size_t> solventIndices) const
{
	const typename IngredientsType::molecules_type& molecules(this->getSource().getMolecules());
	//exception if indices out of range
	if(solventIndices.first>=molecules.size() || solventIndices.second>=molecules.size())
	{
		std::stringstream errormessage;
		errormessage<<"WriteMcs::writeSolventBlock: solvent particle indices "
			<<solventIndices.first<<" "<<solventIndices.second<<" out of range. Size of system is "<<molecules.size()<<"\n";
		throw std::runtime_error(errormessage.str());
	}
	//get box dimensions
	uint32_t boxX=this->getSource().getBoxX();
	uint32_t boxY=this->getSource().getBoxY();
	uint32_t boxZ=this->getSource().getBoxZ();

	//everything is written in a stringstream for convenience first
	std::stringstream solventBlock;
	//now translate all solvent coordinates to integer indices and save them in the map
	//the map stores how many solvents sit on which position (if there is no
	//excluded volume interaction, there could be more than one)
	std::map<int32_t,int32_t> coordinates;
	for(size_t i=solventIndices.first;i<solventIndices.second+1;i++)
	{
		//transform the position of the solvent to a linear index
		int32_t foldedX=foldBack(molecules[i].getX(),boxX);
		int32_t foldedY=foldBack(molecules[i].getY(),boxY);
		int32_t foldedZ=foldBack(molecules[i].getZ(),boxZ);
		int32_t linIndex=foldedZ*boxX*boxY+foldedY*boxX+foldedX;
		coordinates[linIndex]+=1;
	}

	//write the solvent keyword
	solventBlock<<"\nsolvent ";

	 //now write the map
	int32_t oldIndex=0;
	int32_t lineCount=0;
	typedef std::map<int32_t,int32_t>::const_iterator iterator;
	for(iterator i=coordinates.begin();i!=coordinates.end();i++)
	{
		//write the index difference in compressed format
		solventBlock<<compressNumber((i->first)-oldIndex);
		lineCount++;
		//insert new line if line has become too long
		if(lineCount>maxLineLength)
		{
			solventBlock<<"\n"<<"sc ";
			lineCount=0;
		}
		oldIndex=i->first;
		//now write a zero if there are more solvents in the same position
		for(int32_t k=1;k<i->second;k++)
		{
			solventBlock<<compressNumber(0);
			lineCount++;
			//insert new line if line has become too long
			if(lineCount>maxLineLength)
			{
				solventBlock<<"\n"<<"sc ";
				lineCount=0;
			}
		}

	}
	return solventBlock.str();
}

/***********************************************************************/
//!compresses an integer number to char(number+33) for number<95, otherwise the decimal representation is kept
template<class IngredientsType>
std::string WriteMcs<IngredientsType>::compressNumber(int32_t dist) const
{
	std::stringstream compressedNumber;
	int32_t maxSingleDistance=94; //128-33-1
	if(dist<=maxSingleDistance)
		compressedNumber<<char(dist+33); //offset of 33 for zero distance
	else
		compressedNumber<<char(32)<<dist<<char(32); //ASCII 32=space

	return compressedNumber.str();
}

/***********************************************************************/
//!folds a coordinate back into a box
template<class IngredientsType>
int32_t WriteMcs<IngredientsType>::foldBack(int32_t value, uint32_t box) const {

  //make coordinate positive, so that modulo is properly defined
  while (value < 0) {
    value += box;
  }
  //fold back the value
  return (value % box);
}



/***********************************************************************/
/**
 * @class WriteBonds
 * @brief Handles BFM-File-Write \b !bonds.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBonds: public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !bonds into the header of the bfm-file.
  WriteBonds(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  void writeStream(std::ostream& strm);

};

/**********************implementation of members    *******************/
//! Executes the routine to write \b !bonds.
template <class IngredientsType>
void WriteBonds<IngredientsType>::writeStream(std::ostream& strm){

 strm<<"!bonds\n";

 //loop over all monomers once, then check their partners that have a smaller index
 //and write an additional bond, if difference in index is larger than one (because
 //then the bond is not covered by the normal chain order), an if the monomers are
 //connected across a chainstart

 //get reference to molecules
 typename IngredientsType::molecules_type molecules=this->getSource().getMolecules();

 //loop over all monomers
 for(size_t a=1;a< molecules.size();a++){
   //act only if there are bonds
   if(molecules.getNumLinks(a)>0){
     size_t nLinks=molecules.getNumLinks(a);
     //loop over bond partners
     for(size_t b=0;b<nLinks;b++){
       //look only backwards
       if(molecules.getNeighborIdx(a,b)<a){

	 //write bond if indices are more than one apart
	 if(a-molecules.getNeighborIdx(a,b) > 1){
	  strm<<molecules.getNeighborIdx(a,b)+1<<" "<<a+1<<"\n";
#ifdef DEBUG
	  std::cout<<"found additional bond: "<<a<<" "<<molecules.getNeighborIdx(a,b)<<std::endl;
#endif //DEBUG
	 }

       }
     }
   }
 }
 strm<<"\n";
}

/******** write command handling !add_bonds *********************************/
/**
 * @class WriteAddBonds
 * @brief Handles BFM-File-Write \b !add_bonds.
 * 
 * @tparam IngredientsType Ingredients class storing all system information.
 * */
template <class IngredientsType>
class WriteAddBonds: public AbstractWrite<IngredientsType>
{
  	enum BFM_WRITE_TYPE_EXPANDED{
	  C_APPEND=3,	//!< The configuration (excl. header) is append to the file  
	  C_OVERWRITE=4,  //!< The configuration (incl. header) overwrites the existing file
	  C_NEWFILE=5,	//!< The configuration (incl. header) is written to a new file
	  C_APPNOFILE=6,	//!< The file doenst exist
	};
	
public:
	//! constructor
	WriteAddBonds(const IngredientsType& ingredients, int writeType=C_APPEND):
	AbstractWrite<IngredientsType>(ingredients),
	myWriteType(writeType),
	old_ingredients(ingredients)
	{this->setHeaderOnly(false);}
	
	//! writes to the file stream
	void writeStream(std::ostream& strm);
	
	
private:
  	  
	typedef typename IngredientsType::molecules_type::edge_type edge_type;
	//! Storage for added bonds
// 	std::map<std::pair<uint32_t,uint32_t>,edge_type> AddBonds;
	
	//!Storage for a copy of ingredients of the last time step
	IngredientsType old_ingredients;
	
	//! ENUM-type BFM_WRITE_TYPE specify the write-out
	int myWriteType;
	
};
/**********************implementation of members    *******************/
//! Executes the routine to write \b !add_bonds.
template <class IngredientsType>
void WriteAddBonds<IngredientsType>::writeStream(std::ostream& strm){
	switch(myWriteType)
	{ case C_APPNOFILE:
	  case C_NEWFILE: 
	  case C_APPEND: {

		//get a map containing the added bond	
		std::map<std::pair<uint32_t,uint32_t>,edge_type> AddBonds=this->getSource().getMolecules().getEdges();
		//get a map containing the removed bonds
		std::map<std::pair<uint32_t,uint32_t>,edge_type> RemovedBonds=old_ingredients.getMolecules().getEdges();
		
		//erases all bond parnters from the map which are unchanged
		//the rest of RemovedBonds contains only bonds which are removed during
		//the last simulation step, in contrast the rest of the map AddBonds 
		//contains only bonds which are newly formed during the last simulation
		//step
		typename std::map<std::pair<uint32_t,uint32_t>, edge_type>::iterator it;
		for(it=RemovedBonds.begin();it!=RemovedBonds.end();++it){
			if(AddBonds.find(it->first)!=AddBonds.end()){
			  AddBonds.erase(AddBonds.find(it->first));
			}
		}
		
		//write only the bonds that were added since the last update	
		strm<<"!add_bonds\n";
		for(it=AddBonds.begin();it!=AddBonds.end();++it){
			  strm<<it->first.first+1<<" "<<it->first.second+1<<"\n";
		}
		strm<<"\n";
		
		old_ingredients=this->getSource();
		break;}
	  case C_OVERWRITE: {break;} 
	}
	
}
/******** write command handling !remove_bonds *********************************/
/**
 * @class WriteRemoveBonds
 * @brief Handles BFM-File-Write \b !remove_bonds.
 * 
 * @tparam IngredientsType Ingredients class storing all system information.
 * */
template <class IngredientsType>
class WriteRemoveBonds: public AbstractWrite<IngredientsType>
{
  enum BFM_WRITE_TYPE{
	  C_APPEND=3,	//!< The configuration (excl. header) is append to the file  
	  C_OVERWRITE=4,  //!< The configuration (incl. header) overwrites the existing file
	  C_NEWFILE=5,	//!< The configuration (incl. header) is written to a new file
	  C_APPNOFILE=6,	//!< The file doenst exist
	};
	
public:
	//! constructor
	WriteRemoveBonds(const IngredientsType& ingredients, int writeType=C_APPEND):
	AbstractWrite<IngredientsType>(ingredients),
	myWriteType(writeType),
	old_ingredients(ingredients)
	{this->setHeaderOnly(false);}
	

	//! writes to the file stream
	void writeStream(std::ostream& strm);
	
private:
	  
	typedef typename IngredientsType::molecules_type::edge_type edge_type;
	  
	//! Storage for removed bonds
// 	std::map<std::pair<uint32_t,uint32_t>,edge_type> RemovedBonds;
	
	//!Storage for a copy of ingredients of the last time step
	IngredientsType old_ingredients;
	
	//! ENUM-type BFM_WRITE_TYPE specify the write-out
	int myWriteType;

};
/**********************implementation of members    *******************/
//! Executes the routine to write \b !remove_bonds.
template <class IngredientsType>
void WriteRemoveBonds<IngredientsType>::writeStream(std::ostream& strm){
	switch(myWriteType)
	{ case C_APPNOFILE:
	  case C_NEWFILE: 
	  case C_APPEND: {
	  
	      //get a map containing the removed bonds
	      std::map<std::pair<uint32_t,uint32_t>,edge_type> RemovedBonds=old_ingredients.getMolecules().getEdges();
	      //get a map containing the added bond	
	      std::map<std::pair<uint32_t,uint32_t>,edge_type> AddBonds=this->getSource().getMolecules().getEdges();
	      
	      //erases all bond parnters from the map which are unchanged
	      //the rest of RemovedBonds contains only bonds which are removed during
	      //the last simulation step, in contrast the rest of the map AddBonds 
	      //contains only bonds which are newly formed during the last simulation
	      //step
	      typename std::map<std::pair<uint32_t,uint32_t>, edge_type>::iterator it;
	      for(it=AddBonds.begin();it!=AddBonds.end();++it){
		      if(RemovedBonds.find(it->first)!=RemovedBonds.end()){
			RemovedBonds.erase(RemovedBonds.find(it->first));
		      }
	      }
	      
	      //write only the breaks that were removed since the last update
	      strm<<"!remove_bonds\n";
	      for(it=RemovedBonds.begin();it!=RemovedBonds.end();++it){
		      strm<<it->first.first+1<<" "<<it->first.second+1<<"\n";
	      }
	      strm<<"\n";
	      
	      old_ingredients=this->getSource();
	      break;}
	  case C_OVERWRITE: {break;} 
	}


	
}
#endif /* LEMONADE_CORE_MOLECULESWRITE_H */
