/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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

#ifndef FEATURE_SYSTEM_INFORMATION_RING_MELT_H
#define FEATURE_SYSTEM_INFORMATION_RING_MELT_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
/**----------------------------------------------------------------------------
 * @file FeatureSystemInformationRingMelt
 * @brief Provides the read/write functions for the number of rings and 
 * nMonomersPerRing
 * @details The function getNumRings() returns the number of rings and the 
 * function getNumMonomersPerRing() the number the monomers per ring. 
 * There are setter functions for setting these variables. This feature is 
 * purely for convinience and metadata handling.
 * @author Toni Müller
 *---------------------------------------------------------------------------*/

class FeatureSystemInformationRingMelt:public Feature
{
public:

	FeatureSystemInformationRingMelt():nRings(0), nMonomersPerRing(0) {};
	
	virtual ~FeatureSystemInformationRingMelt(){};
	
	//! For all unknown moves: this does nothing
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const {return true; }
	
	//! Overloaded for moves of type MoveScMonomer to check for sinusoidal movement
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const {return true; }

	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);
	
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& filewriter) const;

	//!
	const uint32_t getNumRings() const {return nRings;} 
	//!
	void setNumRings(uint32_t nRings_){nRings=nRings_;}

	//!
	const uint32_t getNumMonomersPerRing() const { return nMonomersPerRing;}
	//!
	void setNumMonomersPerRing(uint32_t nMonomersPerRing_){nMonomersPerRing=nMonomersPerRing_;}
		
protected:
	uint32_t nRings;
	uint32_t nMonomersPerRing;

};

/*****************************************************************/
/**
 * @class ReadNumOfRings
 *
 * @brief Handles BFM-File-Reads \b #!number_of_rings
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNumOfRings: public ReadToDestination<IngredientsType>
{
public:
  ReadNumOfRings(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNumOfRings(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNumOfRings<IngredientsType>::execute()
{
	std::cout<<"reading number of ring...";

	uint32_t Chains = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	Chains = atof(line.c_str());
	std::cout << "#!number_of_rings=" << (Chains) << std::endl;

	ingredients.setNumRings(Chains);
}

/*****************************************************************/
/**
 * @class ReadNumOfMonomersPerRing
 *
 * @brief Handles BFM-File-Reads \b #!number_of_monomers_per_ring#
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNumOfMonomersPerRing: public ReadToDestination<IngredientsType>
{
public:
  ReadNumOfMonomersPerRing(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNumOfMonomersPerRing(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNumOfMonomersPerRing<IngredientsType>::execute()
{
	std::cout<<"reading ring length...";

	uint32_t ringLength = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	ringLength = atof(line.c_str());
	std::cout << "#!number_of_monomers_per_ring=" << (ringLength) << std::endl;

	ingredients.setNumMonomersPerRing(ringLength);
}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_rings
 * * #!number_of_monomers_per_ring
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationRingMelt::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!number_of_rings", new ReadNumOfRings<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!number_of_monomers_per_ring", new ReadNumOfMonomersPerRing<IngredientsType>(fileReader.getDestination()));
    
  
}

/*****************************************************************/
/**
 * @class WriterNumOfRings
 *
 * @brief Handles BFM-File-Write \b #!number_of_rings
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNumOfRings:public AbstractWrite<IngredientsType>
{
public:
	WriterNumOfRings(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
	//  WriteBonds(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

	virtual ~WriterNumOfRings(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNumOfRings<IngredientsType>::writeStream(std::ostream& stream)
{
   //get reference to molecules

	stream<<"#!number_of_rings=" << (this->getSource().getNumRings()) << std::endl<< std::endl;
}

/*****************************************************************/
/**
 * @class WriterNumOfMonomersPerRing
 *
 * @brief Handles BFM-File-Write \b #!number_of_monomers_per_ring
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNumOfMonomersPerRing:public AbstractWrite<IngredientsType>
{
public:
	WriterNumOfMonomersPerRing(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterNumOfMonomersPerRing(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNumOfMonomersPerRing<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!number_of_monomers_per_ring=" << (this->getSource().getNumMonomersPerRing()) << std::endl<< std::endl;
}
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_rings
 * * #!number_of_monomers_per_ring
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationRingMelt::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& filewriter) const
{
    const IngredientsType& source=filewriter.getIngredients_();
    filewriter.registerWrite("#!number_of_rings", new WriterNumOfRings<IngredientsType>(source));
    filewriter.registerWrite("#!number_of_monomers_per_ring", new WriterNumOfMonomersPerRing<IngredientsType>(source));

}

#endif /*FEATURE_SYSTEM_INFORMATION_RING_MELT_H*/
