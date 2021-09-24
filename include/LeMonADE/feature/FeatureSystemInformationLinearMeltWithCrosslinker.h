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

#ifndef FEATURE_SYSTEM_INFORMATION_LINEAR_MELT_WITH_CROSSLINKERH
#define FEATURE_SYSTEM_INFORMATION_LINEAR_MELT_WITH_CROSSLINKERH

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>


class FeatureSystemInformationLinearMeltWithCrosslinker:public Feature
{
public:

	FeatureSystemInformationLinearMeltWithCrosslinker(){
	  nChains=0;
	  nCrossLinkers=0;
	  chainLength=0;
	  functionality=0;
	  nMonomersPerCrossLink=0;
	};
	virtual ~FeatureSystemInformationLinearMeltWithCrosslinker(){};
	
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

	//!return the number of chains 
	uint32_t getNumOfChains() const { return nChains; } 
	
	//!set the number of chains 
	void setNumOfChains(uint32_t nChains_){nChains=nChains_;}

	//! return the number of cross linkers 
	uint32_t getNumOfCrosslinks() const { return nCrossLinkers; }
	
	//!set the number of cross linekrs 
	void setNumOfCrosslinks(uint32_t nCrossLinkers_){nCrossLinkers=nCrossLinkers_;}
	
	//! returns the chain lenght 
	const uint32_t getNumOfMonomersPerChain() const {return chainLength; }
	
	//!set the chain length 
	void setNumOfMonomersPerChain(uint32_t chainLength_){chainLength=chainLength_;}
	
	//! get the functionality of the cross linker 
	uint32_t getFunctionality() const {return functionality;}
	
	//!set the functionality of cross linker
	void setFunctionality(uint32_t functionality_){functionality=functionality_;}
	
	//!get the number of monomres per cross link
	uint32_t getNumOfMonomersPerCrosslink() const { return nMonomersPerCrossLink;}
	
	//!set the number of monomers per cross linker 
	void setNumOfMonomersPerCrosslink(uint32_t nMonomersPerCrossLink_){nMonomersPerCrossLink=nMonomersPerCrossLink_;}
	
private:
	uint32_t nChains;
	uint32_t nCrossLinkers;
	uint32_t chainLength;
	uint32_t functionality;
	uint32_t nMonomersPerCrossLink;
	
};

/*****************************************************************/
/**
 * @class ReadLinearChains
 *
 * @brief Handles BFM-File-Reads \b #!number_of_linear_chains
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadLinearChains: public ReadToDestination<IngredientsType>
{
public:
  ReadLinearChains(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadLinearChains(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadLinearChains<IngredientsType>::execute()
{
	std::cout<<"reading number of Chains...";

	uint32_t Chains = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	Chains = atof(line.c_str());
	std::cout << "#!number_of_linear_chains=" << (Chains) << std::endl;

	ingredients.setNumOfChains(Chains);
}
/*****************************************************************/
/**
 * @class ReadCrossLinkerNumber
 *
 * @brief Handles BFM-File-Reads \b #!number_of_crosslinkers
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadCrossLinkerNumber: public ReadToDestination<IngredientsType>
{
public:
  ReadCrossLinkerNumber(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadCrossLinkerNumber(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadCrossLinkerNumber<IngredientsType>::execute()
{
	std::cout<<"reading number of Crosslinkers...";

	uint32_t CrossLinkers = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	CrossLinkers = atof(line.c_str());
	std::cout << "#!number_of_crosslinkers=" << (CrossLinkers) << std::endl;

	ingredients.setNumOfCrosslinks(CrossLinkers);
}
/*****************************************************************/
/**
 * @class ReadChainLength
 *
 * @brief Handles BFM-File-Reads \b #!chainLength
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadChainLength: public ReadToDestination<IngredientsType>
{
public:
  ReadChainLength(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadChainLength(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadChainLength<IngredientsType>::execute()
{
	std::cout<<"reading chain length...";

	uint32_t chainLength = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	chainLength = atof(line.c_str());
	std::cout << "#!chainLength=" << (chainLength) << std::endl;

	ingredients.setNumOfMonomersPerChain(chainLength);
}

/*****************************************************************/
/**
 * @class ReadFunctionality
 *
 * @brief Handles BFM-File-Reads \b #!functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadFunctionality: public ReadToDestination<IngredientsType>
{
public:
  ReadFunctionality(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadFunctionality(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadFunctionality<IngredientsType>::execute()
{
	std::cout<<"reading functionality...";

	uint32_t func = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	func = atof(line.c_str());
	std::cout << "#!functionality=" << (func) << std::endl;

	ingredients.setFunctionality(func);
}

/*****************************************************************/
/**
 * @class ReadNMonomersPerCrossLink
 *
 * @brief Handles BFM-File-Reads \b #!nMonomersPerCrossLink
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNMonomersPerCrossLink: public ReadToDestination<IngredientsType>
{
public:
  ReadNMonomersPerCrossLink(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNMonomersPerCrossLink(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNMonomersPerCrossLink<IngredientsType>::execute()
{
	std::cout<<"reading number of monomers per cross link...";

	uint32_t crosslinkweight = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	crosslinkweight = atof(line.c_str());
	std::cout << "#!nMonomersPerCrossLink=" << (crosslinkweight) << std::endl;

	ingredients.setNumOfMonomersPerCrosslink(crosslinkweight);
}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_linear_chains
 * * #!number_of_crosslinkers
 * * #!chainLength
 * * #!functionality
 * * #!nMonomersPerCrossLink
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationLinearMeltWithCrosslinker::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!number_of_linear_chains", new ReadLinearChains<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!number_of_crosslinkers", new ReadCrossLinkerNumber<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!chainLength", new ReadChainLength<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!functionality", new ReadFunctionality<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!nMonomersPerCrossLink", new ReadNMonomersPerCrossLink<IngredientsType>(fileReader.getDestination()));
}

/*****************************************************************/
/**
 * @class WriterLinearChains
 *
 * @brief Handles BFM-File-Write \b #!number_of_linear_chains
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterLinearChains:public AbstractWrite<IngredientsType>
{
public:
	WriterLinearChains(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
	//  WriteBonds(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

	virtual ~WriterLinearChains(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterLinearChains<IngredientsType>::writeStream(std::ostream& stream)
{
   //get reference to molecules

	stream<<"#!number_of_linear_chains=" << (this->getSource().getNumOfChains()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterCrossLinkerNumber
 *
 * @brief Handles BFM-File-Write \b #!number_of_linear_chains
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterCrossLinkerNumber:public AbstractWrite<IngredientsType>
{
public:
	WriterCrossLinkerNumber(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterCrossLinkerNumber(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterCrossLinkerNumber<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!number_of_crosslinkers=" << (this->getSource().getNumOfCrosslinks()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterChainLength
 *
 * @brief Handles BFM-File-Write \b #!chainLength
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterChainLength:public AbstractWrite<IngredientsType>
{
public:
	WriterChainLength(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterChainLength(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterChainLength<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!chainLength=" << (this->getSource().getNumOfMonomersPerChain()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterFunctionality
 *
 * @brief Handles BFM-File-Write \b #!functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterFunctionality:public AbstractWrite<IngredientsType>
{
public:
	WriterFunctionality(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterFunctionality(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterFunctionality<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!functionality=" << (this->getSource().getFunctionality()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterNMonomersPerCrossLink
 *
 * @brief Handles BFM-File-Write \b #!nMonomersPerCrossLink
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNMonomersPerCrossLink:public AbstractWrite<IngredientsType>
{
public:
	WriterNMonomersPerCrossLink(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterNMonomersPerCrossLink(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNMonomersPerCrossLink<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!nMonomersPerCrossLink=" << (this->getSource().getNumOfMonomersPerCrosslink()) << std::endl<< std::endl;
}
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_linear_chains
 * * #!number_of_crosslinkers
 * * #!chainLength
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationLinearMeltWithCrosslinker::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& filewriter) const
{
    const IngredientsType& source=filewriter.getIngredients_();
    filewriter.registerWrite("#!number_of_linear_chains", new WriterLinearChains<IngredientsType>(source));
    filewriter.registerWrite("#!number_of_crosslinkers", new WriterCrossLinkerNumber<IngredientsType>(source));
    filewriter.registerWrite("#!chainLength", new WriterChainLength<IngredientsType>(source));
    filewriter.registerWrite("#!functionality", new WriterFunctionality<IngredientsType>(source));
    filewriter.registerWrite("#!nMonomersPerCrossLink", new WriterNMonomersPerCrossLink<IngredientsType>(source));
}

#endif /*FEATURE_SYSTEM_INFORMATION_LINEAR_MELT_WITH_CROSSLINKERH*/
