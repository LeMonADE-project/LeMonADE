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

#ifndef FEATURE_SYSTEM_INFORMATION_TENDOMER_H
#define FEATURE_SYSTEM_INFORMATION_TENDOMER_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>


class FeatureSystemInformationTendomer:public Feature
{
public:

	FeatureSystemInformationTendomer(){
	  nTendomers=0;
	  nCrossLinkers=0;
	  nMonomersPerChain=0;
	  nLabelsPerTendomerArm=0;
	};
	virtual ~FeatureSystemInformationTendomer(){};
	
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
	const uint32_t getNumTendomers() const {return nTendomers;} 
	//!
	void setNumTendomers(uint32_t nTendomers_){nTendomers=nTendomers_;}

	//!
	const uint32_t getNumCrossLinkers() const { return nCrossLinkers;}
	//!
	void setNumCrossLinkers(uint32_t nCrossLinkers_){nCrossLinkers=nCrossLinkers_;}
	
	//!
	const uint32_t getNumMonomersPerChain() const { return nMonomersPerChain;}
	//!
	void setNumMonomersPerChain(uint32_t nMonomersPerChain_){nMonomersPerChain=nMonomersPerChain_;}
	
	//!
	const uint32_t getNumLabelsPerTendomerArm() const { return nLabelsPerTendomerArm;}
	//!
	void setNumLabelsPerTendomerArm(uint32_t nLabelsPerTendomerArm_){nLabelsPerTendomerArm=nLabelsPerTendomerArm_;}
	
protected:
	uint32_t nTendomers;
	uint32_t nCrossLinkers;
	uint32_t nMonomersPerChain;
	uint32_t nLabelsPerTendomerArm;

};

/*****************************************************************/
/**
 * @class ReadNumOfTendomers
 *
 * @brief Handles BFM-File-Reads \b #!number_of_tendomers
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNumOfTendomers: public ReadToDestination<IngredientsType>
{
public:
  ReadNumOfTendomers(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNumOfTendomers(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNumOfTendomers<IngredientsType>::execute()
{
	std::cout<<"reading number of tendomers...";

	uint32_t Chains = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	Chains = atof(line.c_str());
	std::cout << "#!number_of_tendomers=" << (Chains) << std::endl;

	ingredients.setNumTendomers(Chains);
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
	std::cout<<"reading number of crosslinkers...";

	uint32_t CrossLinkers = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	CrossLinkers = atof(line.c_str());
	std::cout << "#!number_of_crosslinkers=" << (CrossLinkers) << std::endl;

	ingredients.setNumCrossLinkers(CrossLinkers);
}
/*****************************************************************/
/**
 * @class ReadNumLabelsPerArm
 *
 * @brief Handles BFM-File-Reads \b #!number_of_labels_per_arm
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNumLabelsPerArm: public ReadToDestination<IngredientsType>
{
public:
  ReadNumLabelsPerArm(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNumLabelsPerArm(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNumLabelsPerArm<IngredientsType>::execute()
{
	std::cout<<"reading number of sliding rings...";

	uint32_t nLabels = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	nLabels = atof(line.c_str());
	std::cout << "#!number_of_labels_per_arm=" << (nLabels) << std::endl;

	ingredients.setNumLabelsPerTendomerArm(nLabels);
}
/*****************************************************************/
/**
 * @class ReadNumOfMonomersPerChain
 *
 * @brief Handles BFM-File-Reads \b #!number_of_monomers_per_chain
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadNumOfMonomersPerChain: public ReadToDestination<IngredientsType>
{
public:
  ReadNumOfMonomersPerChain(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadNumOfMonomersPerChain(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadNumOfMonomersPerChain<IngredientsType>::execute()
{
	std::cout<<"reading chain length...";

	uint32_t chainLength = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	chainLength = atof(line.c_str());
	std::cout << "#!number_of_monomers_per_chain=" << (chainLength) << std::endl;

	ingredients.setNumMonomersPerChain(chainLength);
}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_tendomers
 * * #!number_of_crosslinkers
 * * #!number_of_labels_per_arm
 * * #!number_of_monomers_per_chain
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationTendomer::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!number_of_tendomers", new ReadNumOfTendomers<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!number_of_crosslinkers", new ReadCrossLinkerNumber<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!number_of_labels_per_arm", new ReadNumLabelsPerArm<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!number_of_monomers_per_chain", new ReadNumOfMonomersPerChain<IngredientsType>(fileReader.getDestination()));
    
  
}

/*****************************************************************/
/**
 * @class WriterNumOfTendomers
 *
 * @brief Handles BFM-File-Write \b #!number_of_tendomers
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNumOfTendomers:public AbstractWrite<IngredientsType>
{
public:
	WriterNumOfTendomers(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
	//  WriteBonds(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

	virtual ~WriterNumOfTendomers(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNumOfTendomers<IngredientsType>::writeStream(std::ostream& stream)
{
   //get reference to molecules

	stream<<"#!number_of_tendomers=" << (this->getSource().getNumTendomers()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterCrossLinkerNumber
 *
 * @brief Handles BFM-File-Write \b #!number_of_crosslinkers
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
	stream<<"#!number_of_crosslinkers=" << (this->getSource().getNumCrossLinkers()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterNumLabelsPerArm
 *
 * @brief Handles BFM-File-Write \b #!number_of_labels_per_arm
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNumLabelsPerArm:public AbstractWrite<IngredientsType>
{
public:
	WriterNumLabelsPerArm(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterNumLabelsPerArm(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNumLabelsPerArm<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!number_of_labels_per_arm=" << (this->getSource().getNumLabelsPerTendomerArm()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterNumOfMonomersPerChain
 *
 * @brief Handles BFM-File-Write \b #!number_of_monomers_per_chain
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterNumOfMonomersPerChain:public AbstractWrite<IngredientsType>
{
public:
	WriterNumOfMonomersPerChain(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriterNumOfMonomersPerChain(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterNumOfMonomersPerChain<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!number_of_monomers_per_chain=" << (this->getSource().getNumMonomersPerChain()) << std::endl<< std::endl;
}
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!number_of_tendomers
 * * #!number_of_crosslinkers
 * * #!number_of_labels_per_arm
 * * #!number_of_monomers_per_chain
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationTendomer::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& filewriter) const
{
    const IngredientsType& source=filewriter.getIngredients_();
    filewriter.registerWrite("#!number_of_tendomers", new WriterNumOfTendomers<IngredientsType>(source));
    filewriter.registerWrite("#!number_of_crosslinkers", new WriterCrossLinkerNumber<IngredientsType>(source));
    filewriter.registerWrite("#!number_of_labels_per_arm", new WriterNumLabelsPerArm<IngredientsType>(source));
    filewriter.registerWrite("#!number_of_monomers_per_chain", new WriterNumOfMonomersPerChain<IngredientsType>(source));

}

#endif /*FEATURE_SYSTEM_INFORMATION_TENDOMER_H*/
