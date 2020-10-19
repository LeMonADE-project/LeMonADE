/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2020 by
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


class FeatureSystemInformationDendrimer:public Feature
{
public:

	FeatureSystemInformationDendrimer():	generation(0), spacerLength(0), coreFunctionality(0), branchingPointFunctionality(0){};

	virtual ~FeatureSystemInformationDendrimer(){};
	
	//! For all unknown moves: this does nothing
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredientsngredients,const MoveBase& move) const {return true; }

	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);
	
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& filewriter) const;

	//! getter function for dendrimer generation
	const uint32_t getGeneration() const {return generation;} 
	//! setter function for dendrimer generation
	void setGeneration(uint32_t generation_){generation=generation_;}

	//! getter function for dendrimer spacer length
	const uint32_t getSpacerLength() const { return spacerLength;}
	//! setter function for dendrimer spacer length
	void setSpacerLength(uint32_t spacerLength_){spacerLength=spacerLength_;}
	
	//! getter function for the functionality of branching of the dendrimer focal point
	const uint32_t getCoreFunctionality() const { return coreFunctionality;}
	//! setter function for the functionality of branching of the dendrimer focal point
	void setCoreFunctionality(uint32_t coreFunctionality_){coreFunctionality=coreFunctionality_;}
	
	//! getter function for the branching functionality of the dendrimer branching points except the focal point
	const uint32_t getBranchingPointFunctionality() const { return branchingPointFunctionality;}
	//! setter function for the branching functionality of the dendrimer branching points except the focal point
	void setBranchingPointFunctionality(uint32_t branchingPointFunctionality_){branchingPointFunctionality=branchingPointFunctionality_;}
	
protected:
	uint32_t generation;
	uint32_t spacerLength;
	uint32_t coreFunctionality;
	uint32_t branchingPointFunctionality;

};

/*****************************************************************/
/**
 * @class ReadDendrimerGeneration
 *
 * @brief Handles BFM-File-Reads \b #!dendrimer_generation
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadDendrimerGeneration: public ReadToDestination<IngredientsType>
{
public:
  ReadDendrimerGeneration(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadDendrimerGeneration(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadDendrimerGeneration<IngredientsType>::execute()
{
	std::cout<<"reading dendrimer generation...";

	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	uint32_t dGeneration(std::stoi(line));
	std::cout << "#!dendrimer_generation=" << (dGeneration) << std::endl;

	ingredients.setGeneration(dGeneration);

}
/*****************************************************************/
/**
 * @class ReadDendrimerSpacerLength
 *
 * @brief Handles BFM-File-Reads \b #!dendrimer_spacer_length
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadDendrimerSpacerLength: public ReadToDestination<IngredientsType>
{
public:
  ReadDendrimerSpacerLength(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadDendrimerSpacerLength(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadDendrimerSpacerLength<IngredientsType>::execute()
{
	std::cout<<"reading dendrimer spacer length...";

	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	uint32_t dSpacerLength(std::stoi(line));
	std::cout << "#!dendrimer_spacer_length=" << (dSpacerLength) << std::endl;

	ingredients.setSpacerLength(dSpacerLength);
}
/*****************************************************************/
/**
 * @class ReadDendrimerCoreFunctionality
 *
 * @brief Handles BFM-File-Reads \b #!dendrimer_core_functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadDendrimerCoreFunctionality: public ReadToDestination<IngredientsType>
{
public:
  ReadDendrimerCoreFunctionality(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadDendrimerCoreFunctionality(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadDendrimerCoreFunctionality<IngredientsType>::execute()
{
	std::cout<<"reading dendrimer core functionality...";

	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	uint32_t dCoreFunct(std::stoi(line));
	std::cout << "#!dendrimer_core_functionality=" << (dCoreFunct) << std::endl;

	ingredients.setCoreFunctionality(dCoreFunct);
}
/*****************************************************************/
/**
 * @class ReadDendrimerBranchingPointFunctionality
 *
 * @brief Handles BFM-File-Reads \b #!dendrimer_branching_point_functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadDendrimerBranchingPointFunctionality: public ReadToDestination<IngredientsType>
{
public:
  ReadDendrimerBranchingPointFunctionality(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadDendrimerBranchingPointFunctionality(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadDendrimerBranchingPointFunctionality<IngredientsType>::execute()
{
	std::cout<<"reading dendrimer branching point functionality...";

	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	uint32_t dBraPoFunct(std::stoi(line));
	std::cout << "#!dendrimer_branching_point_functionality=" << (dBraPoFunct) << std::endl;

	ingredients.setBranchingPointFunctionality(dBraPoFunct);
}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!dendrimer_generation
 * * #!dendrimer_spacer_length
 * * #!dendrimer_core_functionality
 * * #!dendrimer_branching_point_functionality
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationDendrimer::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!dendrimer_generation", new ReadDendrimerGeneration<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!dendrimer_spacer_length", new ReadDendrimerSpacerLength<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!dendrimer_core_functionality", new ReadDendrimerCoreFunctionality<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!dendrimer_branching_point_functionality", new ReadDendrimerBranchingPointFunctionality<IngredientsType>(fileReader.getDestination()));
}

/*****************************************************************/
/**
 * @class WriterDendrimerGeneration
 *
 * @brief Handles BFM-File-Write \b #!dendrimer_generation
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterDendrimerGeneration:public AbstractWrite<IngredientsType>
{
public:
	WriterDendrimerGeneration(const IngredientsType& ingredients)
	:AbstractWrite<IngredientsType>(ingredients){this->setHeaderOnly(true);}

	virtual ~WriterDendrimerGeneration(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterDendrimerGeneration<IngredientsType>::writeStream(std::ostream& stream)
{
   //get reference to molecules
	stream<<"#!dendrimer_generation=" << (this->getSource().getGeneration()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterDendrimerSpacerLength
 *
 * @brief Handles BFM-File-Write \b #!dendrimer_spacer_length
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterDendrimerSpacerLength:public AbstractWrite<IngredientsType>
{
public:
	WriterDendrimerSpacerLength(const IngredientsType& ingredients)
	:AbstractWrite<IngredientsType>(ingredients){this->setHeaderOnly(true);}

	virtual ~WriterDendrimerSpacerLength(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterDendrimerSpacerLength<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!dendrimer_spacer_length=" << (this->getSource().getSpacerLength()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterDendrimerCoreFunctionality
 *
 * @brief Handles BFM-File-Write \b #!dendrimer_core_functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterDendrimerCoreFunctionality:public AbstractWrite<IngredientsType>
{
public:
	WriterDendrimerCoreFunctionality(const IngredientsType& ingredients)
	:AbstractWrite<IngredientsType>(ingredients){this->setHeaderOnly(true);}

	virtual ~WriterDendrimerCoreFunctionality(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterDendrimerCoreFunctionality<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!dendrimer_core_functionality=" << (this->getSource().getCoreFunctionality()) << std::endl<< std::endl;
}
/*****************************************************************/
/**
 * @class WriterDendrimerBranchingPointFunctionality
 *
 * @brief Handles BFM-File-Write \b #!dendrimer_branching_point_functionality
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriterDendrimerBranchingPointFunctionality:public AbstractWrite<IngredientsType>
{
public:
	WriterDendrimerBranchingPointFunctionality(const IngredientsType& ingredients)
	:AbstractWrite<IngredientsType>(ingredients){this->setHeaderOnly(true);}

	virtual ~WriterDendrimerBranchingPointFunctionality(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriterDendrimerBranchingPointFunctionality<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!dendrimer_branching_point_functionality=" << (this->getSource().getBranchingPointFunctionality()) << std::endl<< std::endl;
}
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!dendrimer_generation
 * * #!dendrimer_spacer_length
 * * #!dendrimer_core_functionality
 * * #!dendrimer_branching_point_functionality
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureSystemInformationDendrimer::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& filewriter) const
{
    const IngredientsType& source=filewriter.getIngredients_();
    filewriter.registerWrite("#!dendrimer_generation", new WriterDendrimerGeneration<IngredientsType>(source));
    filewriter.registerWrite("#!dendrimer_spacer_length", new WriterDendrimerSpacerLength<IngredientsType>(source));
    filewriter.registerWrite("#!dendrimer_core_functionality", new WriterDendrimerCoreFunctionality<IngredientsType>(source));
    filewriter.registerWrite("#!dendrimer_branching_point_functionality", new WriterDendrimerBranchingPointFunctionality<IngredientsType>(source));

}

#endif /*FEATURE_SYSTEM_INFORMATION_TENDOMER_H*/
