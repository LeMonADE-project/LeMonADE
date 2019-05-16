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

#ifndef LEMONADE_CORE_MOLECULESREAD_H
#define LEMONADE_CORE_MOLECULESREAD_H

/**
 * @file 
 * @brief Reading routines for bfm-Reads !number_of_monomers, 
 * !bonds/!add_bonds, !remove_bonds and !mcs 
 **/

#include <sstream>
#include <limits>
#include <list>
#include <iostream>


#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractRead.h>


/***********************************************************************/
/**
 * @class ReadNrOfMonomers
 * @brief Handles BFM-File-Reads \b !number_of_monomers
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadNrOfMonomers : public ReadToDestination < IngredientsType >
{
  public:
	//! Empty Constructor, but delegates call to the Feature
    ReadNrOfMonomers(IngredientsType& destination):ReadToDestination< IngredientsType > (destination){}

    void execute();
};

/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !number_of_monomers.
 *
 * @throw <std::runtime_error> if number of monomers could not be read.
 **/
template < class IngredientsType >
void ReadNrOfMonomers<IngredientsType>::execute()
{
  std::cout<<"reading number of monomers...";
  int x;
  if( this->getInputStream()>> x)
  {
     this->getDestination().modifyMolecules().resize(x);
     std::cout<<x<<std::endl;
  }
  else throw std::runtime_error("ReadNrOfMonomers<IngredientsType>::execute()\n Could not read number of monomers");
}


/***********************************************************************/
/**
 * @class ReadBonds
 * @brief Handles BFM-File-Reads \b !bonds.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadBonds: public ReadToDestination < IngredientsType >
{
 public:
	//! Empty Constructor, but delegates call to the Feature
	ReadBonds(IngredientsType& destination):ReadToDestination< IngredientsType > (destination){}

	void execute();
};

/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !bonds.
 *
 * @throw <std::runtime_error> if bond could not be read.
 **/
template < class IngredientsType >
void ReadBonds<IngredientsType>::execute()
{

  int a,b;
  int nBonds=0;

  std::string line;
  std::streampos previous;
  //go to next line and save position of get pointer to streampos previous
  this->getInputStream().get();
  previous=this->getInputStream().tellg();

  //get line from file for processing
  getline(this->getInputStream(),line);

  //start processing the line
  while(!line.empty() && this->getInputStream().good()){

    //if the line contains a bfm Read, stop the procedure and set the get pointer back
    if(this->detectRead(line)){
      this->getInputStream().seekg(previous);
      break;
    }

    //initialize a stringstream with the line for easier processing
    std::stringstream stream(line);

    //read bond partners
    stream>>a>>b;
    //throw exception if something went wrong
    if(stream.fail()){

    	std::stringstream errormessage;
      errormessage<<"ReadBonds<IngredientsType>::execute()\n"
		  <<"Could not read bond partners in "<<nBonds+1<<" th bond definition\n";
      throw std::runtime_error(errormessage.str());

    }

    //if still here, add bond to bondset
    this->getDestination().modifyMolecules().connect(a-1,b-1);
    nBonds++;
    //read next line from file
    getline(this->getInputStream(),line);
  }

}

/***********************************************************************/
/**
 * @class ReadRemoveBonds
 * @brief Handles BFM-File-Reads \b !remove_bonds.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadRemoveBonds: public ReadToDestination < IngredientsType >
{
 public:
	//! Empty Constructor, but delegates call to the Feature
	ReadRemoveBonds(IngredientsType& destination):ReadToDestination< IngredientsType > (destination){}

	void execute();
};
/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !remove_bonds.
 *
 * @throw <std::runtime_error> if broken bond could not be read.
 **/

template < class IngredientsType >
void ReadRemoveBonds<IngredientsType>::execute()
{
  int a,b;
  int nBreaks=0;
  
  std::string line;
  std::streampos previous;
  //go to next line and save position of get pointer to streampos previous
  this->getInputStream().get();
  previous=this->getInputStream().tellg();

  //get line from file for processing
  getline(this->getInputStream(),line);

  //start processing the line
  while(!line.empty() && this->getInputStream().good()){

    //if the line contains a bfm Read, stop the procedure and set the get pointer back
    if(this->detectRead(line)){
      this->getInputStream().seekg(previous);
      break;
    }

    //initialize a stringstream with the line for easier processing
    std::stringstream stream(line);

    //read bond partners
    stream>>a>>b;
    //throw exception if something went wrong
    if(stream.fail()){

    	std::stringstream errormessage;
      errormessage<<"ReadRemoveBonds<IngredientsType>::execute()\n"
		  <<"Could not read broken bond partners in "<<nBreaks+1<<" th broken bond definition\n";
      throw std::runtime_error(errormessage.str());

    }
      
    //if still here, add broken bond to bondset
    this->getDestination().modifyMolecules().disconnect(a-1,b-1);
    nBreaks++;
    //read next line from file
    getline(this->getInputStream(),line);
  }

}

/***********************************************************************/
/**
 * @class ReadMcs
 * @brief Handles BFM-File-Reads \b !mcs.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadMcs: public ReadToDestination < IngredientsType >
{
	//! counts the number of monomers in each mcs in the file to check consistency
	int monomerCount;
	//! Holds if first !mcs is read-in, already.
	bool first_call;

	//! Ignores if monomer position is missing in the !mcs (maybe due to solvent)
	bool ignore_missing_monos;

	//! The number of frames/conformations/!mcs-command in file.
	uint32_t nFrames;

	/**
	 * @brief Return if first !mcs is read-in, already.
	 *
	 * @return \a True if first call/occurance of !mcs - \a False otherwise.
	 */
	bool isFirstCall() const {return first_call;}

	//! processes a regular mcs line from a stream
	void processRegularLine(const std::string& line);

	//! processes a line beginning with keyword solvent
	void processSolventLine(const std::string& line,uint32_t& offset);

	//! processes a line beginning with keyword solvent
	void processSolventContinueLine(const std::string& line,uint32_t& offset);


 public:
  //! Empty Constructor, but delegates call to the Feature. Default: no ignoring of missing monomers.
  ReadMcs(IngredientsType& destination):ReadToDestination< IngredientsType > (destination)
  ,monomerCount(0)
  ,first_call(true)
  ,ignore_missing_monos(false)
  ,nFrames(0)
  {}

  /**
   * @brief Set this if missing monomers should ignore by read-in (maybe due to solvent).
   *
   * @details Specify the variable \a ignore_missing_monos if you need this (dangerous) behavior.
   * If not you get a \a runtime_error if monomers are missing
   * @param val True if missing monomers should ignored - False otherwise.
   */
  void setIgnoreMissingMonomers(bool val){ignore_missing_monos = val;}

  void execute();

};

/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !mcs.
 *
 * @throw <std::runtime_error> if Monte-Carlo-Step could not be parsed or if monomers are missing.
 **/
template < class IngredientsType > void ReadMcs< IngredientsType >::execute()
{
  //get writeable reference to the molecules container
  typename IngredientsType::molecules_type& molecules = this->getDestination().modifyMolecules();


  if ( this->isFirstCall() ) {
	std::cout << "ReadMcs:execute() : updating ";
	std::cout << "connectivity and ";
	std::cout << "positions.";
  }


  //some variables needed for reading
  unsigned long mcs;
  std::string line;
  std::streampos previous;

  //read mcs number from file and update data
  this->getInputStream()>>mcs;

  molecules.setAge(mcs);

  //throw exception if there is a reading error or update data otherwise
  if(this->getInputStream().fail()){

    std::stringstream errormessage;
    errormessage<<"ReadMcs<IngredientsType>::execute()\n"
		<<"Could not read mcs number. Previous mcs number was "<<molecules.getAge();
    throw std::runtime_error(errormessage.str());

  }

  molecules.setAge(mcs);

  nFrames++;
  if(nFrames%1000 == 0) std::cout<<"Setting age: "<<mcs<<std::endl;

  // the following lines will check if 'jumps' of the '!mcs'-Read are present
  // if so you should use the FeatureJumps for backward-compability
  // for the new LeMonADe-file there´re no jumps anymore!
  std::string jumpString;
  getline(this->getInputStream(),jumpString);

  //tokenize the jump-string
  std::vector<std::string> jumpStringVector;
  jumpStringVector = this->tokenize2Parameter(jumpString, ' ', ':');

  // if size of jumpStringVector is != 0 it seems there are jumps in the file
  if(jumpStringVector.size() != 0)
  {
	  std::stringstream errormessage;
	  errormessage<<"ReadMcs<IngredientsType>::execute()\n"
	 		<<"It seems this is an old bfm-file with jumps.\n"
	 		<<"Please use FeatureJumps for read-in backward-functionality.\n"
	 		<<"In LeMonADe-files there´re no jumps anymore!\n"
	 		<<"Jumps at !mcs: "<<molecules.getAge();
      throw std::runtime_error(errormessage.str());
  }


  //go on with the all positions etc.
  previous=this->getInputStream().tellg();
  getline(this->getInputStream(),line); //get the first line with positions

  //process input lines in this loop
  while(!line.empty() && !this->getInputStream().fail()){

	  //if the line contains a bfm Read, stop the procedure and set the get pointer back
	  if(this->detectRead(line)){
		  this->getInputStream().seekg(previous);
		  return;
	  }

	  //if the line starts with the solvent keyword, process the compressed solvent format
	  if(line.rfind("solvent ",0)==0)
	  {
		  size_t startMonomerIndex=monomerCount;
		  //the offset counts the current position in a linearized array of box coordinates
		  uint32_t offset=0;
		  processSolventLine(line,offset);
		  getline(this->getInputStream(),line);

		  //if solvent extends over more than one line, process "sc" lines as well
		  while(line.rfind("sc ",0)==0)
		  {
			  processSolventContinueLine(line,offset);
			  getline(this->getInputStream(),line);
		  }
		  size_t stopMonomerIndex=monomerCount;

		  //remember which monomers are to be compressed when writing an output file
		  this->getDestination().setCompressedOutputIndices(startMonomerIndex,stopMonomerIndex-1);
	  }
	  else
	  {
		  //process this line and get the next one from the file
		  processRegularLine(line);
		  getline(this->getInputStream(),line);
	  }
  }

  //if total number of monomers in !mcs differs from previous, throw exception
  if((uint32_t)monomerCount!=molecules.size()){

    std::stringstream errormessage;
    errormessage<<"ReadMcs<IngredientsType>::execute()\n"
		<<"Inconsistent number of monomers in bfm file in mcs number "
		<<molecules.getAge()<<std::endl;

    if ( (this->ignore_missing_monos==true)  || this->isFirstCall() )
    {
		if(nFrames%1000==0)
		{
			std::cout << "Ignoring logic error:\n" << errormessage.str() << std::endl;
			std::cout << "Not all monomers have been initialized.\n" << std::endl;
		}
    }
    else
    {
		throw std::runtime_error(errormessage.str());
    }

  }

  monomerCount=0;
  first_call = false;


}

//! reads coordinates from a line containing a connected chain
template<class IngredientsType>void ReadMcs<IngredientsType>::processRegularLine(const std::string& line )
{
	//put content into stream for ease of processing
	std::stringstream stream(line);
	typename IngredientsType::molecules_type& molecules = this->getDestination().modifyMolecules();
	int x,y,z,bondId;
	//read first three coordinates
	stream>>x>>y>>z;


	if(!stream.fail())
	{
		molecules[monomerCount].setAllCoordinates(x,y,z);
			++monomerCount;
	}

	//throw exception if first coordinates of chain cannot be extracted from file
	if(stream.fail())
	{

		std::stringstream errormessage;
		errormessage<<"ParsingError::ReadMcs<IngredientsType>::execute()\n"
				<<"Could not read chain initial coordinates in mcs "<<molecules.getAge();

		throw std::runtime_error(errormessage.str());

	}

	//ignore spaces
	stream.ignore(1);

	//read the ASCII coded bond vectors of this chain
	while(stream.peek()>0 && stream.good())
	{
		//read next bond from file
		bondId=stream.get();
		x+=this->getDestination().getBondset().getBondVector(bondId).getX();
		y+=this->getDestination().getBondset().getBondVector(bondId).getY();
		z+=this->getDestination().getBondset().getBondVector(bondId).getZ();

		//update molecules
		molecules[monomerCount].setAllCoordinates(x,y,z);


		//set up connections at first monte carlo step
		if ( this->isFirstCall())
		{
			molecules.connect(monomerCount-1,monomerCount);
		}

		//count number of monomers in this mcs
		++monomerCount;
	}

}

//! reads folded coordinates of compressed solvent from a line starting with keyword "solvent"
template<class IngredientsType> void ReadMcs<IngredientsType>::processSolventLine(const std::string& line, uint32_t& offset)
{
	typename IngredientsType::molecules_type& molecules = this->getDestination().modifyMolecules();
	//need these values to turn the linear index in the file back into positions
	uint32_t boxX=this->getDestination().getBoxX();
	uint32_t boxY=this->getDestination().getBoxY();
	int32_t x,y,z;
	//put line into stream for ease of processing
	std::stringstream stream(line);

	stream.ignore(8); //ignore "solvent "

	//now start processing the rest of the stream
	int32_t distance;
	while(stream.peek()>0 && stream.good())
	{
		//read next position
		if(stream.peek()==32) //corresponds to space
		{
			stream>>distance;
			stream.ignore(1);
		}
		else
		{
			distance=int32_t(stream.get())-33;
		}

		offset+=distance;

		//transform linear offset index back into 3d coordinates
		x=offset%boxX;
		y=int32_t(offset/boxX)%boxY;
		z=int32_t(offset/(boxX*boxY));
		molecules[monomerCount].setAllCoordinates(x,y,z);
		//count number of monomers in this mcs
		++monomerCount;
	}

}

//! reads folded coordinates of compressed solvent from a line starting with keyword "sc"
template<class IngredientsType> void ReadMcs<IngredientsType>::processSolventContinueLine(const std::string& line, uint32_t& offset)
{
	typename IngredientsType::molecules_type& molecules = this->getDestination().modifyMolecules();
	//need these values to turn the linear index in the file back into positions
	uint32_t boxX=this->getDestination().getBoxX();
	uint32_t boxY=this->getDestination().getBoxY();
	int32_t x,y,z;
	//put line into stream for ease of processing
	std::stringstream stream(line);

	stream.ignore(3); //ignore "sc "
	//in this line there is no extra offset to be read from the line
	//now start processing the rest of the stream
	int32_t distance;
	while(stream.peek()>0 && stream.good())
	{
		//read next position
		if(stream.peek()==32) //corresponds to space
		{
			stream>>distance;
			stream.ignore(1);
		}
		else
		{
			distance=int32_t(stream.get())-33;
		}

		offset+=distance;

		//transform linear offset back into 3d coordinate
		x=offset%boxX;
		y=int32_t(offset/boxX)%boxY;
		z=int32_t(offset/(boxX*boxY));
		molecules[monomerCount].setAllCoordinates(x,y,z);
		//count number of monomers in this mcs
		++monomerCount;
	}

}



#endif /* LEMONADE_CORE_MOLECULESREAD_H */
