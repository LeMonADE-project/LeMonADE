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

#ifndef LEMONADE_FEATURE_FEATUREFIXEDMONOMERS_H
#define LEMONADE_FEATURE_FEATUREFIXEDMONOMERS_H

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/feature/Feature.h>

#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

/**
 * @file
 * @brief All classes used by FeatureFixedMonomers (FeatureFixedMonomers,MonomerMovableTag,ReadFixedMonomers,WriteFixedMonomers)
 * */


/******************************************************************************
 * class definitions
 * ****************************************************************************/

/**
 * @class MonomerMovableTag
 * @brief Extends monomers by an boolean tag for fixed (false) or movable (true) properties along with getter and setter.
 *
 * @details Initially the movable tag is set to true (movable).
 * This follows the convention in Move and checkMove() that allowed moves are true.
 **/
class MonomerMovableTag
{
public:
  //! constructor sets initial tag to true==moveable
  MonomerMovableTag():movable(true){}

  /**
   * @brief Getting the movable tag of the monomer.
   *
   * @return True if movable. False if fixed.
   */
  bool getMovableTag() const {return movable;}

  /**
   * @brief Setting the movable tag of the monomer with _movable.
   *
   * @param _movable True if movable. False if fixed.
   */
  void setMovableTag(bool _movable){ movable=_movable;}

private:

  //! Private tag if monomer is movable (true)
  bool movable;
};


/**
 * @class ReadFixedMonomers
 *
 * @brief Handles BFM-File-Reads \b #!fixed_monomers.
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadFixedMonomers: public ReadToDestination<IngredientsType>
{
public:
  ReadFixedMonomers(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadFixedMonomers(){}
  virtual void execute();
};

/**
 * @class WriteFixedMonomers
 *
 * @brief Handles BFM-File-Write \b #!fixed_monomers
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteFixedMonomers:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b #!fixed_monomers into the header of the bfm-file.
  WriteFixedMonomers(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
  virtual ~WriteFixedMonomers(){}
  virtual void writeStream(std::ostream& strm);
};

/**
 * @class FeatureFixedMonomers
 *
 * @brief Extends vertex/monomer by an movable tag (MonomerMovableTag) and provides read/write functionality.
 **/
class FeatureFixedMonomers:public Feature
{
public:
  typedef LOKI_TYPELIST_1(MonomerMovableTag) monomer_extensions;

  //! Export the relevant functionality for reading bfm-files to the responsible reader object
  template<class IngredientsType>
  void exportRead(FileImport<IngredientsType>& fileReader);

  //! Export the relevant functionality for writing bfm-files to the responsible writer object
  template<class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

  /**
   * @brief Check move for all unknown moves: this does nothing
   *
   * @details Returns true for all moves other than the ones that have specialized versions of this function.
   * This dummy function is implemented for generality.
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system
   * @param [in] move General move other than MoveLocalBase (MoveLocalSc or MoveLocalBcc).
   * @return true Always!
   */
  template<class IngredientsType>
  bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const
  {
    return true;
  }

  /**
   * @brief Overloaded for MoveLocalBase. See MoveLocalSc and MoveLocalBcc
   *
   * @details Checks if the monomer is able to move and not fixed.
   * Returns \a True if monomer is movable (\a true ) or fixed (\a false ).
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system.
   * @param [in] move A reference to MoveLocalBase.
   * @return if monomer is movable (true) or fixed (false).
   */
  template<class IngredientsType,class LocalMoveType>
  bool checkMove(const IngredientsType& ingredients, const MoveLocalBase<LocalMoveType>& move) const
  {
    //get the number of bond partners of the particle to be moved
    uint32_t monoIndex=move.getIndex();
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    return molecules[monoIndex].getMovableTag();
  }

};


/******************************************************************************
 * member implementations
 * ****************************************************************************/

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!fixed_monomers
 *
 * @param fileReader File importer for the bfm-file
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureFixedMonomers::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("#!fixed_monomers",new ReadFixedMonomers<IngredientsType>(fileReader.getDestination()));
}


/******************************************************************************/
/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * #!fixed_monomers
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureFixedMonomers::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("#!fixed_monomers",new WriteFixedMonomers<IngredientsType>(fileWriter.getIngredients_()));
}

/******************************************************************************/
/**
 * @brief Executes the reading routine to extract \b #!fixed_monomers.
 *
 * @throw <std::runtime_error> attributes and identifier could not be read.
 **/
template < class IngredientsType>
void ReadFixedMonomers<IngredientsType>::execute()
{
	std::cout << " executeReadFixedMonomers"<< std::endl;

  //some variables used during reading
  //counts the number of fixedMonomers lines in the file
  int nfixedMonomers=0;
  int startIndex,stopIndex,isfixedMonomers;
  //contains the latest line read from file
  std::string line;
  //used to reset the position of the get pointer after processing the command
  std::streampos previous;
  //for convenience: get the input stream
  std::istream& source=this->getInputStream();
  //for convenience: get the set of monomers
  typename IngredientsType::molecules_type& molecules=this->getDestination().modifyMolecules();

  //go to next line and save the position of the get pointer into streampos previous
  getline(source,line);
  previous=(source).tellg();

  //read and process the lines containing the bond vector definition
  getline(source,line);

  while(!line.empty() && !((source).fail())){

    //stop at next Read and set the get-pointer to the position before the Read
    if(this->detectRead(line)){
      (source).seekg(previous);
      break;
    }

    //initialize stringstream with content for ease of processing
    std::stringstream stream(line);

    //read vector components
    stream>>startIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
      std::stringstream messagestream;
      messagestream<<"ReadFixedMonomers<IngredientsType>::execute()\n"
		   <<"Could not read first index in fixed_monomers line "<<nfixedMonomers+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

      std::stringstream messagestream;
      messagestream<<"ReadFixedMonomers<IngredientsType>::execute()\n"
		   <<"Wrong definition of fixed_monomers\nCould not find separator \"-\" "
		   <<"in fixed_monomers definition no "<<nfixedMonomers+1;
      throw std::runtime_error(messagestream.str());

    }

    //read bond identifier, throw exception if extraction fails
    stream>>stopIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
      std::stringstream messagestream;
      messagestream<<"ReadFixedMonomers<IngredientsType>::execute()\n"
		   <<"Could not read second index in fixed_monomers line "<<nfixedMonomers+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt ":"
    if(!this->findSeparator(stream,':')){

      std::stringstream messagestream;
      messagestream<<"ReadFixedMonomers<IngredientsType>::execute()\n"
		   <<"Wrong definition of fixed_monomers\nCould not find separator \":\" "
		   <<"in fixed_monomers definition no "<<nfixedMonomers+1;
      throw std::runtime_error(messagestream.str());

    }
    //read the isfixedMonomers tag
    stream>>isfixedMonomers;
    //if extraction worked, save the attributes
    if(!stream.fail()){

      //save attributes
      for(int n=startIndex;n<=stopIndex;n++)
      {
	//use n-1 as index, because bfm-files start counting indices at 1 (not 0)
	molecules[n-1].setMovableTag((isfixedMonomers == 0));
	//std::cout << " monomer: " << (n-1) << "  movable: " << (attribute == 0) << "  vs  " << molecules[n-1].getMovableTag() << std::endl;

      }
      nfixedMonomers++;
      getline((source),line);

    }
    //otherwise throw an exception
    else{

      std::stringstream messagestream;
      messagestream<<"ReadFixedMonomers<IngredientsType>::execute()\n"
		   <<"could not read fixed_monomers in fixed_monomers definition no "<<nfixedMonomers+1;
      throw std::runtime_error(messagestream.str());

    }
  }
}

/******************************************************************************/
//! Executes the routine to write \b #!fixed_monomers.
template < class IngredientsType>
void WriteFixedMonomers<IngredientsType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"#!fixed_monomers\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers=molecules.size();
  //attribute blocks begin with startIndex
  size_t startIndex=0;
  //counter varable
  size_t n=0;
  //attribute to be written (updated in loop below)
  int isfixedMonomers=molecules[0].getMovableTag() ? 0 : 1; //0->movable, 1-> fixed

  //write attibutes (blockwise)
  while(n<nMonomers){
    if((molecules[n].getMovableTag() ? 0 : 1)!=isfixedMonomers)
    {
      strm<<startIndex+1<<"-"<<n<<":"<<isfixedMonomers<<std::endl;
      isfixedMonomers=molecules[n].getMovableTag() ? 0 : 1; //0->movable, 1-> fixed
      startIndex=n;
    }
    n++;
  }
  //write final fixed_monomers: 0->movable, 1-> fixed
  strm<<startIndex+1<<"-"<<nMonomers<<":"<<isfixedMonomers<<std::endl<<std::endl;

}


#endif /* LEMONADE_FEATURE_FEATUREFIXEDMONOMERS_H */
