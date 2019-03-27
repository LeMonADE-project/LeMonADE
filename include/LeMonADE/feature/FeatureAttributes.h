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

#ifndef LEMONADE_FEATURE_FEATUREATTRIBUTES_H
#define LEMONADE_FEATURE_FEATUREATTRIBUTES_H

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBase.h>

/**
 * @file
 * @brief Adding an attribute tag to each monomer.
 * */

/**
 * @class MonomerAttributeTag
 * @brief Extends monomers by an signed integer (int32_t) as tag along with getter and setter\n
 * 		  Initially the tag is set to NULL.
 * */
template<class TagType=int32_t>
class MonomerAttributeTag
{
public:

	//! Standard constructor- initially the tag is set to NULL.
	MonomerAttributeTag():tag(){}

	//! Getting the tag of the monomer.
	TagType getAttributeTag() const {return tag;}

	/**
	 * @brief Setting the tag of the monomer with \para attr.
	 *
	 * @param attr
	 */
	void setAttributeTag(TagType attribute){ tag=attribute;}

private:
	 //! Private variable holding the tag. Default is NULL.
     TagType tag;
};


/*****************************************************************/
/**
 * @class ReadAttributes
 *
 * @brief Handles BFM-File-Reads \b !attributes
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType, class TagType>
// template < class TagType> 
class ReadAttributes: public ReadToDestination<IngredientsType>
{
public:
  ReadAttributes(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadAttributes(){}
  virtual void execute();
};


/*****************************************************************/
/**
 * @class WriteAttributes
 *
 * @brief Handles BFM-File-Write \b !attributes
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType, class TagType>
class WriteAttributes:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !attributes into the header of the bfm-file.
  WriteAttributes(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
  virtual ~WriteAttributes(){}
  virtual void writeStream(std::ostream& strm);
};

/*****************************************************************/
/**
 * @class FeatureAttributes
 * @brief Extends vertex/monomer by an attribute tag (MonomerAttributeTag) and provides read/write functionality.
 **/
template<class TagType=int32_t>
class FeatureAttributes:public Feature
{
public:
  //! This Feature requires a monomer_extensions.
  typedef LOKI_TYPELIST_1(MonomerAttributeTag<TagType>) monomer_extensions;

  //! Export the relevant functionality for reading bfm-files to the responsible reader object
  template<class IngredientsType>
  void exportRead(FileImport<IngredientsType>& fileReader);

  //! Export the relevant functionality for writing bfm-files to the responsible writer object
  template<class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

  //! For all unknown moves: this does nothing
  template<class IngredientsType>
  void applyMove(IngredientsType& ing, const MoveBase& move){};


  //! Overloaded for moves of type MoveAddMonomerBase to set the attribute tag by inserting a monomer
  template<class IngredientsType,class AddMoveType>
  void applyMove(IngredientsType& ing, const MoveAddMonomerBase<AddMoveType, TagType>& move);

};


/******************************************************************************
 * member implementations
 * ****************************************************************************/



/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !attributes
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class TagType>
template<class IngredientsType>
void FeatureAttributes<TagType>::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("!attributes",new ReadAttributes<IngredientsType,TagType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !attributes
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class TagType>
template<class IngredientsType>
void FeatureAttributes<TagType>::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("!attributes",new WriteAttributes<IngredientsType,TagType>(fileWriter.getIngredients_()));
}


/**
 * This function updates the attribute tag given by move itself (most cases Zero).
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move general addmove (MoveAddScMonomer/MoveAddBccMonomer)
 */
template<class TagType>
template<class IngredientsType,class AddMoveType>
void FeatureAttributes<TagType>::applyMove(IngredientsType& ingredients, const MoveAddMonomerBase<AddMoveType, TagType>& move)
{
	ingredients.modifyMolecules()[move.getMonomerIndex()].setAttributeTag(move.getTag());
}


/**
 * @brief Executes the reading routine to extract \b !attributes.
 *
 * @throw <std::runtime_error> attributes and identifier could not be read.
 **/
template < class IngredientsType, class TagType >
// template <class TagType >
void ReadAttributes<IngredientsType,TagType>::execute()
{
  //some variables used during reading
  //counts the number of attribute lines in the file
  int nAttributes=0;
  int startIndex,stopIndex;
  TagType attribute;
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
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read first index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \"-\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }

    //read bond identifier, throw exception if extraction fails
    stream>>stopIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read second index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt ":"
    if(!this->findSeparator(stream,':')){

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \":\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
    //read the attribute tag
    stream>>attribute;
    //if extraction worked, save the attributes
    if(!stream.fail()){

      //save attributes
      for(int n=startIndex;n<=stopIndex;n++)
      {
	//use n-1 as index, because bfm-files start counting indices at 1 (not 0)
	molecules[n-1].setAttributeTag(attribute);
      }
      nAttributes++;
      getline((source),line);

    }
    //otherwise throw an exception
    else{

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"could not read attribute in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
  }
}


//! Executes the routine to write \b !attributes.
template < class IngredientsType, class TagType>
void WriteAttributes<IngredientsType,TagType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"!attributes\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers=molecules.size();
  //attribute blocks begin with startIndex
  size_t startIndex=0;
  //counter varable
  size_t n=0;
  //attribute to be written (updated in loop below)
  TagType attribute=molecules[0].getAttributeTag();

  //write attribute (blockwise)
  while(n<nMonomers){
    if(molecules[n].getAttributeTag()!=attribute)
    {
      strm<<startIndex+1<<"-"<<n<<":"<<attribute<<std::endl;
      attribute=molecules[n].getAttributeTag();
      startIndex=n;
    }
    n++;
  }
  //write final attributes
  strm<<startIndex+1<<"-"<<nMonomers<<":"<<attribute<<std::endl<<std::endl;

}


#endif /* LEMONADE_FEATURE_FEATUREATTRIBUTES_H */
