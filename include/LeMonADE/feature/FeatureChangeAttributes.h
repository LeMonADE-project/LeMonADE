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

#ifndef LEMONADE_FEATURE_FEATURECHANGEATTRIBUTES_H
#define LEMONADE_FEATURE_FEATURECHANGEATTRIBUTES_H


#include <LeMonADE/feature/FeatureAttributes.h>

/*****************************************************************/
/**
 * @class WriteChangeAttributes
 *
 * @brief Handles BFM-File-Write \b !changeattributes
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType, class TagType>
class WriteChangeAttributes:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !attributes into the header of the bfm-file.
  WriteChangeAttributes(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}
  virtual ~WriteChangeAttributes(){}
  virtual void writeStream(std::ostream& strm);
};


/*****************************************************************/
/**
 * @class FeatureChangeAttributes
 * @brief Extends vertex/monomer by an attribute tag (MonomerAttributeTag) and provides read/write functionality as well as the motivication during the simulation.
 **/
template<class TagType=int32_t>
class FeatureChangeAttributes:public FeatureAttributes<TagType>
{
public:
  //! This Feature requires a monomer_extensions.
  typedef LOKI_TYPELIST_1(MonomerAttributeTag<TagType>) monomer_extensions;
  
  typedef LOKI_TYPELIST_1(FeatureAttributes<TagType>) required_features_back;

  //! Export the relevant functionality for reading bfm-files to the responsible reader object
  template<class IngredientsType>
  void exportRead(FileImport<IngredientsType>& fileReader);

  //! Export the relevant functionality for writing bfm-files to the responsible writer object
  template<class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

};



/******************************************************************************
 * member implementations
 * ****************************************************************************/



/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !changeattributes
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class TagType>
template<class IngredientsType>
void FeatureChangeAttributes<TagType>::exportRead(FileImport< IngredientsType >& fileReader)
{
  auto  reader = new ReadAttributes<IngredientsType,TagType>(fileReader.getDestination());
//   fileReader.registerRead("!attributes",reader);
  fileReader.registerRead("!changeattributes",reader);
}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !changeattributes
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class TagType>
template<class IngredientsType>
void FeatureChangeAttributes<TagType>::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
//   fileWriter.registerWrite("!attributes",new WriteAttributes<IngredientsType>(fileWriter.getIngredients_()));
  fileWriter.registerWrite("!changeattributes",new WriteChangeAttributes<IngredientsType,TagType>(fileWriter.getIngredients_()));
}




//! Executes the routine to write \b !changeattributes.
template < class IngredientsType, class TagType>
void WriteChangeAttributes<IngredientsType,TagType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"!changeattributes\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers=molecules.size();
  //attribute blocks begin with startIndex
  size_t startIndex=0;
  //counter varable
  size_t n=0;
  //attribute to be written (updated in loop below)
  TagType attribute=molecules[0].getAttributeTag();

  //write attibutes (blockwise)
  while(n<nMonomers){
    if(molecules[n].getAttributeTag()!=attribute)
    {
      strm<<startIndex+1<<"-"<<n<<":"<<attribute<<std::endl;
      attribute=molecules[n].getAttributeTag();
      startIndex=n;
    }
    n++;
  }
  //write final changed attributes
  strm<<startIndex+1<<"-"<<nMonomers<<":"<<attribute<<std::endl<<std::endl;

}

#endif /* LEMONADE_FEATURE_FEATUREATTRIBUTES_H */

