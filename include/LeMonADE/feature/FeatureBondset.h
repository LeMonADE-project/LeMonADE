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

#ifndef LEMONADE_FEATURE_FEATUREBONDSET_H
#define LEMONADE_FEATURE_FEATUREBONDSET_H

#include <iostream>
#include <sstream>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/utility/FastBondset.h>
#include <LeMonADE/utility/SlowBondset.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>



/*****************************************************************/
/**
 * @class ReadBondset
 *
 * @brief Handles BFM-File-Reads \b !set_of_bondvectors
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template <class IngredientsType>
class ReadBondset : public AbstractRead
{
public:
  ReadBondset(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};


/*****************************************************************/
/**
 * @class WriteBondset
 *
 * @brief Handles BFM-File-Write \b !set_of_bondvectors
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBondset : public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !set_of_bondvectors into the header of the bfm-file.
  WriteBondset(const IngredientsType& system):AbstractWrite<IngredientsType>(system){this->setHeaderOnly(true);}
  void writeStream(std::ostream& strm);
};


/*****************************************************************/
/**
 * @brief Feature adding a set of bond-vectors (Bondset or SlowBondset) to the system.
 *
 **/
template< class BondSetType=FastBondset>
class FeatureBondset : public Feature
{
 public:
	//! Standard constructor (empty)
  FeatureBondset(){}

  //! Standard destructor (empty)
  virtual ~FeatureBondset(){}

  //! Get read/write-reference to the set of bond-vectors for modifications.
  BondSetType& modifyBondset()     {return bondset;};

  //! Get read-only reference to the set of bond-vectors.
  const BondSetType&    getBondset()const{return bondset;};

  /* if you want to use the SlowBondset, comment the two lines above and uncomment
   * the two lines below.*/
// 	SlowBondset& modifyBondset()     {return bondset;};
//  const SlowBondset&    getBondset()const{return bondset;};



 /**
  * @brief Export the relevant functionality for reading bfm-files to the responsible reader object
  *
  * @details The function is called by the Ingredients class when an object of type Ingredients
  * is associated with an object of type FileImport. The export of the Reads is thus
  * taken care automatically when it becomes necessary.\n
  * Registered Read-In Commands:
  * * !set_of_bondvectors
  *
  * @param fileReader File importer for the bfm-file
  * @tparam IngredientsType Features used in the system. See Ingredients.
  **/
  template <class IngredientsType>
  void exportRead(FileImport <IngredientsType>& fileReader)
  {
    fileReader.registerRead("!set_of_bondvectors",new ReadBondset <FeatureBondset> (*this));
  };


 /**
  * @brief Export the relevant functionality for writing bfm-files to the responsible writer object
  *
  * @details The function is called by the Ingredients class when an object of type Ingredients
  * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
  * taken care automatically when it becomes necessary.\n
  * Registered Write-Out Commands:
  * * !set_of_bondvectors
  *
  * @param fileWriter File writer for the bfm-file.
  * @tparam IngredientsType Features used in the system. See Ingredients.
  */
  template <class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
  {
    fileWriter.registerWrite("!set_of_bondvectors", new WriteBondset<FeatureBondset>(*this));
  }

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
   * @details Checks if the new bond for this move of type LocalMoveType is valid.
   * Returns if move is allowed (\a true ) or rejected (\a false ).
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system.
   * @param [in] move A reference to LocalMoveType.
   * @return if move is allowed (true) or rejected (false).
   */
  template<class IngredientsType,class LocalMoveType>
  bool checkMove(const IngredientsType& ingredients, const MoveLocalBase<LocalMoveType>& move) const
  {

      //get the number of bond partners of the particle to be moved
          uint32_t monoIndex=move.getIndex();
          const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

          for (size_t j=0; j< molecules.getNumLinks(monoIndex); ++j){
              if (!bondset.isValidStrongCheck(molecules[molecules.getNeighborIdx(monoIndex,j)]-(molecules[monoIndex]+move.getDir()))) return false;
          }

          return true;
  }
  
  /**
   * @brief Overloaded for MoveConnectBase. 
   *
   * @details Checks if the new bond for this move of type ConnectMoveType is valid.
   * Returns if move is allowed (\a true ) or rejected (\a false ).
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system.
   * @param [in] move A reference to ConnectMoveType.
   * @return if move is allowed (true) or rejected (false).
   */
  template<class IngredientsType,class ConnectMoveType>
  bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<ConnectMoveType>& move) const
  {

	  //get the number of bond partners of the particle to be moved
          uint32_t MonID=move.getIndex();
	  uint32_t partnerID=move.getPartner();
          const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

	  if (!bondset.isValidStrongCheck(molecules[MonID]-molecules[partnerID])) return false;

          return true;
  }
  /**
   * @brief Updates the bond-set lookup table if necessary
   *
   * @details Synchronizes the bond-set look-up table with the rest of the system, and checks for invalid bonds.
   *
   * @param ingredients A reference to the IngredientsType - mainly the system.
   */
  template<class IngredientsType> void synchronize(IngredientsType& ingredients)
  {
    //this function only does something if the bondset has change since the last update
    bondset.updateLookupTable();

    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

    for (size_t i=0; i< molecules.size(); ++i)
    {
     for (size_t j=0; j< molecules.getNumLinks(i); ++j){

	 uint n = molecules.getNeighborIdx(i,j);
	 if (!bondset.isValidStrongCheck(molecules[n]-molecules[i]))
	{
	  std::ostringstream errorMessage;
	  errorMessage << "FeatureBondset::synchronize(): Invalid bond vector between monomer " << i << " at " << molecules[i] << " and " << n << " at " <<  molecules[n] <<  ".\n";throw std::runtime_error(errorMessage.str());
	}
    }
    }
  }

protected:

  //! Stores the set of allowed bond-vectors.
  BondSetType bondset;


  /* if you want to use the SlowBondset, comment the line above and uncomment
   * the line below.*/
  //SlowBondset bondset;
};


/**
 * @brief Executes the reading routine to extract \b !set_of_bondvectors.
 *
 * @throw <std::runtime_error> bond-vector and identifier could not be read.
 **/
template <class LemonadeSystem>
void ReadBondset <LemonadeSystem>::execute()
{
  int nBondVectors=0;
  int x,y,z,identifier;
  std::string line;
  std::streampos previous;

  //go to next line and save the position of the get pointer into streampos previous
  getline(*source,line);
  previous=(*source).tellg();

  //read and process the lines containing the bond vector definition
  getline(*source,line);

  while(!line.empty() && !((*source).fail())){

    //stop at next Read and set the get-pointer to the position before the Read
    if(detectRead(line)){
      (*source).seekg(previous);
      break;
    }

    //initialize stringstream with content for ease of processing
    std::stringstream stream(line);

    //read vector components
    stream>>x>>y>>z;

    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadBondset<LemonadeSystem>::execute()\n"
		   <<"Could not read vector components in vector no "<<nBondVectors+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt ":"
    if(!findSeparator(stream,':')){

    	std::stringstream messagestream;
      messagestream<<"ReadBondset<LemonadeSystem>::execute()\n"
		   <<"Wrong definition of bondvector\nCould not find separator \":\" ";
      throw std::runtime_error(messagestream.str());

    }

    //read bond identifier, throw exception if extraction fails
    stream>>identifier;

    if(!stream.fail()){

      //add bond to bondset
      bfmSystem.modifyBondset().addBond(x,y,z,identifier);
      nBondVectors++;
      getline((*source),line);

    }
    else{

    	std::stringstream messagestream;
      messagestream<<"ReadBondset<LemonadeSystem>::execute()\n"
		   <<"could not read identifier in vector no "<<nBondVectors+1;
      throw std::runtime_error(messagestream.str());

    }
  }

}


//! Executes the routine to write \b !set_of_bondvectors.
template<class LemonadeSystem>
void WriteBondset<LemonadeSystem>::writeStream(std::ostream& strm){
  strm<<"!set_of_bondvectors\n";
  std::map <int32_t, VectorInt3>::const_iterator bondVec;

  //iterate over all bonds in bondset and print to stream, using point's operator <<
  for(bondVec=this->getSource().getBondset().begin();bondVec!=this->getSource().getBondset().end();++bondVec){
    strm<<(bondVec->second)<<":"<<(bondVec->first)<<"\n";
  }
  strm<<"\n";

}

#endif
