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

#ifndef LEMONADE_FEATURE_FEATURECONNECTIONSC_H
#define LEMONADE_FEATURE_FEATURECONNECTIONSC_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>
/**
 * @class MonomerReactivity
 * @brief set the monomer reactive or unreactive 
 * @todo The Implementation of the output stream operator does not work properly. 
 */
class MonomerReactivity
{
  public:

	//! Standard constructor- initially the reactivity is set to false and default numMaxLinks is unconnected monomer.
	MonomerReactivity():reactivity(false),numMaxLinks(0){}

	//! Getting the reactivity of the monomer.
	const bool isReactive() const {return reactivity;}
	
	//! Getting the number of maximum possible bonds for the monomer.
	const uint32_t getNumMaxLinks() const {return numMaxLinks;};

	MonomerReactivity& operator= (const MonomerReactivity source)
	{
	  reactivity=source.isReactive();
	  numMaxLinks=source.getNumMaxLinks();
	  return *this;
	}
	const MonomerReactivity& getMonomerReactivity() const {return *this; }
	
	void setMonomerReactivity(const MonomerReactivity& react)
	{
	  reactivity = react.isReactive();
	  numMaxLinks = react.getNumMaxLinks();
	  /// \todo There should be a check in the function to test against the "default" maximum connectivity for consistency reason.
	}
	const bool operator == (const MonomerReactivity &react) const 
	{
	  if ( react.getNumMaxLinks() != numMaxLinks  ) return false;
	  if ( react.isReactive()   != reactivity ) return false;
	  return true;
	}
	const bool operator!= (const MonomerReactivity &react) const 
	{ 
	  return !(*this == react);
	}
	/**
	 * @brief Setting the reactivity of the monomer with \para reactivity_.
	 *
	 * @param reactivity_ either trueor false
	 */
	void setReactive(bool reactivity_){ reactivity=reactivity_;}
	
	/**
	 * @brief Setting the maximum possible bonds of the monomer with \para NumMaxLinks_.
	 *
	 * @param NumMaxLinks_
	 * \todo There should be a check in the function to test against the "default" maximum connectivity for consistency reason.
	 */	
	void setNumMaxLinks(uint32_t numMaxLinks_){numMaxLinks=numMaxLinks_;}
	
	/**
	* @brief \b Stream \b Out \b operator of the MonomerReactivity
	*
	* @details Streams out the elements chain ID, number of labels and the label
	* 	   ID  of the label separated by space
	*          its implemented inside the class to not break the one definition rule (ODR)
	* 	   (another solution would be to put the Implementation into a cpp file...)
	* @param stream output-stream
	* @param label object of class MonomerReactivity
	* @return output-stream
	**/
	friend std::istream& operator>> (std::istream& stream, MonomerReactivity & Reactivity)
	{
	  int temp;
	  stream >> temp; Reactivity.setReactive(temp); stream.ignore(1);
	  stream >> temp; Reactivity.setNumMaxLinks(temp);

	  return stream;
	}
	
	/**
	* @brief \b Stream \b Out \b operator of the MonomerReactivity
	*
	* @details Streams out the elements chain ID, number of labels and the label
	* 	     ID  of the label separated by space
	*
	* @param stream output-stream
	* @param label object of class MonomerReactivity
	* @return output-stream
	* @todo this operator does override the <<molecules[i] operator
	**/
	friend std::ostream& operator<< (std::ostream& stream, const MonomerReactivity & Reactivity)
	{
		stream
		<< Reactivity.isReactive() << "/"
		<< Reactivity.getNumMaxLinks();
		return stream;
	};
	
private:
     //! Private variable holding the tag. Default is NULL.
     bool reactivity;
     //! Number of maximimum possible links/connections to another not-yet-connected monomers
     uint32_t numMaxLinks;
};

/*****************************************************************/
/**
 * @class ReadReactivity
 *
 * @brief Handles BFM-File-Reads \b !reactivity
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType>
class ReadReactivity: public ReadToDestination<IngredientsType>
{
public:
  ReadReactivity(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadReactivity(){}
  virtual void execute();
};


/*****************************************************************/
/**
 * @class WriteReactivity
 *
 * @brief Handles BFM-File-Write \b !reactivity
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteReactivity:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !reactivity into the header of the bfm-file.
  WriteReactivity(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
  virtual ~WriteReactivity(){}
  virtual void writeStream(std::ostream& strm);
};
/*****************************************************************************/
/**
 * @file
 * @date   2019/02/05
 * @author Toni
 *
 * @class FeatureConnectionSc
 * @brief This Feature add new bonds between monomers.
 *
 * @details Works only in combination with an excluded volume feature
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureConnectionSc : public Feature {
  
public:
	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;
	
	//! this feature will not use any of FeatureExcludedVolumeSc<> but we need excluded volume property
//  	typedef LOKI_TYPELIST_1(FeatureExcludedVolumeSc<>) required_features_back;

	//constructor
	FeatureConnectionSc() :latticeFilledUp(false)
	{connectionLattice.setupLattice();}
	
	/**
	 * Returns true if the underlying lattice is synchronized and all excluded volume condition
	 * (e.g. monomer/vertex occupies lattice edges) is applied.
	 * Returns false if this feature is out-of-sync.
	 *
	 * @return true if this feature is synchronized
	 * 		   false if this feature is out-of-sync.
	 **/
	bool isLatticeFilledUp() const {
		return latticeFilledUp;
	}

	/**
	 * Set's the need of synchronization of this feature e.g. escp. if the underlying lattice needs
	 * to refilled and if all excluded volume condition needs to be updated.
	 *
	 * @param[in] latticeFilledUp Specified if ExVol should be refilled (false) or everything is in-sync (true).
	 *
	 **/
	void setLatticeFilledUp(bool latticeFilledUp) {
		this->latticeFilledUp = latticeFilledUp;
	}
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template<class IngredientsType>
	void exportRead(FileImport<IngredientsType>& fileReader);

	//! Export the relevant functionality for writing bfm-files to the responsible writer object
	template<class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const{ return true;};

	//! check bas connect move - always true 
	template<class IngredientsType, class SpecializedMove> 
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const;

	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalSc& move) const {return true;};
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalScDiag& move) const {return true;};
	
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const {throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalBcc: wrong lattice type ... \n"); return false;};
	
	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const {throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalBcc: wrong lattice type ... \n"); return false;};

	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const {return true;};
	
	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move){};
	
	//!apply move for the scBFM connection move for connection
	template<class IngredientsType, class SpecializedMove> 
	void applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove>& move);
	
	//!apply move for the scBFM local move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);	
	
	//!apply move for the scBFM local diagonal move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalScDiag& move);	

	//!
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);	
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	
	//! Get the lattice value at a certain point
	//! return monomer index between (0; molecules.size()-1) if there's a monomer
	//! return uint32_t(-1)=4294967295 if place is empty
	uint32_t getIdFromLattice(const VectorInt3& pos) const { return connectionLattice.getLatticeEntry(pos)-1;} ;

	//! Get the lattice value at a certain point
	//! return monomer index between (0; molecules.size()-1) if there's a monomer
	//! return uint32_t(-1)=4294967295 if place is empty
	uint32_t getIdFromLattice(const int x, const int y, const int z) const { return connectionLattice.getLatticeEntry(x,y,z)-1;};

protected:

	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;
	//!
	Lattice<uint32_t> connectionLattice;

};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !reactivity
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureConnectionSc::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("!reactivity",new ReadReactivity<IngredientsType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !reactivity
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureConnectionSc::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("!reactivity",new WriteReactivity<IngredientsType>(fileWriter.getIngredients_()));
}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc::checkMove( const IngredientsType& ingredients, const MoveConnectSc& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 * @details  it might make a difference for the speed if the order of statements is switched for different systems parameters
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureConnectionSc ::checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const
{
  
  	if (!latticeFilledUp)
	    throw std::runtime_error("*****FeatureConnectionSc::checkMove....lattice is not populated. Run synchronize!\n");
	uint32_t ID(move.getIndex());
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

	//check for maximum number of bonds for the first monomer
	if ( molecules.getNumLinks(ID) >=  molecules[ID].getNumMaxLinks()) return false;

	uint32_t Neighbor(move.getPartner());

	//check if neighbor is reactive 
	if ( !molecules[Neighbor].isReactive() ) return false;
	
	//check for maximum number of bonds for the second monomer
	if ( molecules.getNumLinks(Neighbor) >= molecules[Neighbor].getNumMaxLinks() ) return false;

	//check if the two monomers are already connected
	if ( molecules.areConnected(ID,Neighbor) ) return false;

	//if still here, then the two monomers are allowed to connect 
	return true;
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveConnectSc& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveConnectBase<SpecializedMove>& move)
{
  //because in the most cases we have an irreversible connection we can erase the lattice entry 
  const typename IngredientsType::molecules_type& molecules=ing.getMolecules();
  uint32_t ID(move.getIndex());
  if ( molecules.getNumLinks(ID) ==  molecules[ID].getNumMaxLinks())
    connectionLattice.setLatticeEntry(molecules[ID].getVector3D(),0);
  uint32_t Neighbor(move.getPartner());
  if ( molecules.getNumLinks(Neighbor) ==  molecules[Neighbor].getNumMaxLinks())
    connectionLattice.setLatticeEntry(molecules[Neighbor].getVector3D(),0);
//connection is made in the move 
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalSc& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveLocalSc& move)
{
  VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
  VectorInt3 direction=move.getDir();
  VectorInt3 oldPlusDir=oldPos+direction;
  connectionLattice.moveOnLattice(oldPos,oldPlusDir);
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveLocalScDiag& move)
{
  VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
  VectorInt3 direction=move.getDir();
  VectorInt3 oldPlusDir=oldPos+direction;
  connectionLattice.moveOnLattice(oldPos,oldPlusDir);
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType, class TagType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveAddMonomerSc<TagType>& move)
{
  uint32_t MonID(move.getMonomerIndex()); 
  VectorInt3 pos=ing.getMolecules()[MonID];
  if (move.isReactive()) 
    connectionLattice.setLatticeEntry(pos,MonID+1 );
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc  ::synchronize(IngredientsType& ingredients)
{

	std::cout << "FeatureConnectionSc::synchronizing lattice occupation...\n";
	fillLattice(ingredients);
	std::cout << "done\n";
}


/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc ::fillLattice(IngredientsType& ingredients)
 * @brief This function populates the lattice directly with positions from molecules.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 * */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc::fillLattice(IngredientsType& ingredients)
{
  
	connectionLattice.setupLattice(ingredients.getBoxX(),ingredients.getBoxY(),ingredients.getBoxZ());
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	//copy the lattice occupation from the monomer coordinates
	for(size_t n=0;n<molecules.size();n++)
	{
		VectorInt3 pos=ingredients.getMolecules()[n];
		if( connectionLattice.getLatticeEntry(pos)!=0 )
		{
			throw std::runtime_error("********** FeatureConnectionSc::fillLattice: multiple lattice occupation ******************");
		}
		else if (ingredients.getMolecules()[n].isReactive() && (molecules.getNumLinks(n) <  molecules[n].getNumMaxLinks()) )
		{
			//here we simply set the monomer id (plus one!) on the lattice site 
			// the offset implies that the index zero is still used for unoccupied
			// with and unreactive monomer
			VectorInt3 pos=molecules[n];
			connectionLattice.setLatticeEntry(pos,n+1);
		}
	}
	latticeFilledUp=true;
}

/**
 * @brief Executes the reading routine to extract \b !attributes.
 *
 * @throw <std::runtime_error> attributes and identifier could not be read.
 **/
template < class IngredientsType >
void ReadReactivity<IngredientsType>::execute()
{
  //some variables used during reading
  //counts the number of reactivity lines in the file
  int nReactiveBlocks=0;
  int startIndex,stopIndex;
  MonomerReactivity reactivity;
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
      messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
		   <<"Could not read first index in reactivity line "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

    	std::stringstream messagestream;
      messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
		   <<"Wrong definition of reactivity\nCould not find separator \"-\" "
		   <<"in reactivity definition no "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());

    }

    //read bond identifier, throw exception if extraction fails
    stream>>stopIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
		   <<"Could not read second index in reactivity line "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt ":"
    if(!this->findSeparator(stream,':')){

    	std::stringstream messagestream;
      messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
		   <<"Wrong definition of reactivity\nCould not find separator \":\" "
		   <<"in reactivity definition no "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());

    }
    //read the attribute tag
    stream>>reactivity;
    //if extraction worked, save the attributes
    if(!stream.fail()){

      //save attributes
      for(int n=startIndex;n<=stopIndex;n++)
      {
	//check if the number of maximum bonds is consistent with the maximum number of bonds given for ingredients
	if ( this->getDestination().getMolecules().getMaxConnectivity() >= reactivity.getNumMaxLinks() ) 
	{
	    //use n-1 as index, because bfm-files start counting indices at 1 (not 0)
	    molecules[n-1].setMonomerReactivity(reactivity);
// 	    std::cout << "idx" << n-1 << " reactivity: " << reactivity.isReactive() << " numMaxBonds" << reactivity.getNumMaxLinks() << std::endl;
	}else 
	{
    	std::stringstream messagestream;
      messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
		   <<"the numMaxBonds for the current monomer is exceeding the max number \n"
		   <<"of allowed connectivity for ingredients in "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());
	}
      }
      nReactiveBlocks++;
      getline((source),line);

    }
    //otherwise throw an exception
    else{

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"could not read reactivity information in reactivity definition no "<<nReactiveBlocks+1;
      throw std::runtime_error(messagestream.str());

    }
  }
}


//! Executes the routine to write \b !reactivity.
template < class IngredientsType>
void WriteReactivity<IngredientsType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"!reactivity\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers=molecules.size();
  //reactivity blocks begin with startIndex
  size_t startIndex=0;
  //counter variable
  size_t n=0;
  //reactivity to be written (updated in loop below)
  MonomerReactivity reactivity=molecules[0].getMonomerReactivity();

  //write reactivity (blockwise)
  while(n<nMonomers){
    if(molecules[n].getMonomerReactivity()!=reactivity)
    {
      strm<<startIndex+1<<"-"<<n<<":"<<reactivity<<std::endl;
      reactivity=molecules[n].getMonomerReactivity();
      startIndex=n;
    }
    n++;
  }
  //write final reactivity
  strm<<startIndex+1<<"-"<<nMonomers<<":"<<reactivity<<std::endl<<std::endl;
}
#endif
