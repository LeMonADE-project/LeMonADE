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

#ifndef LEMONADE_FEATURE_FEATUREREACTIVEBONDS_H
#define LEMONADE_FEATURE_FEATUREREACTIVEBONDS_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/feature/FeatureBreak.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveBreakBase.h>
#include <LeMonADE/updater/moves/MoveBreak.h>
#include <LeMonADE/updater/moves/MoveBreakReactive.h>
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
	 bool isReactive() const {return reactivity;}
	
	//! Getting the number of maximum possible bonds for the monomer.
	 uint32_t getNumMaxLinks() const {return numMaxLinks;};

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
	 bool operator == (const MonomerReactivity &react) const 
	{
	  if ( react.getNumMaxLinks() != numMaxLinks  ) return false;
	  if ( react.isReactive()   != reactivity ) return false;
	  return true;
	}
	 bool operator!= (const MonomerReactivity &react) const 
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
 * @class FeatureReactiveBonds
 * @brief This Feature add new bonds between reactive monomers.
 *
 * @details Works only in combination with an excluded volume feature. 
 * Considers only reactive monomers which are able to form a new bond.
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureReactiveBonds : public Feature {
    typedef std::pair < uint32_t, uint32_t > BondPair;
    typedef uint32_t edge_type;
public:
    //! This Feature requires a monomer_extensions.
    typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;
      typedef LOKI_TYPELIST_2(FeatureBreak, FeatureConnectionSc) required_features_back;
      
      FeatureReactiveBonds():nReactedBonds(0),nReactiveSites(0){};
      
    //! Export the relevant functionality for reading bfm-files to the responsible reader object
    template<class IngredientsType>
    void exportRead(FileImport<IngredientsType>& fileReader); 

    //! Export the relevant functionality for writing bfm-files to the responsible writer object
    template<class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
    
      //! check base move - always true 
    template<class IngredientsType> 
    bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const{return true;};
      
    //! check base connection move  
    template<class IngredientsType, class SpecializedMove> 
    bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const;

      //! check base break move 
    template<class IngredientsType, class SpecializedMove >
    bool checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const ;
      
      //! check movebreakreactive
    template<class IngredientsType >
    bool checkMove(const IngredientsType& ingredients, const MoveBreakReactive& move) const {return true;};
      
      //!apply move for the connection moves 
    template<class IngredientsType >
    void applyMove(IngredientsType& ing, const MoveBase& move){};	
      
      //!apply move for the connection moves 
    template<class IngredientsType, class SpecializedMove>
    void applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove>& move);	
      
      //!apply move for the breaking moves 
    template<class IngredientsType, class SpecializedMove>
    void applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move);	
      
    //! Synchronize, count the number of reactive sites and reacted bonds 
    template<class IngredientsType>
    void synchronize(IngredientsType& ingredients);
    
    //!returns the number of bond made of two reactive monomers 
    uint32_t getNReactedBonds() const { return BondedReactiveMonomers.size(); };

    //!returns the total number of reactive monomers capable to form a bond
    uint32_t getNReactiveSites() const {return nReactiveSites;};

    //!returns the bond table of the reacted reactive monomers.
    const std::map<BondPair,edge_type>& getBondedMonomers()const {return BondedReactiveMonomers;};

    //! returns false if the bond does not exist
    bool checkReactiveBondExists(uint32_t Mon1, uint32_t Mon2) const 
    {
      BondPair edge_key(std::min(Mon1,Mon2),std::max(Mon1,Mon2));
      return (BondedReactiveMonomers.find(edge_key) != BondedReactiveMonomers.end()); 
    }

    //! returns the map of unreacted monomers 
    const std::map<uint32_t,uint32_t>& getUnreactiveMonomers() const {return UnbondedReactiveMonomers;};

    //!returns the number of reactive monomers capable to form a bond 
    uint32_t getNUnreactedMonomers()const{return UnbondedReactiveMonomers.size();}

    //!returns true if the monomer can be connected to another monomer
    bool checkCapableFormingBonds(uint32_t MonID )const { return (UnbondedReactiveMonomers.find(MonID) != UnbondedReactiveMonomers.end()); }

    //!get extent of reaction 
    double getConversion() const {return (double(nReactedBonds*2))/(double(nReactiveSites));}
    
private: 
    //!number of bonds from reactive monomers and the total number of reactive monomers 
    uint32_t nReactedBonds, nReactiveSites;

    //! holds all bonded reactive monomers during simulation
    std::map<BondPair, edge_type> BondedReactiveMonomers;

    //! holds the ids of all reactive monomers which can have another bonds
    std::map<uint32_t,uint32_t> UnbondedReactiveMonomers; 

    //! adds the monomer to the map of unreacted monomers
    void addReactiveMonomer(uint32_t MonID, uint32_t value=0){UnbondedReactiveMonomers[MonID]=value;} // adds a new entry if the 

    //! erase the monomer ID from the map of unreacted monomers 
    void eraseReactiveMonomer(uint32_t MonID, uint32_t value=0){UnbondedReactiveMonomers.erase(UnbondedReactiveMonomers.find(MonID));}

    //!add a bond to the container BondedReactiveMonomers
    void addBondedPair(uint32_t Monomer1, uint32_t Monomer2, uint32_t attribute=0){
      BondPair edge_key(std::min(Monomer1,Monomer2),std::max(Monomer1,Monomer2));
      BondedReactiveMonomers[edge_key]=attribute;
    }

    //!erase a bond from the container BondedReactiveMonomers
    void eraseBondedPair(uint32_t Monomer1, uint32_t Monomer2, uint32_t attribute=0){
      BondPair edge_key(std::min(Monomer1,Monomer2),std::max(Monomer1,Monomer2));
      typename std::map<BondPair, edge_type>::iterator it;
      it=BondedReactiveMonomers.find(edge_key);
      BondedReactiveMonomers.erase(it);
    }
};
///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !reactivity
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureReactiveBonds::exportRead(FileImport< IngredientsType >& fileReader) {
    fileReader.registerRead("!reactivity",new ReadReactivity<IngredientsType>(fileReader.getDestination()));
}
/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !reactivity
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureReactiveBonds::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const{
    fileWriter.registerWrite("!reactivity",new WriteReactivity<IngredientsType>(fileWriter.getIngredients_()));
}
/******************************************************************************/
/**
 * @fn bool FeatureReactiveBonds::checkMove( const IngredientsType& ingredients, const MoveConnectSc& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 * @details  it might make a difference for the speed if the order of statements is switched for different systems parameters
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureReactiveBonds ::checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const {
    uint32_t ID(move.getIndex());
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    //check for maximum number of bonds for the first monomer
    if ( molecules.getNumLinks(ID) >=  molecules[ID].getNumMaxLinks()) return false;
    uint32_t Neighbor(move.getPartner());
    //check if neighbor is reactive 
    if ( !molecules[Neighbor].isReactive() ) return false;
    //check for maximum number of bonds for the second monomer
    if ( molecules.getNumLinks(Neighbor) >= molecules[Neighbor].getNumMaxLinks() ) return false;
    //if still here, then the two monomers are allowed to connect 
    return true;
}
/******************************************************************************/
/**
 * @fn bool FeatureReactiveBonds ::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 * @details  it might make a difference for the speed if the order of statements is switched for different systems parameters
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureReactiveBonds ::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const {
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    //check if neighbor is reactive 
    uint32_t ID(move.getIndex());
    if ( !molecules[ID].isReactive() ) return false;
    //check if neighbor is reactive 
    uint32_t Neighbor(move.getPartner());
    if ( !molecules[Neighbor].isReactive() ) return false;
    //if still here, then the two monomers are allowed to connect 
    return true;
}
/******************************************************************************/
/**
 * @fn void FeatureReactiveBonds ::applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove>& move)
 * @brief Updates the bond table in the feature 
 * @details  
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move of base type MoveConnectBase
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove>
void FeatureReactiveBonds ::applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove>& move) {
    nReactedBonds++;
    auto MonID(move.getIndex()); 
    auto Partner(move.getPartner());
    addBondedPair(MonID,Partner,ing.getMolecules().getAge());
    if (ing.getMolecules()[MonID].getNumMaxLinks()==ing.getMolecules().getNumLinks(MonID) )
      eraseReactiveMonomer(MonID);
    if (ing.getMolecules()[Partner].getNumMaxLinks()==ing.getMolecules().getNumLinks(Partner) )
      eraseReactiveMonomer(Partner);
}
/******************************************************************************/
/**
 * @fn void FeatureReactiveBonds ::applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move)
 * @brief Updates the bond table in the feature 
 * @details  
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move of base type MoveBreakBase
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove>
void FeatureReactiveBonds ::applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move) {
    nReactedBonds--;
    auto MonID(move.getIndex()); 
    auto Partner(move.getPartner());
    eraseBondedPair(MonID,Partner);
    addReactiveMonomer(MonID,ing.getMolecules().getAge());
    addReactiveMonomer(Partner,ing.getMolecules().getAge());
  
}
/******************************************************************************/
/**
 * @fn void FeatureReactiveBonds ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the number number of reactive sites and reacted sites.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureReactiveBonds::synchronize(IngredientsType& ingredients){
    std::cout << "FeatureReactiveBonds::synchronizing ...\n";
    nReactedBonds=0;
    nReactiveSites=0;
    for(size_t i = 0 ; i < ingredients.getMolecules().size(); i++ ){
        if ( ingredients.getMolecules()[i].isReactive() ){
            addReactiveMonomer(i,ingredients.getMolecules().getAge());
			uint32_t NLinks(ingredients.getMolecules().getNumLinks(i));
			uint32_t nIrreversibleBonds=0;
			for (uint32_t n = 0 ; n < NLinks ;n++){
				uint32_t neighbor(ingredients.getMolecules().getNeighborIdx(i,n));
				if( ingredients.getMolecules()[neighbor].isReactive() )
					nReactedBonds++;
				else
				nIrreversibleBonds++;
			}
			if(ingredients.getMolecules()[i].getNumMaxLinks()-ingredients.getMolecules().getNumLinks(i) != 0)
				addReactiveMonomer(i,ingredients.getMolecules().getAge());
			nReactiveSites+=(ingredients.getMolecules()[i].getNumMaxLinks()-nIrreversibleBonds);
      }
    }
    nReactedBonds/=2; // they are counted twice: each monomer counts the bond, but essentially there is only one bond.
	std::cout << "Number of reactive bonds: " << nReactedBonds <<"\n"
				<< "Number of reactive sites: " << nReactiveSites<<"\n"
				<< "Extent of reaction      : " << getConversion()
				<<std::endl;
	if (nReactiveSites != 0 ){
    	BondedReactiveMonomers.clear();
		auto edges=ingredients.getMolecules().getEdges();
		for(auto it=edges.begin();it!=edges.end();++it){
			auto MonID1(it->first.first);
			auto MonID2(it->first.second);
			if(ingredients.getMolecules()[MonID1].isReactive() && ingredients.getMolecules()[MonID2].isReactive())
				addBondedPair(MonID1,MonID2,ingredients.getMolecules().getAge());
		}
      }
    std::cout << "done\n";
}
///////////////////////////////////////////////////////////////////////////////
///////////////////// member definitions for output ///////////////////////////
/******************************************************************************/
/**
 * @brief Executes the reading routine to extract \b !reactivity.
 *
 * @throw <std::runtime_error> reactivity and identifier could not be read.
 **/
/******************************************************************************/
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
			for(int n=startIndex;n<=stopIndex;n++){
				//check if the number of maximum bonds is consistent with the maximum number of bonds given for ingredients
				if ( this->getDestination().getMolecules().getMaxConnectivity() >= reactivity.getNumMaxLinks() ) {
					//use n-1 as index, because bfm-files start counting indices at 1 (not 0)
					molecules[n-1].setMonomerReactivity(reactivity);
				}else{
					std::stringstream messagestream;
					messagestream<<"ReadReactivity<IngredientsType>::execute()\n"
						<<"the numMaxBonds for the current monomer is exceeding the max number \n"
						<<"of allowed connectivity for ingredients in "<<nReactiveBlocks+1;
					throw std::runtime_error(messagestream.str());
				}
			}
			nReactiveBlocks++;
			getline((source),line);
		}else{//otherwise throw an exception
			std::stringstream messagestream;
			messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
				<<"could not read reactivity information in reactivity definition no "<<nReactiveBlocks+1;
			throw std::runtime_error(messagestream.str());
		}
	}
}
//! Executes the routine to write \b !reactivity.
template < class IngredientsType>
void WriteReactivity<IngredientsType>::writeStream(std::ostream& strm){
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
		if(molecules[n].getMonomerReactivity()!=reactivity){
			if( reactivity.isReactive() == true )
				strm<<startIndex+1<<"-"<<n<<":"<<reactivity<<std::endl;
			reactivity=molecules[n].getMonomerReactivity();
			startIndex=n;
		}
		n++;
	}
	//write final reactivity
	if( reactivity.isReactive() == true )
		strm<<startIndex+1<<"-"<<nMonomers<<":"<<reactivity<<std::endl;
	strm<<std::endl;
}
#endif
