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

#ifndef LEMONADE_FEATURE_FEATURELABEL_H
#define LEMONADE_FEATURE_FEATURELABEL_H
#include <string>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>

#include "FeatureSystemInformationTendomer.h"
#include "../updater/moves/MoveLabelBase.h"
#include "../updater/moves/MoveLabelAlongChain.h"

/**
 * @class MonomerLabel
 * @brief set the monomer labeled or not 
 */
class MonomerLabel
{
  public:
	//! Standard constructor 
	MonomerLabel():label(0){}

	//! is labeled 
	const uint32_t getLabel() const {return label;}
	//! set label on/off
	void setLabel(uint32_t label_){ label=label_;}
	
private:
     //! Private variable holding the tag. Default is NULL.
     uint32_t label;

};

/*****************************************************************/
/**
 * @class ReadLabel
 *
 * @brief Handles BFM-File-Reads \b !label
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType>
class ReadLabel: public ReadToDestination<IngredientsType>
{
public:
  ReadLabel(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadLabel(){}
  virtual void execute();
};


/*****************************************************************/
/**
 * @class WriteLabel
 *
 * @brief Handles BFM-File-Write \b !label
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteLabel:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !label into the header of the bfm-file.
  WriteLabel(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}
  virtual ~WriteLabel(){}
  virtual void writeStream(std::ostream& strm);
};
/*****************************************************************************/
/**
 * @file
 * @date   2019/05/09
 * @author Toni
 *
 * @class FeatureLabel
 * @brief This Feature handles labeled monomers.
 *
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////
class FeatureLabel : public Feature {
  
public:
	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerLabel) monomer_extensions;
	typedef LOKI_TYPELIST_1(FeatureSystemInformationTendomer) required_features_front;
	//constructor
	FeatureLabel() :latticeFilledUp(false),nArms(2),nTendomers(0), nMonomersPerChain(0)
	{labelLattice.setupLattice();}
	
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
	
	//! check move for label move 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLabelAlongChain& move) const ;
	
	//! check add move is always true 
	template<class IngredientsType, class TagType, class SpecializedMove>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBase<SpecializedMove,TagType >& move) const{return true;};

	//! check connection move is always true 
	template<class IngredientsType, class SpecializedMove>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove >& move) const{return true;};
	
	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move){};

	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLabelAlongChain& move);
	
	//!apply move for the all AddMoves
	template<class IngredientsType, class TagType, class SpecializedMove>
	void applyMove(IngredientsType& ing, const MoveAddMonomerBase<SpecializedMove, TagType >& move);
	
	//!apply move for the all AddMoves
	template<class IngredientsType, class SpecializedMove>
	void applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove >& move);
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	
	//getter for the label connected to another label. the return value starts at 1! but the ID starts at 0 
	inline uint32_t getLabelPartner(uint32_t ID) const 
	{
	  auto pair(ConnectedLabel.find(ID+1)); 
	  if ( pair != ConnectedLabel.end() ) return pair->second;
	  else return 0;
	}
	
protected:

	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;
	//!
	Lattice<uint32_t> labelLattice;
	
private:
	uint32_t nArms, nTendomers, nMonomersPerChain;

  	//! get the lattice coordinates from the id 
	std::map<uint32_t,VectorInt3> IDToCoordiantes;
	
	//! key:ID (starting from one) value: neighbor(starting from 1, if 0 = no neighbor) 
	std::map<uint32_t,uint32_t> ConnectedLabel;

};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !label
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureLabel::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("!label",new ReadLabel<IngredientsType>(fileReader.getDestination()));
}
/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !label
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureLabel::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("!label",new WriteLabel<IngredientsType>(fileWriter.getIngredients_()));
}
/******************************************************************************/
/**
 * @fn bool FeatureLabel::checkMove( const IngredientsType& ingredients, const MoveLabelAlongChain& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 * @details  it might make a difference for the speed if the order of statements is switched for different systems parameters
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType> 
bool FeatureLabel::checkMove(const IngredientsType& ingredients, const MoveLabelAlongChain& move) const
{
  if (!latticeFilledUp)
      throw std::runtime_error("*****FeatureLabel::checkMove....lattice is not populated. Run synchronize!\n");
  uint32_t ID(move.getIndex());
  //check if the monomer has a label
  if(ingredients.getMolecules()[ID].getLabel()==0) return false;
  //either the chain ends are reached or the monomer is already occupied 
  if ( labelLattice.getLatticeEntry( IDToCoordiantes.at(ID)+VectorInt3(0,0,move.getDir()) ) > 0 ) return false; 
  //if still here, then the two monomers are allowed to connect 
  return true;
}
/******************************************************************************/
/**
 * @fn void FeatureLabel ::applyMove(IngredientsType& ing, const MoveLabelAlongChain& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move  MoveLabelAlongChain.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureLabel::applyMove(IngredientsType& ing,const MoveLabelAlongChain& move)
{
  uint32_t MonID(move.getIndex()); 
  uint32_t NewID(MonID+move.getDir());
  //update lattice 
  labelLattice.moveOnLattice(IDToCoordiantes.at(MonID),IDToCoordiantes.at(NewID));
  uint32_t OtherLabel(getLabelPartner(MonID)); 
  if (OtherLabel>0){
    MonID++;NewID++;
    //update ConnectionTable
    ConnectedLabel.at(MonID)=0;
    ConnectedLabel.at(OtherLabel)=NewID;
    auto it = ConnectedLabel.lower_bound(NewID);
    if (it != ConnectedLabel.end() && it->first == NewID) {
      //key already exists
      ConnectedLabel.at(NewID)=OtherLabel;
    } else {
      ConnectedLabel.emplace_hint(it, NewID, OtherLabel);
    }
  }
}
/******************************************************************************/
/**
 * @fn void FeatureLabel ::applyMove(IngredientsType& ing, const MoveAddMonomerBase<SpecializedMove, TagType >& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType, class TagType, class SpecializedMove>
void FeatureLabel  ::applyMove(IngredientsType& ing,const MoveAddMonomerBase<SpecializedMove, TagType >& move)
{
  uint32_t MonID(move.getMonomerIndex()); 
  if ( move.getLabel() != 0 ) 
  {
    ing.modifyMolecules()[MonID].setLabel(move.getLabel());
    uint32_t tendomerID(floor(MonID/(1.0*nArms*nMonomersPerChain)));
    uint32_t IdOnTendomer(MonID%(nArms*nMonomersPerChain));
    uint32_t ArmID(floor(IdOnTendomer/(1.0*nMonomersPerChain)));
    uint32_t SegmentalID(IdOnTendomer%nMonomersPerChain);
    VectorInt3 pos(IdOnTendomer,ArmID,SegmentalID+1);
    labelLattice.setLatticeEntry(pos,move.getLabel());
  }
}
/******************************************************************************/
/**
 * @fn void FeatureLabel ::applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove >& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move for connection 
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove>
void FeatureLabel  ::applyMove(IngredientsType& ing,const MoveConnectBase<SpecializedMove >& move)
{
  uint32_t MonID(move.getIndex()); 
  if ( ing.getMolecules()[MonID].getLabel() != 0 ) 
  {
    for(size_t j =0; j < ing.getMolecules().getNumLinks(MonID);j++)
    {
      uint32_t Neighbor(ing.getMolecules().getNeighborIdx(MonID,j));
      uint32_t NeighborLabel(ing.getMolecules()[Neighbor].getLabel());
      
      if( (NeighborLabel>0) && (NeighborLabel != ing.getMolecules()[MonID].getLabel()) )
      {
// 	std::cout << "Connected Labels are : " << MonID << " " << Neighbor <<std::endl;
	ConnectedLabel[MonID+1]=Neighbor+1;
	ConnectedLabel[Neighbor+1]=MonID+1; 
      }
    }
  }
}
/******************************************************************************/
/**
 * @fn void FeatureLabel ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureLabel  ::synchronize(IngredientsType& ingredients)
{

  nTendomers = ingredients.getNumTendomers();
  nMonomersPerChain = ingredients.getNumMonomersPerChain();
  
  std::cout << "FeatureLabel::synchronizing lattice occupation...\n";
  fillLattice(ingredients);
  std::cout << "done\n";
  std::cout << "FeatureLabel::synchronizing ConnectedLabel occupation...\n";
  ConnectedLabel.clear();
  for (size_t i =0; i < ingredients.getMolecules().size();i++)
  {
    uint32_t MonLabel(ingredients.getMolecules()[i].getLabel());
    if (MonLabel>0)
    {
      for(size_t j =0; j < ingredients.getMolecules().getNumLinks(i);j++)
      {
	uint32_t Neighbor(ingredients.getMolecules().getNeighborIdx(i,j));
	uint32_t NeighborLabel(ingredients.getMolecules()[Neighbor].getLabel());
	if( (NeighborLabel>0) && (NeighborLabel != MonLabel)  &&  (Neighbor>i))
	{
// 	  std::cout <<"Add " << i << " and " << Neighbor << " to ConnectedLabel\n";
	  ConnectedLabel[i+1]=Neighbor+1;
	  ConnectedLabel[Neighbor+1]=i+1; 
	}
      }
    }
  }
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
void FeatureLabel::fillLattice(IngredientsType& ingredients)
{
  labelLattice.setupLattice(nTendomers,nArms,2+nMonomersPerChain);
  const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
  //copy the lattice occupation from the monomer coordinates
  //number of tendomers 
  for(uint32_t n=0;n<nTendomers;n++)
    //number of arms
    for(uint32_t a=0;a<nArms;a++)
      //number of monomers per arm
      for(uint32_t m=0;m<nMonomersPerChain+2;m++)
      {
	  VectorInt3 pos(n,a,m);
	  //marks the chain ends and cannot be occupied 
	  if(m==0 || m == nMonomersPerChain+1)
	    labelLattice.setLatticeEntry(pos,1);
	  else
	  {
	    uint32_t ID(m+n*nMonomersPerChain*nArms+a*nMonomersPerChain-1);
	    int32_t label;
	    if(ID <molecules.size())
	      label=molecules[ID].getLabel();
	    else 
	      label=0;
	    if( labelLattice.getLatticeEntry(pos)!=0 )
	    {
		    throw std::runtime_error("********** FeatureLabel::fillLattice: multiple lattice occupation ******************");
	    }
	    else if (label > 0)
	    {
		    // here we simply set the monomer label (plus one!) on the lattice site 
		    // the offset implies that the index zero is still used for unoccupied
		    labelLattice.setLatticeEntry(pos,label+1);
	    }
	    else if (label == 0)
	    {
		    //here we simply set the monomer label on the lattice site for the 
		    // unoccupied monomers 
		    labelLattice.setLatticeEntry(pos,label);
	    }
	    IDToCoordiantes[ID]=pos;
// 	    std::cout << ID << " " << IDToCoordiantes[ID] << " " << label <<std::endl; 
	  }
// 	  std::cout << "labelLattice["<<n<<","<<a<<","<<m<<"]="<<labelLattice.getLatticeEntry(n,a,m)<<std::endl;
      }
  latticeFilledUp=true;
}

/**
 * @brief Executes the reading routine to extract \b !label.
 *
 * @throw <std::runtime_error> label and identifier could not be read.
 **/
template < class IngredientsType >
void ReadLabel<IngredientsType>::execute()
{
    std::cout<<"reading labels of monomers ...\n";
  //some variables used during reading
  //counts the number of attribute lines in the file
  int nGroupLabels=0;
  int index(0);
  uint32_t label1, label2;
  std::string sequence1,sequence2;
  uint32_t groupLabel;
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

    //read sequence along the chain components
    getline(stream,sequence1,':');
    
    //throw exception, if extraction fails
    if(stream.fail()){
      std::stringstream messagestream;
      messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"Could not read first index in groupLabels line "<<nGroupLabels+1;
      throw std::runtime_error(messagestream.str());
    }

    //read the value of the label 
    stream>>label1;
    
    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

        std::stringstream messagestream;
      messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"Wrong definition of groupLabels\nCould not find separator \"-\" "
                   <<"in attribute definition no "<<nGroupLabels+1;
      throw std::runtime_error(messagestream.str());
    }
    //next tendomer arm:
    //read sequence along the chain components
    getline(stream,sequence2,':');

     //throw exception, if extraction fails
    if(stream.fail()){
      std::stringstream messagestream;
      messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"Could not read first index in groupLabels line "<<nGroupLabels+1;
      throw std::runtime_error(messagestream.str());
    }
    
    //read the value of the label 
    stream>>label2;
    
    //if extraction worked, save the labels
    if(!stream.fail())
    {

      for(auto i=0;i<sequence1.length();i++,index++){
	if (sequence1.compare(i,1,"1") == 0 )
	  molecules[index].setLabel(label1);
	else if (sequence1.compare(i,1,"0") == 0 )
	  molecules[index].setLabel(0);
	else {
	  std::stringstream messagestream;
	  messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"Could not read if monomer is occupied or not\n"
		   <<"with label at position"<< i <<"\n";
	  throw std::runtime_error(messagestream.str());
	}
      }
      for(auto i=0;i<sequence2.length();i++,index++){
	if (sequence2.compare(i,1,"1") == 0 )
	  molecules[index].setLabel(label2);
	else if (sequence2.compare(i,1,"0") == 0 )
	  molecules[index].setLabel(0);
	else {
	  std::stringstream messagestream;
	  messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"Could not read if monomer is occupied or not\n"
		   <<"with label at position"<< i <<"\n";
	  throw std::runtime_error(messagestream.str());
	}
      }

      nGroupLabels++;
      getline((source),line);

    }else
    {	
      //otherwise throw an exception
      std::stringstream messagestream;
      messagestream<<"ReadLabel<IngredientsType>::execute()\n"
                   <<"could not read groupLabel in groupLabel definition number "<<nGroupLabels+1;
      throw std::runtime_error(messagestream.str());
    }
  }
}


//! Executes the routine to write \b !label.
template < class IngredientsType>
void WriteLabel<IngredientsType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"!label\n";
  //counter variable
  size_t n=0;
  //
  size_t nTendomers(0);
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();
  size_t nMonomers=molecules.size();
  while(n<nMonomers && nTendomers<this->getSource().getNumTendomers())
  {
    uint32_t groupLabel;
    std::stringstream num; 
    for(uint32_t i=0; i < this->getSource().getNumMonomersPerChain();i++)
    {
      if(molecules[n].getLabel() > 0 )
      {
	strm<<"1";
	groupLabel=molecules[n].getLabel();
      }
      else 
	strm<<"0";
      n++;
    }
    strm<<":"<<groupLabel<<"-";
    for(uint32_t i=0; i < this->getSource().getNumMonomersPerChain();i++)
    {
      if(molecules[n].getLabel() > 0 )
      {
	strm<<"1";
	groupLabel=molecules[n].getLabel();
      }
      else 
	strm<<"0";
      n++;
    }
    strm<<":"<<groupLabel;
    strm<<"\n";
    nTendomers++;
  }
  strm<<std::endl;
  
}
#endif /* LEMONADE_FEATURE_FEATURELABEL_H */
