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

#ifndef FEATURE_SPRINGPOTENTIAL_TWOGROUPS_H
#define FEATURE_SPRINGPOTENTIAL_TWOGROUPS_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>

/**
 * @file
 * @brief Applying a harmonic spring potential to the center of mass positions of two groups.
 * */

/**
 * @class MonomerSpringPotentialGroupTag
 * @brief Extends monomers by an unsigned integer (uint32_t) as group tag along with getter and setter\n
 * 		  Initially the tag is set to 0, that is FeatureSpringPotentialTwoGroups::UNAFFECTED.
 * */

class MonomerSpringPotentialGroupTag{
public:
		//! constructor setting the group tag to 0 = unaffected by default
      MonomerSpringPotentialGroupTag():tagGroup(0){}
		//! getter of the group Tag
      uint32_t getMonomerGroupTag() const {return tagGroup;}
		/**
			 * @brief Setting the group tag of the monomer with \para tagGroup_.
			 * 
			 * @detail The tag must be of value 0 (=unaffected), 1 (=groupA) or 2 (=groupB).
			 *
			 * @param tagGroup_
			 */
      void setMonomerGroupTag(uint32_t tagGroup_){
			if( tagGroup_ == 0 || tagGroup_ == 1 || tagGroup_ == 2 ){
				tagGroup = tagGroup_;
			}else{
				throw std::runtime_error("FeatureSpringPotentialTwoGroups::setMonomerGroupTag not of value 0, 1 or 2\n");
			}
		}

private:
		//! Private variable holding the group tag. Default is 0.
    uint32_t tagGroup;
};


/*****************************************************************/
/**
 * @class FeatureSpringPotentialTwoGroups
 * @brief Extends vertex/monomer by an group tag (MonomerSpringPotentialGroupTag). Provides read/write functionality 
 * Implements the harmonic potential as external potential applied to the center of mass of the two groups. 
 **/
class FeatureSpringPotentialTwoGroups:public Feature
{
public:
	FeatureSpringPotentialTwoGroups(): equilibrium_length(0.0),spring_constant(0.0) {};
	virtual ~FeatureSpringPotentialTwoGroups(){};
	
	//! This Feature require Feature Boltzmann afterwards to evaluate the potential energy change
	typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

	//! This Feature requires a monomer_extensions: MonomerSpringPotentialGroupTag
	typedef LOKI_TYPELIST_1(MonomerSpringPotentialGroupTag) monomer_extensions;

	//! define an enum for the group identification
	enum SPRING_GROUP_ID{
	  UNAFFECTED=0,       
	  GROUPA=1,      
	  GROUPB=2
	};
	
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move)const;
	
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const;

	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, MoveLocalScDiag& move) const;
	
	//! getter function for the harmonic potential spring length r0 in V(r)=k/2(r-r0)^2
	double getEquilibriumLength() const{
		return equilibrium_length;
	}

	//! setter function for the harmonic potential spring length r0 in V(r)=k/2(r-r0)^2
	void setEquilibriumLength(double equilibriumLength) {
		equilibrium_length = equilibriumLength;
	}

	//! getter function for the harmonic potential spring constant k in V(r)=k/2(r-r0)^2
	double getSpringConstant() const{
		return spring_constant;
	}

	//! setter function for the harmonic potential spring constant k in V(r)=k/2(r-r0)^2
	void setSpringConstant(double springConstant){
		spring_constant = springConstant;
	}

	//! helper function to calculate the center of mass of an arbitrary monomer group
	template<class IngredientsType>
	VectorDouble3 getGroupCenterOfMass(const IngredientsType& ingredients,const std::vector<uint32_t>& group) const;

	template<class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);
	
	template<class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

private:
	//! equilibrium length r0 in harmonic potential V(r)=k/2(r-r0)^2
	double equilibrium_length;

	//! spring constant k in harmonic potential V(r)=k/2(r-r0)^2
	double spring_constant;

	//! contains the indices of the monomers of type affectedMonomerType
	std::vector<uint32_t> affectedMonomerGroup0;

	//! contains the indices of the monomers of type affectedMonomerType
	std::vector<uint32_t> affectedMonomerGroup1;

};


/*****************************************************************/
/**
 * @class ReadVirtualSpringConstant 
 *
 * @brief Handles BFM-File-Read \b #!spring_potential_constant
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
class ReadVirtualSpringConstant:public ReadToDestination<IngredientsType>
{
public:
	ReadVirtualSpringConstant(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){};

	virtual ~ReadVirtualSpringConstant(){};
	virtual void execute();
};

template<class IngredientsType>
void ReadVirtualSpringConstant<IngredientsType>::execute()
{
	std::cout<<"reading VirtualSpringConstant...";

	double springConstant = 0.0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	springConstant = atof(line.c_str());
	std::cout << "#!spring_potential_constant=" << (springConstant) << std::endl;

	ingredients.setSpringConstant(springConstant);
}

/*****************************************************************/
/**
 * @class ReadVirtualSpringLength 
 *
 * @brief Handles BFM-File-Read \b #!spring_potential_length
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadVirtualSpringLength: public ReadToDestination<IngredientsType>
{
public:
	ReadVirtualSpringLength(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadVirtualSpringLength(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadVirtualSpringLength<IngredientsType>::execute()
{
	std::cout<<"reading VirtualSpringLength...";

	double springLength = 0.0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	springLength = atof(line.c_str());
	std::cout << "#!spring_potential_length=" << (springLength) << std::endl;

	ingredients.setEquilibriumLength(springLength);
}

/*****************************************************************/
/**
 * @class ReadSpringPotentialGroups 
 *
 * @brief Handles BFM-File-Read \b !spring_potential_groups
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadSpringPotentialGroups: public ReadToDestination<IngredientsType>
{
public:
  ReadSpringPotentialGroups(IngredientsType& ingredients):ReadToDestination<IngredientsType>(ingredients){}
  virtual ~ReadSpringPotentialGroups(){}
  virtual void execute();
};

template < class IngredientsType>
void ReadSpringPotentialGroups<IngredientsType>::execute()
{
	std::cout<<"reading SpringPotentialGroups ...";
  //some variables used during reading
  //counts the number of attribute lines in the file
  int nGroupTags=0;
  int startIndex,stopIndex;
  uint32_t groupTag;
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
      messagestream<<"ReadSpringPotentialGroups<IngredientsType>::execute()\n"
                   <<"Could not read first index in groupTags line "<<nGroupTags+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character is not "-"
    if(!this->findSeparator(stream,'-')){

        std::stringstream messagestream;
      messagestream<<"ReadSpringPotentialGroups<IngredientsType>::execute()\n"
                   <<"Wrong definition of groupTags\nCould not find separator \"-\" "
                   <<"in attribute definition no "<<nGroupTags+1;
      throw std::runtime_error(messagestream.str());
    }

    //read bond identifier, throw exception if extraction fails
    stream>>stopIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
        std::stringstream messagestream;
      messagestream<<"ReadSpringPotentialGroups<IngredientsType>::execute()\n"
                   <<"Could not read second index in groupTags line "<<nGroupTags+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character is not ":"
    if(!this->findSeparator(stream,':')){

        std::stringstream messagestream;
      messagestream<<"ReadSpringPotentialGroups<IngredientsType>::execute()\n"
                   <<"Wrong definition of groupTag\nCould not find separator \":\" "
                   <<"in groupTags definition number "<<nGroupTags+1;
      throw std::runtime_error(messagestream.str());
    }
    //read the attribute tag
    stream>>groupTag;
    //if extraction worked, save the attributes
    if(!stream.fail()){

      //save attributes
      for(int n=startIndex;n<=stopIndex;n++)
      {
        //use n-1 as index, because bfm-files start counting indices at 1 (not 0)
        molecules[n-1].setMonomerGroupTag(groupTag);
      }
      nGroupTags++;
      getline((source),line);

		}else{	//otherwise throw an exception
        std::stringstream messagestream;
      messagestream<<"ReadSpringPotentialGroups<IngredientsType>::execute()\n"
                   <<"could not read groupTag in groupTag definition number "<<nGroupTags+1;
      throw std::runtime_error(messagestream.str());
    }
  }
}

/**
 * @brief perform the file reading using the FileImport class
 * 
 * @detail the following read commands are supported:
 * #!spring_potential_constant, #!spring_potential_length, !spring_potential_groups
 * 
 **/
template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!spring_potential_constant", new ReadVirtualSpringConstant<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("#!spring_potential_length", new ReadVirtualSpringLength<IngredientsType>(fileReader.getDestination()));
    fileReader.registerRead("!spring_potential_groups", new ReadSpringPotentialGroups<IngredientsType>(fileReader.getDestination()));
}

/***************************************************************************************/
/*****************************************************************/
/**
 * @class WriteVirtualSpringConstant
 *
 * @brief Handles BFM-File-Write \b #!spring_potential_constant
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteVirtualSpringConstant:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringConstant(const IngredientsType& ing):AbstractWrite<IngredientsType>(ing){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringConstant(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringConstant<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!spring_potential_constant=" << (this->getSource().getSpringConstant()) << std::endl<< std::endl;
}

/*****************************************************************/
/**
 * @class WriteVirtualSpringLength
 *
 * @brief Handles BFM-File-Write \b #!spring_potential_length
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteVirtualSpringLength:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringLength(const IngredientsType& ing):AbstractWrite<IngredientsType>(ing){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringLength(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringLength<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!spring_potential_length=" << (this->getSource().getEquilibriumLength()) << std::endl<< std::endl;
}

/*****************************************************************/
/**
 * @class WriteSpringPotentialGroups
 *
 * @brief Handles BFM-File-Write \b !spring_potential_groups
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteSpringPotentialGroups:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !attributes into the header of the bfm-file.
  WriteSpringPotentialGroups(const IngredientsType& ing)
    :AbstractWrite<IngredientsType>(ing){this->setHeaderOnly(true);}
  virtual ~WriteSpringPotentialGroups(){}
  virtual void writeStream(std::ostream& strm);
};

//! Function to write out the spring potential groups by monomer index, similar to FeatureAttributes
template < class IngredientsType>
void WriteSpringPotentialGroups<IngredientsType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)

  //write bfm command
  strm<<"!spring_potential_groups\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers = molecules.size();
  //groupTag blocks begin with startIndex
  size_t startIndex=0;
  //counter varable
  size_t n=0;
  //groupTag to be written (updated in loop below)
  uint32_t groupTag = molecules[0].getMonomerGroupTag();

  //write groupTags (blockwise)
  while(n<nMonomers){
    if(molecules[n].getMonomerGroupTag()!=groupTag)
    {
			if(groupTag != FeatureSpringPotentialTwoGroups::UNAFFECTED){
				strm<<startIndex+1<<"-"<<n<<":"<<groupTag<<std::endl;
			}
      groupTag = molecules[n].getMonomerGroupTag();
      startIndex=n;
    }
    n++;
  }
  //write final groupTags
	if(groupTag != FeatureSpringPotentialTwoGroups::UNAFFECTED){
  	strm<<startIndex+1<<"-"<<nMonomers<<":"<<groupTag<<std::endl;
	}
	strm<<std::endl;
}

//! perform the file writing using the AnalyzerWriteBfmFile<IngredientsType>
template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
	fileWriter.registerWrite("#!spring_potential_constant",new WriteVirtualSpringConstant<IngredientsType>(fileWriter.getIngredients_()));
	fileWriter.registerWrite("#!spring_potential_length",new WriteVirtualSpringLength<IngredientsType>(fileWriter.getIngredients_()));
	fileWriter.registerWrite("!spring_potential_groups",new WriteSpringPotentialGroups<IngredientsType>(fileWriter.getIngredients_()));
}

/***************************************************************************************/
/************* private member functions of FeatureSpringPotentialTwoGroups *************/

/**
 * This function performs the Monte-Carlo check for the basic move.
 * This is called for unknown move types.
 * It always return true.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] the basic move: MoveBase
 */
template<class IngredientsType>
bool FeatureSpringPotentialTwoGroups::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}

/**
 * This function performs the Monte-Carlo check for the MoveLocalSc.
 * The center of mass differences induced by the moveand the resulting potential difference is calculated.
 * It passes the corresponding move probability to FeatureBoltzmann. 
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move the standard simple cubic lattice move: MoveLocalSc
 */
template<class IngredientsType>
bool FeatureSpringPotentialTwoGroups::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
	//Index of moved Monomer is monoIndex
	uint32_t monoIndex=move.getIndex();

	//if the moved monomer has the same attribute as the affectedMonomerType,
	//get the z position of the group after a hypothetical move
	//otherwise return true right away
	int32_t moveGroupTag=ingredients.getMolecules()[move.getIndex()].getMonomerGroupTag();

	//this is the COM of the group where the monomer which shall be moves belongs to 
	VectorDouble3 COM_position_old;
	//this is the COM of the group where the monomer does not belong to 
	VectorDouble3 COM_position_not_moved;
	//number of monomers which belong to the group of the potentially moved monomer
	double size_moved_group = 0.0;
	
	VectorDouble3 projected_move = move.getDir();

	if(moveGroupTag == GROUPA)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		size_moved_group=affectedMonomerGroup0.size();
	}
	else if(moveGroupTag == GROUPB)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		size_moved_group=affectedMonomerGroup1.size();
	}
	else 
	  return true;

	//the rest happens only if we have not returned yet, i.e. if distance has been calculated
	
	//distance between the two COM if the move is not applied
	double rel_length_old( (COM_position_old-COM_position_not_moved).getLength() );
	//distance between the two COM if the move would be applied
	double rel_length_moved( (COM_position_old-COM_position_not_moved+projected_move/size_moved_group).getLength() );

	/* calculate the potential difference
	 * V=k/2*(|R_COM|-R_0)^2
	 * R_COM=R1-R2
	 * dV=V(R_COM(unmoved))-V(R_COM(moved))
	 * the simplified equation below assumes a step length of 1 !!!
	 */
	
// 	double dV = spring_constant*(COM_position_old*projected_move/size_moved_group);
// 	dV += 0.5*spring_constant/(size_moved_group*size_moved_group);
// 	dV -= spring_constant*(COM_position_not_moved*projected_move/size_moved_group);
// 	dV += spring_constant*equilibrium_length*(rel_length_old-rel_length_moved);
	//seems to be shorter...
	double 	dV  = 0.5/(size_moved_group*size_moved_group);
		dV += equilibrium_length*(rel_length_old-rel_length_moved);
		dV += projected_move/size_moved_group*(COM_position_old-COM_position_not_moved);
		dV *= spring_constant;
	
	//calculate the transition probability
	//Metropolis: zeta = exp (-dV)
	double prob=exp(-dV);

	//std::cout << "prob: " <<  prob << std::endl;
	move.multiplyProbability(prob);

	return true;

}
/**
 * This function performs the Monte-Carlo check for the MoveLocalScDiag.
 * The center of mass differences induced by the moveand the resulting potential difference is calculated.
 * It passes the corresponding move probability to FeatureBoltzmann. 
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move the simple cubic lattice move: MoveLocalScDiag
 */
template<class IngredientsType>
bool FeatureSpringPotentialTwoGroups::checkMove(const IngredientsType& ingredients, MoveLocalScDiag& move) const
{
  	//Index of moved Monomer is monoIndex
	uint32_t monoIndex=move.getIndex();

	//if the moved monomer has the same attribute as the affectedMonomerType,
	//get the z position of the group after a hypothetical move
	//otherwise return true right away
	int32_t moveGroupTag=ingredients.getMolecules()[move.getIndex()].getMonomerGroupTag();

	//this is the COM of the group where the monomer which shall be moves belongs to 
	VectorDouble3 COM_position_old;
	//this is the COM of the group where the monomer does not belong to 
	VectorDouble3 COM_position_not_moved;
	//number of monomers which belong to the group of the potentially moved monomer
	double size_moved_group = 0.0;
	
	VectorDouble3 projected_move = move.getDir();

	if(moveGroupTag == GROUPA)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		size_moved_group=affectedMonomerGroup0.size();
	}
	else if(moveGroupTag == GROUPB)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		size_moved_group=affectedMonomerGroup1.size();
	}
	else 
	  return true;

	//the rest happens only if we have not returned yet, i.e. if distance has been calculated
	
	//distance between the two COM if the move is not applied
	double rel_length_old( (COM_position_old-COM_position_not_moved).getLength() );
	//distance between the two COM if the move would be applied
	double rel_length_moved( (COM_position_old-COM_position_not_moved+projected_move/size_moved_group).getLength() );

	/* calculate the potential difference
	 * V=k/2*(|R_COM|-R_0)^2
	 * R_COM=R1-R2
	 * dV=V(R_COM(unmoved))-V(R_COM(moved))
	 */
	double 	dV  = 0.5*projected_move.getLength()*projected_move.getLength()/(size_moved_group*size_moved_group);
		dV += equilibrium_length*(rel_length_old-rel_length_moved);
		dV += projected_move/size_moved_group*(COM_position_old-COM_position_not_moved);
		dV *= spring_constant;
	
	//calculate the transition probability
	//Metropolis: zeta = exp (-dV)
	double prob=exp(-dV);

	//std::cout << "prob: " <<  prob << std::endl;
	move.multiplyProbability(prob);

	return true;
  
}
/**
 * Performs the synchronize for the utilities of the feature:
 *   Clear the monomer groups and refill them by reading the monomer Groups Tag.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] the standard simple cubic lattice move: MoveLocalSc
 */
template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::synchronize(IngredientsType& ingredients)
{
	// delete the old content
	affectedMonomerGroup0.clear();
	affectedMonomerGroup1.clear();

	//sort the monomers into groups
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		if(ingredients.getMolecules()[n].getMonomerGroupTag()==GROUPA)
		{
			affectedMonomerGroup0.push_back(n);
		}

		if(ingredients.getMolecules()[n].getMonomerGroupTag()==GROUPB)
		{
			affectedMonomerGroup1.push_back(n);
		}
	}

	std::cout<<"FeatureSpringPotentialTwoGroups::synchronize()...affected group size 1 ="<<affectedMonomerGroup0.size()<<
			"...affected group size 2 ="<<affectedMonomerGroup1.size()<<std::endl;

}

template<class IngredientsType>
VectorDouble3 FeatureSpringPotentialTwoGroups::getGroupCenterOfMass(const IngredientsType& ingredients, const std::vector<uint32_t>& group) const
{
	//the default for index is 0, the default for direction is 0,0,0. the index 0 points to a particle,
	//but since one has to explicitly specify index if one changes direction to anything other than 0,0,0
	//there is no danger in using the function without explicit arguments

	int32_t sumX=0;
	int32_t sumY=0;
	int32_t sumZ=0;
	for(size_t n=0;n<group.size();n++)
	{
			sumX+=ingredients.getMolecules()[group[n]].getX();
			sumY+=ingredients.getMolecules()[group[n]].getY();
			sumZ+=ingredients.getMolecules()[group[n]].getZ();

	}
	return VectorDouble3(sumX,sumY,sumZ)/(group.size());
}


#endif /*FEATURE_SPRINGPOTENTIAL_TWOGROUPS_H_H*/
