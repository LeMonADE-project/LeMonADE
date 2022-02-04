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


#ifndef FEATURE_NN_INTERACTION_H
#define FEATURE_NN_INTERACTION_H

/**
 * @file
 * @date 2016/06/18, 2022/02/01
 * @author Hauke Rabbel, Toni Mueller
 * @brief Definition and implementation of class template FeatureNNInteractionSc
**/
#include <iostream>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/feature/FeatureNNInteractionReadWrite.h>

/**
 * @file
 * @brief Adding an interaction tag to each monomer.
 * */

/**
 * @class MonomerInteractionTag
 * @brief Extends monomers by an signed integer (int32_t) as tag along with getter and setter\n
 * 		  Initially the tag is set to NULL.
 * */
class MonomerInteractionTag
{
public:
	//! Standard constructor- initially the tag is set to NULL.
	MonomerInteractionTag():tag(0){}
	//! Getting the tag of the monomer.
	uint8_t getInteractionTag() const {return tag;}
	/**
	 * @brief Setting the tag of the monomer with \para attr.
	 *
	 * @param attr
	 */
	void setInteractionTag(uint8_t interaction){ tag=interaction;}
private:
	//! Private variable holding the tag. Default is NULL.
	uint8_t tag;
};

/*****************************************************************/
/**
 * @class ReadInteractionTags
 *
 * @brief Handles BFM-File-Reads \b !interactionTag
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType>
class ReadInteractionTags: public ReadToDestination<IngredientsType>
{
public:
  ReadInteractionTags(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadInteractionTags(){}
  virtual void execute();
};
/*****************************************************************/
/**
 * @class WriteInteractionTags
 *
 * @brief Handles BFM-File-Write \b !interactionTag
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteInteractionTags:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !interactionTag into the header of the bfm-file.
	WriteInteractionTags(const IngredientsType& i)
		:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
	virtual ~WriteInteractionTags(){}
	virtual void writeStream(std::ostream& strm);
};

/**
 * @brief Executes the reading routine to extract \b !interactionTag.
 *
 * @details interactionTag must be uint32_t for reading. For uint8_t the input 
 * stream  would be interpreted as char and not as number! 
 * 
 * @throw <std::runtime_error> interactionTag and identifier could not be read.
 **/
template < class IngredientsType>
void ReadInteractionTags<IngredientsType>::execute()
{
	//some variables used during reading
	//counts the number of interactionTag lines in the file
	int nAttributes=0;
	int startIndex,stopIndex;
	uint32_t interactionTag;
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
			messagestream<<"ReadInteractionTags<IngredientsType>::execute()\n"
				<<"Could not read first index in interactionTag line "<<nAttributes+1;
			throw std::runtime_error(messagestream.str());
		}
		//throw exception, if next character isnt "-"
		if(!this->findSeparator(stream,'-')){
			std::stringstream messagestream;
			messagestream<<"ReadInteractionTags<IngredientsType>::execute()\n"
				<<"Wrong definition of interactionTag\nCould not find separator \"-\" "
				<<"in interactionTag definition no "<<nAttributes+1;
			throw std::runtime_error(messagestream.str());
		}
		//read bond identifier, throw exception if extraction fails
		stream>>stopIndex;

		//throw exception, if extraction fails
		if(stream.fail()){
			std::stringstream messagestream;
		messagestream<<"ReadInteractionTags<IngredientsType>::execute()\n"
			<<"Could not read second index in interactionTag line "<<nAttributes+1;
		throw std::runtime_error(messagestream.str());
		}

		//throw exception, if next character isnt ":"
		if(!this->findSeparator(stream,':')){
			std::stringstream messagestream;
			messagestream<<"ReadInteractionTags<IngredientsType>::execute()\n"
				<<"Wrong definition of interactionTag\nCould not find separator \":\" "
				<<"in interactionTag definition no "<<nAttributes+1;
			throw std::runtime_error(messagestream.str());
		}
		//read the interactionTag tag
		stream>>interactionTag;

		//if extraction worked, save the interactionTag
		if(!stream.fail()){
			//save interactionTag
			for(int n=startIndex;n<=stopIndex;n++)
			{
				//use n-1 as index, because bfm-files start counting indices at 1 (not 0)
				molecules[n-1].setInteractionTag(static_cast<uint8_t>(interactionTag));
			}
			nAttributes++;
			getline((source),line);
		}
		//otherwise throw an exception
		else{
			std::stringstream messagestream;
			messagestream<<"ReadInteractionTags<IngredientsType>::execute()\n"
				<<"could not read interactionTag in interactionTag definition no "<<nAttributes+1;
			throw std::runtime_error(messagestream.str());
		}
	}
}
//! Executes the routine to write \b !interactionTag.
template < class IngredientsType>
void WriteInteractionTags<IngredientsType>::writeStream(std::ostream& strm)
{
	//for all output the indices are increased by one, because the file-format
	//starts counting indices at 1 (not 0)
	//write bfm command
	strm<<"!interactionTag\n";
	//get reference to monomers
	const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();
	size_t nMonomers=molecules.size();
	//interactionTag blocks begin with startIndex
	size_t startIndex=0;
	//counter varable
	size_t n=0;
	//interactionTag to be written (updated in loop below)
	
	// uint32_t interactionTag=static_cast<uint32_t>(molecules[0].getInteractionTag());
	uint32_t interactionTag=molecules[0].getInteractionTag();
	//write interactionTag (blockwise)
	while(n<nMonomers){
		// if(static_cast<uint32_t>(molecules[n].getInteractionTag())!=interactionTag)
		if( molecules[n].getInteractionTag()!=interactionTag )
		{
			strm<<startIndex+1<<"-"<<n<<":"<<interactionTag<<std::endl;
			// interactionTag=static_cast<uint32_t>(molecules[n].getInteractionTag());
			interactionTag=molecules[n].getInteractionTag();
			startIndex=n;
		}
		n++;
	}
	//write final interactionTag
	strm<<startIndex+1<<"-"<<nMonomers<<":"<<interactionTag<<std::endl<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
/////////////IMPLEMENT FeatureNNInteractionSc /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/**
 * @class FeatureNNInteractionSc
 * @brief Provides interaction of monomers on distances d<=sqrt(6) for standard BFM
 *
 * @tparam FeatureLatticeType Underlying lattice feature, e.g. FeatureLattice or
 * FeatureLatticePowerOfTwo (template template parameter)
 *
 * @details
 * The interaction energy can be set for pairs of monomer-types A,B, where
 * the type can be any integer between 1 and 255.
 * The feature automatically adds FeatureExcludedVolumeSc<FeatureLatticeType<uint8_t>
 * to the system. Given an energy E in kT between  two types, the interaction potential
 * as a function of the distance d is:
 * - inf d<2 (implicitly through excluded volume)
 * - 4*E d=2
 * - 2*E d=sqrt(5)
 * - 1*E d=sqrt(6)
 * - 0   d>sqrt(6)
 * .
 * Usage: In the feature list defining Ingredients use this feature as
 * FeatureNNInteractionSc<FeatureLattice> (arbitrary lattices), or as
 * FeatureNNInteractionSc<FeatureLatticePowerOfTwo> (2**n lattices)
 * The feature adds the bfm-file command !nn_interaction A B E
 * for monomers of types A B with interaction energy of E in kT.
**/

class FeatureNNInteractionSc:public Feature
{
public:
	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerInteractionTag) monomer_extensions;

private:
	//! Type for the underlying lattice, used as template parameter for FeatureLatticeType<...>
	typedef uint8_t lattice_value_type;

	//! Interaction energies between monomer types. Max. type=255 given by max(uint8_t)=255
	double interactionTable[256][256];

	//! Lookup table for exp(-interactionTable[a][b])
	double probabilityLookup[256][256];

	//! Returns this feature's factor for the acceptance probability for the given Monte Carlo move
	template<class IngredientsType>
	double calculateAcceptanceProbability(const IngredientsType& ingredients,
						const MoveLocalSc& move) const;

	//! Occupies the lattice with the interactionTag tags of all monomers
	template<class IngredientsType>
	void fillLattice(IngredientsType& ingredients);

	//! Access to array probabilityLookup with extra checks in Debug mode
	double getProbabilityFactor(int32_t typeA,int32_t typeB) const;

	//! Lattice storing the species Tag on the lattice 
	Lattice<uint8_t> interactionLattice;

protected:

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;

public:
	
	FeatureNNInteractionSc();
	~FeatureNNInteractionSc(){};
	
	//This feature adds interaction energies, so it requires FeatureBoltzmann
	typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

	//! check for all Monte Carlo moves without special check functions (always true)
	template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;

	//! check for standard sc-BFM local move
	template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;

	//! check move for bcc-BFM local move. always throws std::runtime_error
	template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveLocalBcc& move) const;

	//! apply function for all Monte Carlo moves without special apply functions (does nothing)
	template<class IngredientsType>
    void applyMove(const IngredientsType& ing, const MoveBase& move){}
	
	//! apply function for sc-BFM local move
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move );

	//! apply function for bcc-BFM local move (always throws std::runtime_error)
	template<class IngredientsType>
    void applyMove(const IngredientsType& ing, const MoveLocalBcc& move);

	//! apply function for adding a monomer in sc-BFM
	template<class IngredientsType>
    void applyMove(IngredientsType& ing, const MoveAddMonomerSc<int32_t>& move);

	//note: apply function for sc-BFM local move is not necessary, because
	//job of moving lattice entries is done by the underlying FeatureLatticeType

	//! guarantees that the lattice is properly occupied with monomer interactionTag
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

	//!adds interaction energy between two types of monomers
	void setNNInteraction(int32_t typeA,int32_t typeB,double energy);

	//!returns the interaction energy between two types of monomers
	double getNNInteraction(int32_t typeA,int32_t typeB) const;

	//!export bfm-file read command !nn_interaction
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);

	//!export bfm-file write command !nn_interaction
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

	//! Get the lattice value at a certain point for the interaction species 
	uint32_t getInteractionLatticeEntry(const VectorInt3& pos) const { return interactionLattice.getLatticeEntry(pos);} ;

	//! Get the lattice value at a certain point for the interaction species 
	uint32_t getInteractionLatticeEntry(const int x, const int y, const int z) const { return interactionLattice.getLatticeEntry(x,y,z);};

};

//////////////////  IMPLEMENTATION OF MEMBERS //////////////////////////////////
/**
 * @brief Constructor
 **/
FeatureNNInteractionSc::FeatureNNInteractionSc():
latticeFilledUp(false)
{
	interactionLattice.setupLattice();
	//initialize the energy and probability lookups with default values
	for(size_t n=0;n<256;n++)
    {
      	for(size_t m=0;m<256;m++)
		{	
			interactionTable[m][n]=0.0;
			probabilityLookup[m][n]=1.0;
        }
    }
}

/**
 * @details Returns true for all moves other than the ones that have specialized versions of this function.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move other than moves with specialized functions.
 * @return true (always)
 **/
template<class IngredientsType>
bool FeatureNNInteractionSc::checkMove(const IngredientsType& ingredients,
							 const MoveBase& move) const
{
    return true;
}

/**
 * @details calculates the factor for the acceptance probability of the move
 * arising from the contact interactions and adds it to the move.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalSc
 * @return true (always)
 **/
template<class IngredientsType>
bool FeatureNNInteractionSc::checkMove(const IngredientsType& ingredients,
							 MoveLocalSc& move) const
{
	//add the probability factor coming from this feature, then return true,
	//because the total probability is evaluated by FeatureBoltzmann at the end
	double prob=calculateAcceptanceProbability(ingredients,move);
	move.multiplyProbability(prob);
	return true;
}

/**
 * @details Because moves of type MoveLocalBcc must not be used with this
 * feature, this function always throws an exception when called. The function
 * is only implemented for savety purposes.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalBcc
 * @throw std::runtime_error
 * @return false always throws exception before returning
 **/
template<class IngredientsType>
bool FeatureNNInteractionSc::checkMove(const IngredientsType& ingredients,
							 const MoveLocalBcc& move) const
{
	//throw exception in case someone accidentaly uses a bcc-BFM move with this feature
	std::stringstream errormessage;
	errormessage<<"FeatureNNInteractionSc::checkMove(...):\n";
	errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
	throw std::runtime_error(errormessage.str());
	return false;
}

/**
 * @details When a new monomer is added to the system, the lattice sites occupied
 * by this monomer must be filled with the interactionTag tag. This function takes care
 * of this.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveAddMonomerSc
 * @throw std::runtime_error if interactionTag tag is not in range [1,255]
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::applyMove(IngredientsType& ing,
							 const MoveAddMonomerSc<int32_t>& move)
{
    //get the position and interactionTag tag of the monomer to be inserted
    VectorInt3 pos=move.getPosition();
    VectorInt3 dx(1,0,0);
    VectorInt3 dy(0,1,0);
    VectorInt3 dz(0,0,1);
    uint8_t type(move.getInteractionTag());

    //the feature is based on a uint8_t lattice, thus the max type must not
    //exceed the max value of uint8_t (255)
    if(int32_t(type) != move.getInteractionTag() )
	{
		std::stringstream errormessage;
		errormessage<<"FeatureNNInteractionSc::applyMove(MoveAddMonomerSc)\n";
		errormessage<<"Trying to add monomer with type "<<int32_t(type)<<">maxType=255\n";
		throw std::runtime_error(errormessage.str());
	}
	else if (int32_t(type)!=0){
    //update lattice
    interactionLattice.setLatticeEntry(pos,type);
    interactionLattice.setLatticeEntry(pos+dx,type);
    interactionLattice.setLatticeEntry(pos+dy,type);
    interactionLattice.setLatticeEntry(pos+dx+dy,type);
    interactionLattice.setLatticeEntry(pos+dz,type);
    interactionLattice.setLatticeEntry(pos+dz+dx,type);
    interactionLattice.setLatticeEntry(pos+dz+dy,type);
    interactionLattice.setLatticeEntry(pos+dz+dx+dy,type);
	}
}
/**
 * @details When a new monomer is added to the system, the lattice sites occupied
 * by this monomer must be filled with the interactionTag tag. This function takes care
 * of this.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalSc
 * @throw std::runtime_error if interactionTag tag is not in range [1,255]
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::applyMove(IngredientsType& ing,
							const MoveLocalSc& move )
{
	//get old position and direction of the move
	VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
	VectorInt3 direction=move.getDir();

	/*get two directions perpendicular to vector directon of the move*/
	VectorInt3 perp1,perp2;
  	/* first perpendicular direction is either (0 1 0) or (1 0 0)*/
	int32_t x1=((direction.getX()==0) ? 1 : 0);
	int32_t y1=((direction.getX()!=0) ? 1 : 0);
	perp1.setX(x1);
	perp1.setY(y1);
	perp1.setZ(0);

	/* second perpendicular direction is either (0 0 1) or (0 1 0)*/
	int32_t y2=((direction.getZ()==0) ? 0 : 1);
	int32_t z2=((direction.getZ()!=0) ? 0 : 1);
	perp2.setX(0);
	perp2.setY(y2);
	perp2.setZ(z2);
	if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) oldPos-=direction;
	direction*=2;
    VectorInt3 oldPlusDir=oldPos+direction;
	//change lattice occupation accordingly
    interactionLattice.moveOnLattice(oldPos,oldPlusDir);
    interactionLattice.moveOnLattice(oldPos+perp1,oldPlusDir+perp1);
    interactionLattice.moveOnLattice(oldPos+perp2,oldPlusDir+perp2);
    interactionLattice.moveOnLattice(oldPos+perp1+perp2,oldPlusDir+perp1+perp2);
}
/**
 * @details Because moves of type MoveLocalBcc must not be used with this
 * feature, this function always throws an exception when called. The function
 * is only implemented for savety purposes.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalBcc
 * @throw std::runtime_error
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::applyMove(const IngredientsType& ing,
							 const MoveLocalBcc& move)
{
	//throw exception in case someone accidentaly uses a bcc-BFM move with this feature
	std::stringstream errormessage;
	errormessage<<"FeatureNNInteractionSc::applyMove(...):\n";
	errormessage<<"attempting to use bcc-BFM move, which is not allowed\n";
	throw std::runtime_error(errormessage.str());
}

/**
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::synchronize(IngredientsType& ingredients)
{
    //refill the lattice with interactionTag tags
    //caution: this overwrites, what is currently written on the lattice
    fillLattice(ingredients);
}

/**
 * @details If not compiled with DEBUG flag this function only returns the content
 * of the lookup table probabilityLookup. If compiled with DEBUG flag it checks
 * that the interactionTag tags typeA, typeB are within the allowed range.
 * @param typeA monomer interactionTag type in range [1,255]
 * @param typeB monomer interactionTag type in range [1,255]
 * @throw std::runtime_error In debug mode, if types are not in range [1,255]
 **/
inline double FeatureNNInteractionSc::getProbabilityFactor(int32_t typeA,
									     int32_t typeB) const
{
#ifdef DEBUG
  //extra checks only in debug mode, because this is very frequently called
  //and this costs performance
  if(typeA<0 || typeA>255 || typeB<0 || typeB>255){
    std::stringstream errormessage;
    errormessage<<"***FeatureNaNInteractionSc::getInteraction(typeA,typeB)***\n";
    errormessage<<"probability undefined between types "<<typeA<<" and "<<typeB<<std::endl;
    errormessage<<"types are out of the allowed range";
    throw std::runtime_error(errormessage.str());
  }
#endif /*DEBUG*/
  	return probabilityLookup[typeA][typeB];
}

/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * - !contactInteraction
 * .
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileReader File importer for the bfm-file
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::exportRead(FileImport< IngredientsType >& fileReader)
{
	fileReader.registerRead("!nn_interaction",new ReadNNInteraction<IngredientsType>(fileReader.getDestination()));
	fileReader.registerRead("!interactionTag",new ReadInteractionTags<IngredientsType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * - !contact_interaction
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureNNInteractionSc::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
	fileWriter.registerWrite("!nn_interaction",new WriteNNInteraction<IngredientsType>(fileWriter.getIngredients_()));
	fileWriter.registerWrite("!interactionTag",new WriteInteractionTags<IngredientsType>(fileWriter.getIngredients_()));
}

/**
 * @details occupies the lattice with the interactionTag tags of the monomers
 * as this is required to determine the contact interactions in this feature.
 * An additional check is performed asserting that the tags are in the range [1,255]
 *
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @throw std::runtime_error In case a monomer has interactionTag tag not in [1,255]
 **/
template<class IngredientsType>
void FeatureNNInteractionSc::fillLattice(IngredientsType& ingredients)
{
	interactionLattice.setupLattice(ingredients.getBoxX(),ingredients.getBoxY(),ingredients.getBoxZ());
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    for(size_t n=0;n<molecules.size();n++)
    {
        VectorInt3 pos=molecules[n];
		uint8_t interactionTag(molecules[n].getInteractionTag());

		if(interactionTag!=molecules[n].getInteractionTag()){
			std::stringstream errormessage;
			errormessage<<"***FeatureNNInteractionSc::fillLattice()***\n";
			errormessage<<"type "<<interactionTag<<" is out of the allowed range";
			throw std::runtime_error(errormessage.str());
		}
        interactionLattice.setLatticeEntry(pos,interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(1,0,0),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(0,1,0),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(1,1,0),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(0,0,1),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(1,0,1),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(0,1,1),interactionTag);
        interactionLattice.setLatticeEntry(pos+VectorInt3(1,1,1),interactionTag);
    }
}


/**
 * @details The function calculates the factor for the acceptance probability
 * for the local move given as argument. The calculation is based on the lattice
 * entries in the vicinity of the monomer to be moved. If the move is accepted,
 * 12 new contacts can potentially be made, and 12 contacts are lost. Thus a number
 * of 24 lattice positions around the monomer have to be checked.
 *
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move reference to the local move for which the calculation is performed
 * @return acceptance probability factor for the move arising from nearest neighbor contacts
 **/
template<class IngredientsType>
double FeatureNNInteractionSc::calculateAcceptanceProbability(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{

    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();

    double prob=1.0;
    int32_t monoType=(int32_t)(ingredients.getMolecules()[move.getIndex()].getInteractionTag());

    /*get two directions perpendicular to vector directon of the move*/
    VectorInt3 perp1,perp2;
    /* first perpendicular direction is either (0 1 0) or (1 0 0)*/
    int32_t x1=((direction.getX()==0) ? 1 : 0);
    int32_t y1=((direction.getX()!=0) ? 1 : 0);
    perp1.setX(x1);
    perp1.setY(y1);
    perp1.setZ(0);

    /* second perpendicular direction is either (0 0 1) or (0 1 0)*/
    int32_t y2=((direction.getZ()==0) ? 0 : 1);
    int32_t z2=((direction.getZ()!=0) ? 0 : 1);
    perp2.setX(0);
    perp2.setY(y2);
    perp2.setZ(z2);

    //the probability is calculated by going through all possible lattice sites
    //at which the contacts may have changed. At every site the type of the
    //monomer sitting there is retrieved from the lattice. the additional
    //factor for the probability (exp(-deltaE/kT)) is retrieved from the
    //lookup using getProbabilityFactor. For new contacts this factor is multiplied
    //with the probability, for contacts taken away the probability is devided.
    VectorInt3 actual=oldPos;

    //first check front,i.e newly acquired contacts
    if(direction.getX()>0 || direction.getY()>0 || direction.getZ()>0) actual+=direction;
    actual+=direction;

    actual-=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp2+direction;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));

    //now check back side (contacts taken away)
    double prob_div=1.0;
    actual=oldPos;
    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) actual-=direction;
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual=actual+perp2-direction;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(interactionLattice.getLatticeEntry(actual)));

    prob/=prob_div;
    return prob;

}

/**
 * @param typeA monomer interactionTag tag in range [1,255]
 * @param typeB monomer interactionTag tag in range [1,255]
 * @param interaction energy between typeA and typeB
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 **/
void FeatureNNInteractionSc::setNNInteraction(int32_t typeA,
								     int32_t typeB,
								     double energy)
{
    if(0<typeA && typeA<=255 && 0<typeB && typeB<=255)
	{
        interactionTable[typeA][typeB]=energy;
        interactionTable[typeB][typeA]=energy;
        probabilityLookup[typeA][typeB]=exp(-energy);
        probabilityLookup[typeB][typeA]=exp(-energy);
        std::cout<<"set interation between types ";
		std::cout<<typeA<<" and "<<typeB<<" to "<<energy<<"kT\n";
	}
    else
	{
		std::stringstream errormessage;
		errormessage<<"FeatureNNInteractionSc::setNNInteraction(typeA,typeB,energy).\n";
		errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
		throw std::runtime_error(errormessage.str());
	}
}

/**
 * @param typeA monomer interactionTag tag in range [1,255]
 * @param typeB monomer interactionTag tag in range [1,255]
 * @throw std::runtime_error In case typeA or typeB exceed range [1,255]
 * @return interaction energy per nearest neighbor contact for typeA,typeB
 **/
double FeatureNNInteractionSc::getNNInteraction(int32_t typeA,
								       int32_t typeB) const
{

    if(0<typeA && typeA<=255 && 0<typeB && typeB<=255)
        return interactionTable[typeA][typeB];
    else
    {
		std::stringstream errormessage;
		errormessage<<"FeatureNNInteractionSc::getNNInteraction(typeA,typeB).\n";
		errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
		throw std::runtime_error(errormessage.str());
    }

}
#endif /*FEATURE_CONTACT_INTERACTION_H*/
