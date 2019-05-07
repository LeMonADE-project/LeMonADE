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
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>

class FeatureSpringPotentialTwoGroups:public Feature
{
public:
	FeatureSpringPotentialTwoGroups(){};
	virtual ~FeatureSpringPotentialTwoGroups(){};
	
	typedef LOKI_TYPELIST_2(FeatureAttributes<>, FeatureBoltzmann) required_features_back;
	
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move)const;
	
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;
	
	double getEquilibriumLength() const{
		return equilibrium_length;
	}

	void setEquilibriumLength(double equilibriumLength) {
		equilibrium_length = equilibriumLength;
	}

	double getSpringConstant() const{
		return spring_constant;
	}

	void setSpringConstant(double springConstant){
		spring_constant = springConstant;
	}

	int32_t getAffectedMonomerType0() const {
		return affectedMonomerType0;
	}

	void setAffectedMonomerType0(int32_t affectedMonomerType0) {
		this->affectedMonomerType0 = affectedMonomerType0;
	}

	int32_t getAffectedMonomerType1() const {
		return affectedMonomerType1;
	}
	void setAffectedMonomerType1(int32_t affectedMonomerType1) {
		this->affectedMonomerType1 = affectedMonomerType1;
	}

	template<class IngredientsType>
	VectorDouble3 getGroupCenterOfMass(const IngredientsType& ingredients,const std::vector<uint32_t>& group) const;

	template<class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);
	
	template<class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

private:
	//! equlibrium length r0 in harmonic potential V(r)=k/2(r-r0)^2
	double equilibrium_length;

	//! spring constant k in harmonic potential V(r)=k/2(r-r0)^2
	double spring_constant;
	
	//attribute tag of the monomer type affected by the potential
	int32_t affectedMonomerType0;

	//contains the indices of the monomers of type affectedMonomerType
	std::vector<uint32_t> affectedMonomerGroup0;

	//attribute tag of the monomer type affected by the potential
	int32_t affectedMonomerType1;

	//contains the indices of the monomers of type affectedMonomerType
	std::vector<uint32_t> affectedMonomerGroup1;

};



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
	std::cout << "#!virtual_spring_constant=" << (springConstant) << std::endl;

	ingredients.setSpringConstant(springConstant);
}

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
	std::cout << "#!virtual_spring_length=" << (springLength) << std::endl;

	ingredients.setEquilibriumLength(springLength);
}

template < class IngredientsType>
class ReadVirtualSpringType0: public ReadToDestination<IngredientsType>
{
public:
	ReadVirtualSpringType0(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadVirtualSpringType0(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadVirtualSpringType0<IngredientsType>::execute()
{
	std::cout<<"reading VirtualSpringType0...";

	int32_t springType0 = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	springType0 = atoi(line.c_str());
	std::cout << "#!virtual_spring_type0=" << (springType0) << std::endl;

	ingredients.setAffectedMonomerType0(springType0);
}

template < class IngredientsType>
class ReadVirtualSpringType1: public ReadToDestination<IngredientsType>
{
public:
	ReadVirtualSpringType1(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadVirtualSpringType1(){}
  virtual void execute();
};



template<class IngredientsType>
void ReadVirtualSpringType1<IngredientsType>::execute()
{
	std::cout<<"reading VirtualSpringType1...";

	int32_t springType1 = 0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	springType1 = atoi(line.c_str());
	std::cout << "#!virtual_spring_type1=" << (springType1) << std::endl;

	ingredients.setAffectedMonomerType1(springType1);
}

template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!virtual_spring_constant", new ReadVirtualSpringConstant<FeatureSpringPotentialTwoGroups>(*this));
    fileReader.registerRead("#!virtual_spring_length", new ReadVirtualSpringLength<FeatureSpringPotentialTwoGroups>(*this));
    fileReader.registerRead("#!virtual_spring_type0", new ReadVirtualSpringType0<FeatureSpringPotentialTwoGroups>(*this));
    fileReader.registerRead("#!virtual_spring_type1", new ReadVirtualSpringType1<FeatureSpringPotentialTwoGroups>(*this));

}

/***************************************************************************************/


template <class IngredientsType>
class WriteVirtualSpringConstant:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringConstant(const IngredientsType& i):AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringConstant(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringConstant<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!virtual_spring_constant=" << (this->getSource().getSpringConstant()) << std::endl<< std::endl;
}

/*****************************************************************/
/**
 * @class WriteVirtualSpringLength
 *
 * @brief Handles BFM-File-Write \b #!VirtualSpringLength
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteVirtualSpringLength:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringLength(const IngredientsType& i):AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringLength(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringLength<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!virtual_spring_length=" << (this->getSource().getEquilibriumLength()) << std::endl<< std::endl;
}


/*****************************************************************/
/**
 * @class WriteVirtualSpringType0
 *
 * @brief Handles BFM-File-Write \b #!VirtualSpringType0
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteVirtualSpringType0:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringType0(const IngredientsType& i):AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringType0(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringType0<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!virtual_spring_type0=" << (this->getSource().getAffectedMonomerType0()) << std::endl<< std::endl;
}


template <class IngredientsType>
class WriteVirtualSpringType1:public AbstractWrite<IngredientsType>
{
public:
	WriteVirtualSpringType1(const IngredientsType& i):AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteVirtualSpringType1(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteVirtualSpringType1<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!virtual_spring_type1=" << (this->getSource().getAffectedMonomerType1()) << std::endl<< std::endl;
}



template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
	fileWriter.registerWrite("#!virtual_spring_constant",new WriteVirtualSpringConstant<FeatureSpringPotentialTwoGroups>(*this));
	fileWriter.registerWrite("#!virtual_spring_length",new WriteVirtualSpringLength<FeatureSpringPotentialTwoGroups>(*this));
	fileWriter.registerWrite("#!virtual_spring_type0",new WriteVirtualSpringType0<FeatureSpringPotentialTwoGroups>(*this));
	fileWriter.registerWrite("#!virtual_spring_type1",new WriteVirtualSpringType1<FeatureSpringPotentialTwoGroups>(*this));
}

template<class IngredientsType>
bool FeatureSpringPotentialTwoGroups::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}


template<class IngredientsType>
bool FeatureSpringPotentialTwoGroups::checkMove(const IngredientsType& ingredients, MoveLocalSc& move)const
{
	//Index of moved Monomer is monoIndex
	uint32_t monoIndex=move.getIndex();
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

	//if the moved monomer has the same attribute as the affectedMonomerType,
	//get the z position of the group afer a hypothetical move
	//otherwise return true right away
	int32_t attributeTagMoveObject=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

	VectorDouble3 COM_position_old;
	VectorDouble3 COM_position_not_moved;
	double size_moved_group = 0.0;


	VectorDouble3 projected_move = move.getDir();

	if(attributeTagMoveObject==affectedMonomerType0)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		size_moved_group=affectedMonomerGroup0.size();
	}
	else if(attributeTagMoveObject==affectedMonomerType1)
	{
		COM_position_old=getGroupCenterOfMass(ingredients,affectedMonomerGroup1);
		COM_position_not_moved=getGroupCenterOfMass(ingredients,affectedMonomerGroup0);
		size_moved_group=affectedMonomerGroup1.size();
	}
	else return true;

	//the rest happens only if we have not returned yet, i.e. if distance has been calculated

	VectorDouble3 rel_length_old(COM_position_old-COM_position_not_moved);
	VectorDouble3 rel_length_moved(COM_position_old-COM_position_not_moved+projected_move/size_moved_group);

	// calculate the potential difference
	double dV = spring_constant*(COM_position_old*projected_move/size_moved_group);
	dV += 0.5*spring_constant/(size_moved_group*size_moved_group);
	dV -= spring_constant*(COM_position_not_moved*projected_move/size_moved_group);
	dV += spring_constant*equilibrium_length*(rel_length_old.getLength()-rel_length_moved.getLength());

	//calculate the transition probability
	//Metropolis: zeta = exp (-dV)
	double prob=exp(-dV);

	//std::cout << "prob: " <<  prob << std::endl;
	move.multiplyProbability(prob);

	return true;

}

//finds the indices of the affected monomers and checks if the
//center of mass is in a valid position, i.e. that the potential is not
//infinite
template<class IngredientsType>
void FeatureSpringPotentialTwoGroups::synchronize(IngredientsType& ingredients)
{
	// delete the old content
	affectedMonomerGroup0.clear();
	affectedMonomerGroup1.clear();

	//sort the monomers into groups
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		if(ingredients.getMolecules()[n].getAttributeTag()==affectedMonomerType0)
		{
			affectedMonomerGroup0.push_back(n);
		}

		if(ingredients.getMolecules()[n].getAttributeTag()==affectedMonomerType1)
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
