#ifndef FEATURE_PAIR_INTERACTION_H
#define FEATURE_PAIR_INTERACTION_H


#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddScMonomer.h>


template<class Potential>
class FeaturePairInteraction:public Feature
{

public:
	typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;
	typedef LOKI_TYPELIST_2(FeatureBox,FeatureAttributes) required_features_front;

	FeaturePairInteraction()
	:verletSkin(0.0),useVerletList(false),verletListIsOutdated(true){};
	
	virtual ~FeaturePairInteraction(){}
	
	//configureing the potential ///////////////////////////////////////////
	Potential& modifyPairPotential(){return potential;}
	const Potential& getPairPotential() const {return potential;}
	
	// configure the verlet list ///////////////////////////////////////////
	void setUseVerletList(bool val){useVerletList=val;}
	void setVerletSkin(double s){verletSkin=s;}
	double getVerletSkin() const {return verletSkin;} 
	
	//for debugging
	size_t getVerletListLength() const {return verletList.size();}
	
	//check- and apply-functions for different types of moves /////////////
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;
	
	template<class IngredientsType,class MoveType>
	bool checkMove(const IngredientsType& ingredients, MoveLocalBase<MoveType>& move);
	
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddScMonomer& move) const;
	
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move){}
	
	template<class IngredientsType,class MoveType>
	void applyMove(IngredientsType& ingredients,const MoveLocalBase<MoveType>& move);
	
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	
	//export read and write functionality //////////////////////////////////
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader){}
	
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const{}
	
	
private:
	
	// stuff for potential /////////////////////////////////////////////////
	Potential potential;
	
	// stuff for verlet list ///////////////////////////////////////////////
	template<class IngredientsType>
	void updateVerletList(const IngredientsType& ingredients);
	
	bool useVerletList;
	bool verletListIsOutdated;
	double verletSkin;
	
	std::vector<size_t> verletList;
	std::map<size_t,size_t> verletNeigborsMap;
	std::map<size_t,size_t> verletIndexMap;
	std::map<size_t,VectorInt3> movedDistanceSinceListUpdate;
	
	
	
};


///////////////////////////////////////////////////////////////////////////////


template<class Potential>
template<class IngredientsType>
bool FeaturePairInteraction<Potential>::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}

template<class Potential>
template<class IngredientsType>
bool FeaturePairInteraction<Potential>::checkMove(const IngredientsType& ingredients, const MoveAddScMonomer& move) const
{
	VectorInt3 myPos=move.getPosition();
	int32_t myType=move.getType();
	
	for(size_t n=0;n<ingredients.getMolecules().size();n++){
		VectorInt3 monoPos=ingredients.getMolecules()[n];
		int32_t monoType=ingredients.getMolecules()[n].getAttributeTag();
		if(potential.getEnergy(myPos,monoPos,myType,monoType,ingredients)==std::numeric_limits< double >::infinity()) return false;
	}
	
	//if still here we reuturn true
	return true;
}

template<class Potential>
template<class IngredientsType,class MoveType>
bool FeaturePairInteraction<Potential>::checkMove(const IngredientsType& ingredients, MoveLocalBase<MoveType>& move)
{
	
	if (useVerletList==false)
	{
		
		double exp_oldEnergy=1.0;
		double exp_newEnergy=1.0;
		
		size_t myIdx=move.getIndex();
		VectorInt3 myPos=ingredients.getMolecules()[myIdx];
		VectorInt3 myPosPlusDelta=ingredients.getMolecules()[myIdx]+move.getDir();
		int32_t myType=ingredients.getMolecules()[myIdx].getAttributeTag();
		
		for(size_t monoIdx=0;monoIdx<ingredients.getMolecules().size();monoIdx++)
		{
			
			if(monoIdx!=myIdx)
			{
				int32_t  monoType=ingredients.getMolecules()[monoIdx].getAttributeTag();
				
 				VectorInt3 monoPos=ingredients.getMolecules()[monoIdx];
				
				if(potential.getEnergy(myPosPlusDelta,monoPos,myType,monoType,ingredients)==std::numeric_limits< double >::infinity()) return false;
			
				exp_oldEnergy*=potential.getExpEnergy(myPos,monoPos,myType,monoType,ingredients);
				exp_newEnergy*=potential.getExpEnergy(myPosPlusDelta,monoPos,myType,monoType,ingredients);
				
			}
		}
		
		double moveProb=exp_newEnergy/exp_oldEnergy;
		move.multiplyProbability( moveProb);
		return true;
	}
	else //if useVerletList==true
	{
		if(verletListIsOutdated==true) 
		{
			updateVerletList(ingredients);
			verletListIsOutdated=false;
		}
		
		double exp_oldEnergy=1.0;
		double exp_newEnergy=1.0;
		
		double oldEnergy=0.0;
		double newEnergy=0.0;
		size_t myIdx=move.getIndex();
		VectorInt3 myPos=ingredients.getMolecules()[myIdx];
		VectorInt3 myPosPlusDelta=ingredients.getMolecules()[myIdx]+move.getDir();	
		int32_t myType=ingredients.getMolecules()[myIdx].getAttributeTag();

		
		size_t numberOfNeighbors=verletNeigborsMap[myIdx];
		size_t verletListIndex=verletIndexMap[myIdx];
		size_t neighborIndex;
		int32_t neighborType;
		
		for(size_t n=0;n<numberOfNeighbors;n++)
		{
			neighborIndex=verletList[verletListIndex+n];
			neighborType=ingredients.getMolecules()[neighborIndex].getAttributeTag();
			
			VectorInt3 neighborPos=ingredients.getMolecules()[neighborIndex];
			
			if(potential.getEnergy(myPosPlusDelta,neighborPos,myType,neighborType,ingredients)==std::numeric_limits< double >::infinity()) return false;
			
			exp_oldEnergy*=potential.getExpEnergy(myPos,neighborPos,myType,neighborType,ingredients);
			
			exp_newEnergy*=potential.getExpEnergy(myPosPlusDelta,neighborPos,myType,neighborType,ingredients);
			
		}
	
		double moveProb=exp_newEnergy/exp_oldEnergy;
		
		move.multiplyProbability( moveProb);
		return true;
	}
}


template<class Potential>
template<class IngredientsType,class MoveType>
void FeaturePairInteraction<Potential>::applyMove(IngredientsType& ing, const MoveLocalBase<MoveType>& move)
{
	
	if(useVerletList==true)
	{	  
		
		movedDistanceSinceListUpdate[move.getIndex()]+=move.getDir();
		
		double distanceMovedSquared=double(movedDistanceSinceListUpdate[move.getIndex()]*movedDistanceSinceListUpdate[move.getIndex()]);
		if(distanceMovedSquared >= verletSkin*verletSkin/4.0)
		{
			verletListIsOutdated=true;
		}
		
	}
}

template<class Potential>
template<class IngredientsType>
void FeaturePairInteraction<Potential>::synchronize(IngredientsType& ingredients)
{
	
	potential.synchronize(ingredients);
	
	
	if(useVerletList==true)
	{
		std::cout<<"FeaturePairInteraction::synchronize()...start setting up Verlet list...\n"; 
		updateVerletList(ingredients);
		verletListIsOutdated=false;
		std::cout<<"FeaturePairInteraction::synchronize()...done setting up Verlet list...\n"; 
	}
	else
	{
		verletIndexMap.clear();
		verletNeigborsMap.clear();
		movedDistanceSinceListUpdate.clear();
		verletList.clear();
	}
	
	for(size_t n=0;n<ingredients.getMolecules().size();n++){
		for(size_t m=0;m<n;m++){
			VectorInt3 posN=ingredients.getMolecules()[n];
			VectorInt3 posM=ingredients.getMolecules()[m];
			int32_t attN=ingredients.getMolecules()[n].getAttributeTag();
			int32_t attM=ingredients.getMolecules()[m].getAttributeTag();
			if(potential.getExpEnergy(posN,posM,attN,attM,ingredients)==0.0){
				std::stringstream errormessage;
				errormessage<<"FeaturePairInteraction: invalid positions\n";
				errormessage<<"mono a "<<n<<" mono pos "<<posN<<" of type "<<attN <<std::endl;
				errormessage<<"mono b "<<m<<" mono pos "<<posM<<" of type "<<attM <<std::endl;
				
				throw std::runtime_error(errormessage.str());
			}
		}
		
	}
	
}


template<class Potential>
template<class IngredientsType>
void FeaturePairInteraction<Potential>::updateVerletList(const IngredientsType& ingredients)
{
	
	//reset all list related containers
	verletIndexMap.clear();
	verletNeigborsMap.clear();
	movedDistanceSinceListUpdate.clear();
	verletList.clear();
	
	
	VectorInt3 distanceVector, posN,posM;
	
	double verletListRadiusSquared=(potential.getRange()+verletSkin)*(potential.getRange()+verletSkin);
	
	size_t neighborCount=0;
	
	for(size_t n=0;n<ingredients.getMolecules().size();n++)
	{
		posN=ingredients.getMolecules()[n];
	
		for(size_t m=0;m<ingredients.getMolecules().size();m++)
		{
			if(m!=n)
			{
				
				posM=ingredients.getMolecules()[m];
				
				distanceVector=Lemonade::calcDistanceVector3D(posN,posM,ingredients);
		
				if(distanceVector*distanceVector <= verletListRadiusSquared)
				{
		
					verletList.push_back(m);
					neighborCount++;
				}
		
			}
		}
		
		verletNeigborsMap.insert(std::make_pair(n,neighborCount));
		
		if(neighborCount>0)
		{
			verletIndexMap.insert(std::make_pair(n,verletList.size()-neighborCount));
		}	
		
		neighborCount=0;
		
		movedDistanceSinceListUpdate.insert(std::make_pair(n,VectorInt3(0,0,0)));
	}
	
}



#endif // FEATURE
