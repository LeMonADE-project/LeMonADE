#ifndef PAIR_POTENTIALS_H
#define PAIR_POTENTIALS_H

#include<cmath>

#include <LeMonADE/utility/DistanceCalculation.h>

class PairPotentialSquare{
public:
	PairPotentialSquare():rangeSquared(0)
	{
		energy.resize(10,std::vector<double>(10,0.0));
		probability.resize(10,std::vector<double>(10,1.0));
		
	}
	
	void setInteraction(int32_t typeA,int32_t typeB,double nrg)
	{
		energy.at(typeA-1).at(typeB-1)=nrg;
		energy.at(typeB-1).at(typeA-1)=nrg;
	}
	
	double getInteraction(int32_t typeA,int32_t typeB) const {return energy.at(typeA-1).at(typeB-1);}
	
	void setRange(double r){rangeSquared=int32_t(r*r);}
	double getRange() const {return std::sqrt(rangeSquared);}
	
	template<class IngredientsType> double getEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		
		if(distVector*distVector<=rangeSquared){
			return energy.at(typeA-1).at(typeB-1);
		}
		
		else return 0.0;
	}
	
	template<class IngredientsType> double getExpEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		
		if(distVector*distVector<=rangeSquared){
			return probability.at(typeA-1).at(typeB-1);
		}
		
		else return 1.0;
		
	}
	
	template<class IngredientsType> void synchronize(const IngredientsType& ing)
	{
		std::cout<<"PairPotentialSquareNNI range set to sqrt("<<rangeSquared<<")\n";
		for(size_t n=0;n<energy.size();n++){
			for(size_t m=0;m<energy.size();m++){
				if(energy[n][m]!=0.0) probability[n][m]=std::exp(-energy[n][m]);
			}
		}
	}
	
private:
	std::vector<std::vector<double> > energy;
	std::vector<std::vector<double> > probability;
	
	int32_t rangeSquared;
	
	
};

class PairPotentialSquareEV{
public:
	PairPotentialSquareEV():rangeSquared(0)
	{
		energy.resize(10,std::vector<double>(10,0.0));
		probability.resize(10,std::vector<double>(10,1.0));
		
	}
	
	void setInteraction(int32_t typeA,int32_t typeB,double nrg)
	{
		energy.at(typeA-1).at(typeB-1)=nrg;
		energy.at(typeB-1).at(typeA-1)=nrg;
	}
	
	double getInteraction(int32_t typeA,int32_t typeB) const {return energy.at(typeA-1).at(typeB-1);}
	
	void setRange(double r){rangeSquared=int32_t(r*r);}
	double getRange() const {return std::sqrt(rangeSquared);}
	
	template<class IngredientsType> double getEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		int32_t dist2=distVector*distVector;
		
		if(dist2<4){
			return std::numeric_limits< double >::infinity();
		}
		else if(dist2<=rangeSquared){
			return energy.at(typeA-1).at(typeB-1);
		}
		
		else return 0.0;
	}
	
	template<class IngredientsType> double getExpEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		int32_t dist2=distVector*distVector;
		
		if(dist2<4){
			return 0.0;
		}
		else if(dist2<=rangeSquared){
			return probability.at(typeA-1).at(typeB-1);
		}
		
		else return 1.0;
		
	}
	
	template<class IngredientsType> void synchronize(const IngredientsType& ing)
	{
		std::cout<<"PairPotentialSquareNNI range set to sqrt("<<rangeSquared<<")\n";
		for(size_t n=0;n<energy.size();n++){
			for(size_t m=0;m<energy.size();m++){
				if(energy[n][m]!=0.0) probability[n][m]=std::exp(-energy[n][m]);
			}
		}
	}
	
private:
	std::vector<std::vector<double> > energy;
	std::vector<std::vector<double> > probability;
	
	int32_t rangeSquared;
	
	
};

class PairPotentialSquareNNI{
public:
	PairPotentialSquareNNI()
	:rangeSquared(0),nnEnergy(0.0)
	,expNNEnergy(1.0),exp2NNEnergy(1.0),exp4NNEnergy(1.0)
	{
		energy.resize(10,std::vector<double>(10,0.0));
		probability.resize(10,std::vector<double>(10,1.0));
		
	}
	
	void setInteraction(int32_t typeA,int32_t typeB,double nrg)
	{
		energy.at(typeA-1).at(typeB-1)=nrg;
		energy.at(typeB-1).at(typeA-1)=nrg;
	}
	
	double getInteraction(int32_t typeA,int32_t typeB) const {return energy.at(typeA-1).at(typeB-1);}
	
	void setRange(double r){rangeSquared=int32_t(r*r);}
	double getRange() const {return std::sqrt(rangeSquared);}
	
	void setNNRepulsion(double nrg){nnEnergy=nrg;}
	double getNNRepultion()const{return nnEnergy;}
	
	template<class IngredientsType> double getEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		int32_t dist2=distVector*distVector;
		
		if(dist2<4){
			return std::numeric_limits< double >::infinity();
		}
		else if(dist2==4){
			return 4*nnEnergy;
		}
		else if(dist2==5){
			return 2*nnEnergy;
		}
		else if(dist2==6){
			return nnEnergy;
		}
		else if(dist2<=rangeSquared){
			return energy.at(typeA-1).at(typeB-1);
		}
		
		else return 0.0;
	}
	
	template<class IngredientsType> double getExpEnergy(VectorInt3 posA,VectorInt3 posB, int32_t typeA, int32_t typeB, const IngredientsType& ing) const
	{
		VectorInt3 distVector=Lemonade::calcDistanceVector3D(posA,posB,ing);
		int32_t dist2=distVector*distVector;
		
		if(dist2<4){
			return 0.0;
		}
		else if(dist2==4){
			return exp4NNEnergy;
		}
		else if(dist2==5){
			return exp2NNEnergy;
		}
		else if(dist2==6){
			return expNNEnergy;
		}
		else if(dist2<=rangeSquared){
			return probability.at(typeA-1).at(typeB-1);
		}
		
		else return 1.0;
		
	}
	
	template<class IngredientsType> void synchronize(const IngredientsType& ing)
	{
		std::cout<<"PairPotentialSquareNNI range set to sqrt("<<rangeSquared<<")\n";
		for(size_t n=0;n<energy.size();n++){
			for(size_t m=0;m<energy.size();m++){
				if(energy[n][m]!=0.0) probability[n][m]=std::exp(-energy[n][m]);
			}
		}
		expNNEnergy=std::exp(-nnEnergy);
		exp2NNEnergy=std::exp(-2.0*nnEnergy);
		exp4NNEnergy=std::exp(-4.0*nnEnergy);
	}
	
private:
	std::vector<std::vector<double> > energy;
	std::vector<std::vector<double> > probability;
	
	int32_t rangeSquared;
	double nnEnergy;
	double expNNEnergy;
	double exp2NNEnergy;
	double exp4NNEnergy;
};
#endif /*PAIR_POTENTIALS_H*/