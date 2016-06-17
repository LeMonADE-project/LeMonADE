#ifndef LATTICE_OCCUPATION_ANALYZER_H
#define LATTICE_OCCUPATION_ANALYZER_H

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>

template<class IngredientsType> 
class LatticeOccupationAnalyzer:public IngredientsAnalyzer<IngredientsType>
{
public:
	LatticeOccupationAnalyzer(const IngredientsType& ing)
		:IngredientsAnalyzer<IngredientsType>(ing)
		{}
	
	virtual ~LatticeOccupationAnalyzer(){}
	
	
	virtual bool execute();
	virtual void initialize();
	virtual void cleanup(){}
};


/************ implementation ******************/

template<class IngredientsType>
bool LatticeOccupationAnalyzer<IngredientsType>::execute()
{
	int32_t boxX=this->getIngredients().getBoxX();
	int32_t boxY=this->getIngredients().getBoxY();
	int32_t boxZ=this->getIngredients().getBoxZ();
	VectorInt3 pos;
	int32_t latticeOccupation=0;
	for(int32_t x=0;x<boxX;x++)
	{
		for(int32_t y=0;y<boxY;y++)
		{
			for(int32_t z=0;z<boxZ;z++)
			{
				pos.setAllCoordinates(x,y,z);
				if(this->getIngredients().getLatticeEntry(pos)>0)
					latticeOccupation++;
			}
		}
	}
	
    std::cout<<"********LatticeOccupationAnalyzer *******************\n";
	std::cout<<"lattice occupation "<<latticeOccupation
		 <<"... number of monomers "<<this->getIngredients().getMolecules().size();
		 
        if(8*this->getIngredients().getMolecules().size()==latticeOccupation)
		std::cout<<"...ok"<<std::endl;
	else
	{
		std::cout<<"...error"<<std::endl;
		throw std::runtime_error("***************** LatticeOccupationAnalyzer:*************************\n lattice occupation does not match number of monomers in the system!!!");
	}
				
}

template<class IngredientsType>
void LatticeOccupationAnalyzer<IngredientsType>::initialize()
{
    std::cout<<"LatticeOccupationAnalyzer:checking initial lattice occupation......";
	execute();
}

#endif
