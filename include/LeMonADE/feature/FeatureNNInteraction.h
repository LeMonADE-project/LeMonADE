#ifndef FEATURE_NN_INTERACTION_H
#define FEATURE_NN_INTERACTION_H

#include<limits>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include <LeMonADE/feature/NNInteractionReadWrite.h>



/******************************************************************************
 *this feature implements nearest neighbour contact interactions.
 *in the bfm file, as well as in the code, the interactions are defined as
 *for each nearest neighbor contact between two types A,B in kT.
 *bfm-command example:
 *#!interaction 1 2 1.6
 *defines an interaction of 1.6 kT per contact for monomer types 1,2
 *
 *as default the feature is implemented for use with ExcludedVolumeSC (writes
 *8 points per monomer on the lattice), but code for ExcludedVolume (writing only
 *1 point) is also there. just un-comment the marked regions.
 *
 *the feature provides two #! bfm-commands: #!interaction_base and #!interaction
 *
 *****************************************************************************/
template<template<typename> class FeatureLatticeType>
class FeatureNNInteractionSc:public Feature
{

private:
  typedef uint8_t lattice_value_type;

  double interactionTable[256][256];
    double probabilityLookup[256][256];

    template<class IngredientsType>
    double calculateAcceptanceProbability(const IngredientsType& ingredients, const MoveLocalSc& move) const;

    template<class IngredientsType>
    void fillLattice(IngredientsType& ingredients);

    double getProbabilityFactor(int32_t typeA,int32_t typeB) const;


public:

    FeatureNNInteractionSc();
    ~FeatureNNInteractionSc();

    
    //and this feature re-initializes it (overwriting the values from excluded volume)
    typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

    //FeatureExcludedVulumeSc needs to be in front, because FeatureNNInteractionSc
    //re-initializes the lattice and overwrites what FeatureExcludedVolumeSc has written
    typedef LOKI_TYPELIST_2(FeatureAttributes,FeatureExcludedVolumeSc<FeatureLatticeType<lattice_value_type> >) required_features_front;


    //check- and apply-functions for different types of moves
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;

    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;

    template<class IngredientsType>
    void applyMove(IngredientsType& ing, const MoveBase& move){}

    template<class IngredientsType>
    void applyMove(IngredientsType& ing, const MoveAddScMonomer& move);

    //synchronize initializes the lattice with the attribute-tags of the monomers
    template<class IngredientsType>
    void synchronize(IngredientsType& ingredients);

    //add interaction between two types of monomers
    void setContactInteraction(int32_t typeA,int32_t typeB,double energy);

    //returns the interaction between two types of monomers in units of interactionBase
    double getContactInteraction(int32_t typeA,int32_t typeB) const;

    //export read and write functionality
    template <class IngredientsType>
    void exportRead(FileImport <IngredientsType>& fileReader);

    template <class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

};





/*******************************************************************************
 *implementation of template and inline members
 *(non-template members in FeatureNNInteraction.cpp)
 ******************************************************************************/
template<template<typename> class LatticeClassType>
inline double FeatureNNInteractionSc<LatticeClassType>::getProbabilityFactor(int32_t typeA,int32_t typeB) const
{
#ifdef DEBUG

  if(typeA<0 || typeA>255 || typeB<0 || typeB>255){
    std::stringstream errormessage;
    errormessage<<"***FeatureNaNInteractionSc::getInteraction(typeA,typeB)***\n";
    errormessage<<"probability undefined between types "<<typeA<<" and "<<typeB<<std::endl;
    errormessage<<"types are out of the allowed range";

    throw std::runtime_error(errormessage.str());
}
#endif

	return probabilityLookup[typeA][typeB];

}


/*******************************************************************************
 *exportRead and exportWrite: export of command #!nn_interaction
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureNNInteractionSc<LatticeClassType>::exportRead(FileImport< IngredientsType >& fileReader)
{
  typedef FeatureNNInteractionSc<LatticeClassType> my_type;
  fileReader.registerRead("!nn_interaction",new ReadInteraction<my_type>(*this));
}

template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureNNInteractionSc<LatticeClassType>::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  typedef FeatureNNInteractionSc<LatticeClassType> my_type;
  fileWriter.registerWrite("!nn_interaction",new WriteInteraction<my_type>(*this));
}


/*******************************************************************************
 *checkMove for unknown moves: does nothing
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureNNInteractionSc<LatticeClassType>::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
    return true;
}

/*******************************************************************************
 *checkMove for MoveLocalSc: gets the transition probability and adds it to
 *the move. the metropolis step is performed later in FeatureBoltzmann (in
 *case other features also add some probability)
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
bool FeatureNNInteractionSc<LatticeClassType>::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
    double prob=calculateAcceptanceProbability(ingredients,move);
    move.multiplyProbability(prob);
    return true;
}



/*******************************************************************************
 *synchronize:
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureNNInteractionSc<LatticeClassType>::synchronize(IngredientsType& ingredients)
{

    //refill the lattice with attribute tags
    //caution: this overwrites, what is currently written on the lattice
    fillLattice(ingredients);

}


/*******************************************************************************
 *applyMove for MoveAddScMonomer: updates the lattice with attribute tag of the
 *added monomer
 *works only when 8 positions are written on the lattice per monomer. for
 *the other case use the function in the lower part of this file
 ******************************************************************************/
//////////////// version for use with FeatureExcludedVolumeSc //////////////////
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureNNInteractionSc<LatticeClassType>::applyMove(IngredientsType& ing, const MoveAddScMonomer& move)
{
    //get the position and attribute tag of the monomer to be inserted
    VectorInt3 pos=move.getPosition();
    VectorInt3 dx(1,0,0);
    VectorInt3 dy(0,1,0);
    VectorInt3 dz(0,0,1);
    int32_t type=move.getType();

    //update lattice
    ing.setLatticeEntry(pos,type);
    ing.setLatticeEntry(pos+dx,type);
    ing.setLatticeEntry(pos+dy,type);
    ing.setLatticeEntry(pos+dx+dy,type);
    ing.setLatticeEntry(pos+dz,type);
    ing.setLatticeEntry(pos+dz+dx,type);
    ing.setLatticeEntry(pos+dz+dy,type);
    ing.setLatticeEntry(pos+dz+dx+dy,type);
}

/*******************************************************************************
 *fillLattice: fills the lattice with attribute tags of the monomers.
 *this is the version for use with ExcludedVolumeSC (it writes 8 positions
 *per monomer)
 *if you want to use a version where only one position is written per monomer,
 *comment this function and uncomment the appropriate one in the section at the
 *end of this file.
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
void FeatureNNInteractionSc<LatticeClassType>::fillLattice(IngredientsType& ingredients)
{
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    for(size_t n=0;n<molecules.size();n++)
    {
        VectorInt3 pos=molecules[n];
	lattice_value_type attribute=lattice_value_type(molecules[n].getAttributeTag());

	if(int32_t(attribute)!=molecules[n].getAttributeTag()){
	  std::stringstream errormessage;
	  errormessage<<"***FeatureNaNInteractionSc::fillLattice()***\n";
	  errormessage<<"type "<<attribute<<" is out of the allowed range";

	  throw std::runtime_error(errormessage.str());
	}

        ingredients.setLatticeEntry(pos,attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,0),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,0,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,1),attribute);
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,1),attribute);

    }

}

/*******************************************************************************
 *calculateAcceptanceProb:
 *calculates the probability of acceptance for the move based on change of contacts
 *with surrounding monomers
 *this function works only when 8 positions are written per monomer. For other
 *case comment it out and use the function in the lower section of this file.
 ******************************************************************************/
template<template<typename> class LatticeClassType>
template<class IngredientsType>
double FeatureNNInteractionSc<LatticeClassType>::calculateAcceptanceProbability(
    const IngredientsType& ingredients,
    const MoveLocalSc& move) const
{
    
    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();

    double prob=1.0;
    int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

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
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+direction;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    //now check back side (contacts taken away)
    double prob_div=1.0;
    actual=oldPos;
    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) actual-=direction;
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2-direction;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,int32_t(ingredients.getLatticeEntry(actual)));

    prob/=prob_div;
    return prob;

}


/*******************************************************************************
 *constructor: initializes the lookup arrays
 ******************************************************************************/
template<template<typename> class LatticeClassType>
FeatureNNInteractionSc<LatticeClassType>::FeatureNNInteractionSc()
{
    for(size_t n=0;n<256;n++)
    {
        for(size_t m=0;m<256;m++)
        {
            interactionTable[m][n]=0.0;
            probabilityLookup[m][n]=1.0;
        }
    }
}

/*******************************************************************************
 *destructor
 ******************************************************************************/
template<template<typename> class LatticeClassType>
FeatureNNInteractionSc<LatticeClassType>::~FeatureNNInteractionSc(){}

template<template<typename> class LatticeClassType>
void FeatureNNInteractionSc<LatticeClassType>::setContactInteraction(int32_t typeA,int32_t typeB,double energy)
{
    if(0<typeA && typeA<=255 && 0<typeB && typeB<=255)
    {
        interactionTable[typeA][typeB]=energy;
        interactionTable[typeB][typeA]=energy;
        probabilityLookup[typeA][typeB]=exp(-energy);
        probabilityLookup[typeB][typeA]=exp(-energy);
        std::cout<<"set interation between types "<<typeA<<" and "<<typeB<<" to "<<energy<<"kT\n";
    }
    else
    {
      std::stringstream errormessage;
      errormessage<<"FeatureNNInteractionSc::setContactInteraction(typeA,typeB,energy).\n";
      errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
      throw std::runtime_error(errormessage.str());
    }
}

//returns the interaction between two types of monomers in units of interactionBase
//int32_t FeatureNNInteraction::getInteraction(uint32_t typeA,uint32_t typeB) const
template<template<typename> class LatticeClassType>
double FeatureNNInteractionSc<LatticeClassType>::getContactInteraction(int32_t typeA,
								       int32_t typeB) const
{

    if(0<typeA && typeA<=255 && 0<typeB && typeB<=255)
        return interactionTable[typeA][typeB];
    else
    {
      std::stringstream errormessage;
      errormessage<<"FeatureNNInteractionSc::getContactInteraction(typeA,typeB).\n";
      errormessage<<"typeA "<<typeA<<" typeB "<<typeB<<": Types out of range\n";
      throw std::runtime_error(errormessage.str());
    }

}


#endif /*FEATURE_NN_INTERACTION_H*/
