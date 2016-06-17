#ifndef FEATURE_NN_INTERACTION_H
#define FEATURE_NN_INTERACTION_H


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
class FeatureNNInteraction:public Feature
{
public:

    FeatureNNInteraction();
    ~FeatureNNInteraction();

    //excluded volume needs to be in front, because it pre-initializes the lattice
    //and this feature re-initializes it (overwriting the values from excluded volume)
    typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

    ///////// for use with FeatureExcludedVolumeSc /////////////////////////////
    typedef LOKI_TYPELIST_2(FeatureAttributes,FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<uint8_t> >) required_features_front;
    ///////// end version with Feature ExcludedVolumeSc ////////////////////////

    ///////////for use with FeatureExcludedVolume //////////////////////////////
    //typedef LOKI_TYPELIST_2(FeatureAttributes,FeatureExcludedVolume<uint8_t>) required_features_front;
    ///////// end version with Feature ExcludedVolume //////////////////////////

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
    void setInteraction(uint32_t typeA,uint32_t typeB,double energy);

    //returns the interaction between two types of monomers in units of interactionBase
    double getInteraction(uint32_t typeA,uint32_t typeB) const;

    //export read and write functionality
    template <class IngredientsType>
    void exportRead(FileImport <IngredientsType>& fileReader);

    template <class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

private:

    double interactionTable[11][11];
    double probabilityLookup[11][11];

    template<class IngredientsType>
    double calculateAcceptanceProb(const IngredientsType& ingredients, const MoveLocalSc& move) const;

    template<class IngredientsType>
    void fillLattice(IngredientsType& ingredients);

    //defined inline in the header (turns out to be faster)
    double getProbabilityFactor(uint32_t typeA,uint32_t typeB) const;

//////////// for use with FeatureExcludedVolume uncomment the block below///////
//
//    static const VectorInt3 contactShell_1[24];
//    static const VectorInt3 contactShell_2[24];
//    static const VectorInt3 contactShell_4[6];
///////// end version with Feature ExcludedVolume //////////////////////////////
};





/*******************************************************************************
 *implementation of template and inline members
 *(non-template members in FeatureNNInteraction.cpp)
 ******************************************************************************/
inline double FeatureNNInteraction::getProbabilityFactor(uint32_t typeA,uint32_t typeB) const
{
#if DEBUG
	if(0<typeA && typeA<=10 && 0<typeB && typeB<=10)
#endif


	return probabilityLookup[typeA][typeB];

#if DEBUG
	else
		throw std::runtime_error("***FeatureNNInteraction::getInteraction()...trying to get undefined probability***\n");
#endif

}


/*******************************************************************************
 *exportRead and exportWrite: export of command #!nn_interaction
 ******************************************************************************/
template<class IngredientsType>
void FeatureNNInteraction::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!nn_interaction",new ReadInteraction<FeatureNNInteraction>(*this));
}

template<class IngredientsType>
void FeatureNNInteraction::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
    fileWriter.registerWrite("#!nn_interaction",new WriteInteraction<FeatureNNInteraction>(*this));
}


/*******************************************************************************
 *checkMove for unknown moves: does nothing
 ******************************************************************************/
template<class IngredientsType>
bool FeatureNNInteraction::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
    return true;
}

/*******************************************************************************
 *checkMove for MoveLocalSc: gets the transition probability and adds it to
 *the move. the metropolis step is performed later in FeatureBoltzmann (in
 *case other features also add some probability)
 ******************************************************************************/
template<class IngredientsType>
bool FeatureNNInteraction::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
    double prob=calculateAcceptanceProb(ingredients,move);
    move.multiplyProbability(prob);
    return true;
}



/*******************************************************************************
 *synchronize:
 ******************************************************************************/
template<class IngredientsType>
void FeatureNNInteraction::synchronize(IngredientsType& ingredients)
{

    //refill the lattice with attribute tags
    //caution: this overwrites, what is currently written on the lattice
    fillLattice(ingredients);

}

//////////////// version for use with FeatureExcludedVolumeSc //////////////////

/*******************************************************************************
 *applyMove for MoveAddScMonomer: updates the lattice with attribute tag of the
 *added monomer
 *works only when 8 positions are written on the lattice per monomer. for
 *the other case use the function in the lower part of this file
 ******************************************************************************/
template<class IngredientsType>
void FeatureNNInteraction::applyMove(IngredientsType& ing, const MoveAddScMonomer& move)
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
template<class IngredientsType>
void FeatureNNInteraction::fillLattice(IngredientsType& ingredients)
{
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    for(size_t n=0;n<molecules.size();n++)
    {
        VectorInt3 pos=molecules[n];

        ingredients.setLatticeEntry(pos,molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,0),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,0),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,0),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(0,0,1),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(1,0,1),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(0,1,1),molecules[n].getAttributeTag());
        ingredients.setLatticeEntry(pos+VectorInt3(1,1,1),molecules[n].getAttributeTag());

    }

}

/*******************************************************************************
 *calculateAcceptanceProb:
 *calculates the probability of acceptance for the move based on change of contacts
 *with surrounding monomers
 *this function works only when 8 positions are written per monomer. For other
 *case comment it out and use the function in the lower section of this file.
 ******************************************************************************/
template<class IngredientsType>
double FeatureNNInteraction::calculateAcceptanceProb(const IngredientsType& ingredients,  const MoveLocalSc& move) const
{
    //stupid implementation for now.
    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();



    double prob=1.0;
    uint32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();

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
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+direction;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));

    //now check back side (contacts taken away)
    double prob_div=1.0;
    actual=oldPos;
    if(direction.getX()<0 || direction.getY()<0 || direction.getZ()<0) actual-=direction;
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2+perp1;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual-perp1-perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp1;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual=actual+perp2-direction;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual+=perp1;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));
    actual-=perp2;
    prob_div*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(actual)));

    prob/=prob_div;
    return prob;

}

///////////////// end version with FeatureExcludedVolumeSc ////////////////////


////////////////  version with FeatureExcludedVolume  ////////////////////////////
//
//
// template<class IngredientsType>
// double FeatureNNInteraction::calculateAcceptanceProb(const IngredientsType& ingredients,  const MoveLocalSc& move) const
// {
//    //stupid implementation for now.
//    VectorInt3 oldPos=ingredients.getMolecules()[move.getIndex()];
//    VectorInt3 newPos=oldPos+move.getDir();
//
//    double prob=1.0;
//    double prob_divide=1.0;
//    uint32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();
//
//
//    for(size_t n=0;n<24;n++)
//    {
//        prob*=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(newPos+contactShell_1[n])));
//        prob_divide/=getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(oldPos+contactShell_1[n])));
//    }
//    for(size_t n=0;n<24;n++)
//    {
//        prob*=pow(getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(newPos+contactShell_2[n]))),2.0);
//        prob_divide/=pow(getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(oldPos+contactShell_2[n]))),2.0);
//    }
//    for(size_t n=0;n<6;n++)
//    {
//        prob*=pow(getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(newPos+contactShell_4[n]))),4.0);
//        prob_divide/=pow(getProbabilityFactor(monoType,uint32_t(ingredients.getLatticeEntry(oldPos+contactShell_4[n]))),4.0);
//    }
//
//
//
//    return prob/prob_divide;
//
// }
//
//
//
//
// template<class IngredientsType>
// void FeatureNNInteraction::fillLattice(IngredientsType& ingredients)
// {
//        const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
//        for(size_t n=0;n<molecules.size();n++)
//        {
//          ingredients.setLatticeEntry(molecules[n],uint8_t( molecules[n].getAttributeTag()));
//
//        }
// }
//
//
// template<class IngredientsType>
// void FeatureNNInteraction::applyMove(IngredientsType& ing, const MoveAddScMonomer& move)
// {
//     //get the position and attribute tag of the monomer to be inserted
//     VectorInt3 pos=move.getPosition();
//     int32_t type=move.getType();
//
//     //update lattice
//     ing.setLatticeEntry(pos,type);
//
// }
//////////////// end version with FeatureExcludedVolume/////////////////////////

/******************************************************************************
 *implementation of read and write classes
 ******************************************************************************/


#endif /*FEATURE_NN_INTERACTION_H*/
