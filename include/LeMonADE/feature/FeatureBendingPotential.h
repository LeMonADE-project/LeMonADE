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


#ifndef FEATURE_BENDING_POTENTIAL_H
#define FEATURE_BENDING_POTENTIAL_H

/**
 * @file
 * @date 2019/08/11
 * @author Ankush Checkervarty
 * @brief Definition and implementation of class template FeatureBendingPotential
**/

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureBendingPotentialReadWrite.h>



namespace BendingPotentials
{
    /*This namespace is created so different functional form could be created or used 
     *for kind of bending potential which is supposed to be implemented. Here for example,
     *a simple harmonic potential has been shown.
     * Functional form: <br>
     * 
     * U(\theta)=\theta**2.
     * 
     */
    
    float simpleHarmonic(VectorInt3& bondVec1, VectorInt3& bondVec2){
        
        float magnitude,cosine;
        
        magnitude =bondVec1*bondVec1;
        magnitude *= bondVec2*bondVec2;
        magnitude =sqrt(magnitude);
        
        cosine= (bondVec1*bondVec2)/magnitude;

        return acos(cosine)*acos(cosine);
    }
}


/**
 * @class FeatureBendingPotential
 * @brief Provides implementation of bending potential on the BFM polymer chains.
 *
 * @details
 * Different potential strength can be set for a bonded monomer type by using
 * setBendingPotential(type,energy). This feature only works for a linear chain
 * of a monomer type, thus, the monomer to be moved could be attached to either 
 * one or two of the same monomer type. 
 * The linear chain could bonded to chains other monomer type. For example, 
 * in a lipid(see UpdaterLipidsCreator). Even though tail monomers are of different 
 * type than head monomer and head monomers are bonded to tail monomers. One can 
 * set different (or None at default) bending potential strength for head and tail monomers.  
 * \b IMPORTANT:<br>
 * (i) It is important to create solvent(or non bonded monomers) after the bonded objects.
 * This reduces memory requirements. This feature throws runtime_error if opposite is done.<br>
 * (ii)The type can be any integer between 1 and 10.<br>
 * The feature adds the bfm-file command !bending_potential A energy
 * for monomers of type A with bending potential of E in kT
 **/

class FeatureBendingPotential:public Feature
{

private:
    
    //! total number of non solvent monomers in the system    
    int32_t numNonSolvent;
    
    //! bending potential strength of the type of the monomers. type range [1:10]
    double bpStrengthTable[11];

    //! Lookup table for exp(-bpStrengthTable[a])
    double probabilityLookup[11][180][180];

    //! Lookup to tag chains ends and their neigbours.
    std::vector <int32_t> chainsEnds;

    //! Returns this feature's factor for the acceptance probability for the given Monte Carlo move and type.
    template<class IngredientsType> 
    double calculateAcceptanceProbability(const IngredientsType& ingredients,
                        const MoveLocalSc& move,
                        int32_t monoType) const;
                                            
    //! Used to fill in the values in chainsEnds.
    void tagNiegbsandEnds(int32_t index, int32_t sameAttNeigbIndex);                                     

public:

    FeatureBendingPotential();
    ~FeatureBendingPotential(){}


    //This feature adds interaction energies, so it requires FeatureBoltzmann
    typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

    //FeatureAttributes needs to be in front to decide which type will 
    //have which bending potential strength.
    //FeatureBondset needs to be in front. As it required to fill probabilityLookup table.
    typedef LOKI_TYPELIST_2(
        FeatureAttributes<int32_t>,
        FeatureBondset<>)
        required_features_front;


    //! check for all Monte Carlo moves without special check functions (always true)
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;

    //! check for standard sc-BFM local move
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;
    
    //! This is simple function to convert bond vector to integer ID in one to one mapping. 
    uint32_t bondVectorToIndex(const VectorInt3& bondVector) const;                                     

    //!adds bending potential strength to type.
    template<class IngredientsType>
    void setBendingPotential(IngredientsType& ingredients, int32_t type,double energy);

    //!sets up the chain ends tags in chainsEnds.
    template<class IngredientsType>
    void setChainEnds(IngredientsType& ingredients);

    //!returns the bending potential strength for type.
    double getBendingPotential(int32_t type) const;

    //!export bfm-file read command !bending_potential
    template<class IngredientsType>
    void exportRead(FileImport <IngredientsType>& fileReader);

    //!export bfm-file write command !bending_potential
    template<class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;

};


//////////////////  IMPLEMENTATION OF MEMBERS //////////////////////////////////

/**
 * @brief Constructor
 **/

FeatureBendingPotential::FeatureBendingPotential()
{
  //initialize the bpStrengthTable and probability lookups with default values
  for(size_t n=0;n<11;n++)
    {
      bpStrengthTable[n]=0;
      for(size_t m=0;m<180;m++)
          for(size_t o=0;o<180;o++)
              probabilityLookup[n][m][o]=1.0;
    }
}

/**
 * @details calculates the factor for the acceptance probability of the move
 * arising from the bending potential and adds it to the move.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalSc
 * @return true (always)
 **/

template<class IngredientsType>
bool FeatureBendingPotential::checkMove(const IngredientsType& ingredients,
							 MoveLocalSc& move) const
{
  //add the probability factor coming from this feature, then return true,
  //because the total probability is evaluated by FeatureBoltzmann at the end.
    
  //if monoType has zero strength then returns true without making any change in 
  //probability
  int32_t monoType=ingredients.getMolecules()[move.getIndex()].getAttributeTag();
  if(!bpStrengthTable[monoType]) return true;
  
  double prob=calculateAcceptanceProbability(ingredients,move,monoType);
  move.multiplyProbability(prob);
  return true;
}


/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * - !bending_potential
 * .
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileReader File importer for the bfm-file
 **/
template<class IngredientsType>
void FeatureBendingPotential::exportRead(FileImport< IngredientsType >& fileReader)
{
//   typedef FeatureBendingPotential my_type;
  fileReader.registerRead("!bending_potential",new ReadBendingPotential<IngredientsType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * - !bending_potential
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureBendingPotential::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
//   typedef FeatureBendingPotential my_type;
  fileWriter.registerWrite("!bending_potential",new WriteBendingPotential<IngredientsType>(fileWriter.getIngredients_()));
}


/**
 * @details The function calculates the factor for the acceptance probability
 * for the local move given as argument. The calculation is based on whether the 
 * choosen monomers is at middle of the chain, last monomer away from start of 
 * the chain monomer, one monomer away from end of the chain monomer or any end monomers.
 * 
 * @tparam IngredientsType The type of the system including all features.
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move reference to the local move for which the calculation is performed.
 * @param [in] monoType attribute of the monomer moved. 
 * 
 * @return acceptance probability factor for the move arising from bending potential interactions.
 **/

template<class IngredientsType>
double FeatureBendingPotential::calculateAcceptanceProbability(
    const IngredientsType& ingredients,
    const MoveLocalSc& move,
    int32_t monoType) const
{

    int32_t index=move.getIndex();
    
    VectorInt3 presentPos=ingredients.getMolecules()[move.getIndex()];
    VectorInt3 direction=move.getDir();
    VectorInt3 bondvector1,bondvector2;
    
    double prob=1.0;

    VectorInt3 futurePos = presentPos + direction;
    /**distance from the start of the chain and end of the chain*/
    int32_t distfromSOC=3&chainsEnds[index];
    int32_t distfromEOC=(3<<2)&chainsEnds[index];
    distfromEOC = distfromEOC>>2;
    
    float probBeforeMove,probAfterMove;    
    
    /**Generally, a move effects all the bending angles of two monomer in postive
     * direction and two monomer negative direction to moved monomer. However, special
     * care has to be taken if the monomer is end(or near) of the chain or start of chain.*/

    //check if it is the end monomer.
    if(distfromSOC&&distfromEOC){
        VectorInt3 presentPosm_1=ingredients.getMolecules()[index-1];
        VectorInt3 presentPosp_1=ingredients.getMolecules()[index+1];
        
        bondvector1=presentPosp_1-futurePos;
        bondvector2=futurePos-presentPosm_1;
        
        probAfterMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];
        
        bondvector1=presentPosp_1-presentPos;
        bondvector2=presentPos-presentPosm_1;

        probBeforeMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];
            
        prob *= probAfterMove/probBeforeMove;
    }

  //check end monomer lies in negative direction.
   if(distfromSOC>1){
       
       VectorInt3 presentPosm_1=ingredients.getMolecules()[index-1];
       VectorInt3 presentPosm_2=ingredients.getMolecules()[index-2];
        
       bondvector1=futurePos-presentPosm_1;
       bondvector2=presentPosm_1-presentPosm_2;
       
       probAfterMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];

       bondvector1=presentPos-presentPosm_1;
        
       probBeforeMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];
       
       prob *= probAfterMove/probBeforeMove;
}  
     
  //check end monomer lies in positive direction.
   if(distfromEOC>1){
       
       VectorInt3 presentPosp_1=ingredients.getMolecules()[index+1];
       VectorInt3 presentPosp_2=ingredients.getMolecules()[index+2];
       
       bondvector1=presentPosp_2-presentPosp_1;
       bondvector2=presentPosp_1-futurePos;

       probAfterMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];
       
       bondvector2=presentPosp_1-presentPos;

       probBeforeMove=probabilityLookup[monoType][bondVectorToIndex(bondvector1)][bondVectorToIndex(bondvector2)];
       
       prob *= probAfterMove/probBeforeMove;
  }
  
  return prob;

}
/**
 * @details This function fills in tags in chainsEnds. Using bit shift operations it fills
 * first 2 bits with distance from start of the chain and second 2 bits distance from end 
 * of the chain.
 * 
 * @tparam IngredientsType The type of the system including all features.
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * 
 **/

template<class IngredientsType>
void FeatureBendingPotential::setChainEnds(IngredientsType& ingredients){
  //First find out all the bonded monomers to initiate the size of chain ends. 
  //Throw error if solvent monomers are created before objects.  
  numNonSolvent=0;
  bool isSolventFirst=true;
  for(size_t i=0;i<ingredients.getMolecules().size();i++)
      if(ingredients.getMolecules().getNumLinks(i)){
          numNonSolvent++;
          isSolventFirst=false;}
      else if(isSolventFirst)
          throw std::runtime_error("FeatureBendingPotential::Solvent should not created before the bonded objects.");

  //Initiate all the chain ends with all distance from the 
  //start or end of the chain to be 2. As any number >1 would act same.        
  chainsEnds.resize(numNonSolvent,2+(2<<2));
  
  for(size_t i=0;i<ingredients.getMolecules().size();i++){
      
      int32_t numLinks = ingredients.getMolecules().getNumLinks(i);
      //If it's the chain end, then tag ends and neigbours.
      if(numLinks==1)
          tagNiegbsandEnds(i,ingredients.getMolecules().getNeighborIdx(i, 0));
      //If it's not the end, find if which bonded neigbours is the same as monomer 
      //it self and then tag neigbours.      
      else if(numLinks>1){
          
          int32_t numSameConnection=0;
          int32_t sameAttNeigbIndex=0;
          
          for(int32_t neigb=0;neigb<numLinks;neigb++){
              
              int32_t neigbIdx = ingredients.getMolecules().getNeighborIdx(i, neigb);
              
              int32_t attributeSelf = ingredients.getMolecules()[i].getAttributeTag();          
              int32_t attributeNeigb = ingredients.getMolecules()[neigbIdx].getAttributeTag();
              
              if(attributeSelf==attributeNeigb){
                  sameAttNeigbIndex=neigbIdx;
                  numSameConnection++;}
          }
          
          if(numSameConnection==1)
              tagNiegbsandEnds(i,sameAttNeigbIndex);
      }
  }
  
}  
/**
 * @details Used by setChainEnds to tag end and bonded neigbour distance from the closest 
 * chain end.
 * 
 * @param [in] index index of the chain end monomer to be tagged.
 * @param [in] sameAttNeigbIndex index of bonded neigbours with the same attribute.
 * 
 **/

void FeatureBendingPotential::tagNiegbsandEnds(int32_t index, int32_t sameAttNeigbIndex){
    //check if the it's end of the chain or
    //start of the chain.
    if(sameAttNeigbIndex>index){
        int32_t tempoOldVal = chainsEnds[index]&(3<<2);
        chainsEnds[index] = (0&3);
        chainsEnds[index] += tempoOldVal;

        tempoOldVal = chainsEnds[index+1]&(3<<2);
        chainsEnds[index+1] = (1&3);
        chainsEnds[index+1] += tempoOldVal;
    }
    else{
        int32_t tempoOldVal = chainsEnds[index]&3;
        chainsEnds[index] = (0&3)<<2;
        chainsEnds[index] += tempoOldVal;

        tempoOldVal = chainsEnds[index-1]&3;
        chainsEnds[index-1] = (1&3)<<2;
        chainsEnds[index-1] += tempoOldVal;
    }
}
/**
 * @param type monomer attribute tag in range [1,10].
 * @param energy magnitude of bending potential energy.
 * @throw std::runtime_error In case type exceed range [1,255].
 **/
template<class IngredientsType>
void FeatureBendingPotential::setBendingPotential(IngredientsType& ingredients,
                                                  int32_t type,
                                                  double energy)
{
    VectorInt3 bondVec1,bondVec2;
    if(0<type && type<=10)
      {
        bpStrengthTable[type]=energy;
        
        std::map <int32_t,VectorInt3>::const_iterator it,it2;
        //go over all the bondvector pair possible
        for(it=ingredients.getBondset().begin();it!=ingredients.getBondset().end();it++)
            for(it2=ingredients.getBondset().begin();it2!=ingredients.getBondset().end();it2++){
                bondVec1=it->second;
                bondVec2=it2->second;
                //find the value from bending potential form
                float potentialVal=BendingPotentials::simpleHarmonic(bondVec1,bondVec2);
                probabilityLookup[type][bondVectorToIndex(bondVec1)][bondVectorToIndex(bondVec2)]=exp(-potentialVal*energy);
            }
            
        std::cout<<"set bending potential for types ";
        std::cout<<type<<" to "<<energy<<"kT\n";
      }
    else
      {
          std::stringstream errormessage;
          errormessage<<"FeatureBendingPotential::setBendingPotential(type,energy).\n";
          errormessage<<"type "<<type<<": Type out of range\n";
          throw std::runtime_error(errormessage.str());
      }
}

/**
 * @param type monomer attribute tag in range [1,10]
 * @throw std::runtime_error In case type exceed range [1,10]
 * @return bending potential energy for type monomer
 **/
double FeatureBendingPotential::getBendingPotential(int32_t type) const
{

    if(0<type && type<=10)
        return bpStrengthTable[type];
    else
    {
      std::stringstream errormessage;
      errormessage<<"FeatureBendingPotential::getBendingPotential(type).\n";
      errormessage<<"type "<<type<<": Type out of range\n";
      throw std::runtime_error(errormessage.str());
    }

}
/**
 * @details Translates a bond-vector into the corresponding lookup table index.
 * The method used here is based on brute force-like way to find out which linear 
 * combination of bond vector components to create a unique index. This method saves
 * memory compare to simple bitshift operation based conversion where 512 size array is 
 * required, here only 180 size array is required.
 * 
 * @param bondVector Reference to VectorInt3 as bond-vector to put on fast look-up-table.
 * @return A unique number representing the bond-vector.
 */

inline uint32_t FeatureBendingPotential::bondVectorToIndex(const VectorInt3& bondVector) const
{
	return 17* bondVector.getX() + 12*bondVector.getY() - 20*bondVector.getZ() + 86;
}


#endif 
