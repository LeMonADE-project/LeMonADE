
template<class IngredientsType>
class TesterFunctions
{
public:
    //!Function to get total number of monomers in system with a given attribute.     
    int32_t getNumMonWithAttribute(IngredientsType& ingredients, int32_t attribute ){
        
        int32_t numMon=0;
        for(int32_t i=0;i<ingredients.getMolecules().size();i++){
            if(ingredients.getMolecules()[i].getAttributeTag()==attribute) numMon++;
        }
        return numMon;
    }
    
    //!Function to get total number of monomers in system with a given number of links.     
    int32_t getNumMonWithLinks(IngredientsType& ingredients, int32_t links ){
        
        int32_t numMon=0;
        for(int32_t i=0;i<ingredients.getMolecules().size();i++){
            if(ingredients.getMolecules().getNumLinks(i)==links) numMon++;
        }
        return numMon;
    }
    
    //!Function to get total number of non solvent monomers.     
    int32_t getTotalNonsolventMonomers(IngredientsType& ingredients){
        
        int32_t NS_monomers=0; 
        for(int32_t i=0;i<ingredients.getMolecules().size();i++)
            if(ingredients.getMolecules()[i].getAttributeTag()!=5)
                NS_monomers++;
        
        return NS_monomers;
    };
    
    //!Function to get total number of  solvent monomers.     
    int32_t getTotalsolventMonomers(IngredientsType& ingredients){
        
        int32_t S_monomers=0; 
        for(int32_t i=0;i<ingredients.getMolecules().size();i++)
            if(ingredients.getMolecules()[i].getAttributeTag()==5)
                S_monomers++;
        
        return S_monomers;
    };

    //!Function to get total number of monomers in system with a given attribute.     
    std::vector<VectorFloat3> getOrientations(IngredientsType& ingredients, int32_t tail_length, bool ringHead=false){
        
        std::vector<VectorFloat3> Orientations;
        int headlinks;
        if(ringHead)
            headlinks=4;
        else
            headlinks=3;
        
        for(int32_t i=0;i<ingredients.getMolecules().size();i++){

            if((ingredients.getMolecules().getNumLinks(i)==headlinks)&&ingredients.getMolecules().getNumLinks(i+tail_length)==1){
                
                VectorInt3 headpos=ingredients.getMolecules()[i];
            
                VectorInt3 tailend1pos=ingredients.getMolecules()[i+tail_length];
                VectorInt3 tailend2pos=ingredients.getMolecules()[i+2*tail_length];
                
                VectorFloat3 avgEndspos =tailend1pos+tailend2pos;
                avgEndspos /= 2;
                
                Orientations.push_back(headpos-avgEndspos);
        }
       }
        return Orientations;
    };
    
};
