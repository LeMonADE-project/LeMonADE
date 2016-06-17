#include <LeMonADE/feature/FeatureNNInteraction.h>

/*******************************************************************************
 *constructor: initializes the lookup arrays
 ******************************************************************************/
FeatureNNInteraction::FeatureNNInteraction()
{
    for(size_t n=0;n<11;n++)
    {
        for(size_t m=0;m<11;m++)
        {
            interactionTable[m][n]=0.0;
            probabilityLookup[m][n]=1.0;
        }
    }
}

/*******************************************************************************
 *destructor
 ******************************************************************************/
FeatureNNInteraction::~FeatureNNInteraction(){}


void FeatureNNInteraction::setInteraction(uint32_t typeA,uint32_t typeB,double energy)
{
    if(0<typeA && typeA<=10 && 0<typeB && typeB<=10)
    {
        interactionTable[typeA][typeB]=energy;
        interactionTable[typeB][typeA]=energy;
        probabilityLookup[typeA][typeB]=exp(-energy);
        probabilityLookup[typeB][typeA]=exp(-energy);
        std::cout<<"set interation between types "<<typeA<<" and "<<typeB<<" to "<<energy<<"kT\n";
    }
    else
    {
        throw std::runtime_error("******FeatureNNInteraction::addInteraction(). This feature supports only attribute tags t with 0<t<=10!\n");
    }
}

//returns the interaction between two types of monomers in units of interactionBase
//int32_t FeatureNNInteraction::getInteraction(uint32_t typeA,uint32_t typeB) const
double FeatureNNInteraction::getInteraction(uint32_t typeA,uint32_t typeB) const
{

    if(0<typeA && typeA<=10 && 0<typeB && typeB<=10)
        return interactionTable[typeA][typeB];

    else
        throw std::runtime_error("***FeatureNNInteraction::getInteraction()...trying to get undefined interaction constant***\n");

}


/////////// for version with FeatureExcludedVolume /////////////////////////////
//const VectorInt3 FeatureNNInteraction::contactShell_1[24]=
//{
//   VectorInt3(2,1,1),
//   VectorInt3(-2,1,1),
//   VectorInt3(2,-1,1),
//   VectorInt3(2,1,-1),
//   VectorInt3(-2,-1,1),
//   VectorInt3(-2,1,-1),
//   VectorInt3(2,-1,-1),
//   VectorInt3(-2,-1,-1),
//   VectorInt3(1,2,1),
//   VectorInt3(-1,2,1),
//   VectorInt3(1,-2,1),
//   VectorInt3(1,2,-1),
//   VectorInt3(-1,-2,1),
//   VectorInt3(-1,2,-1),
//   VectorInt3(1,-2,-1),
//   VectorInt3(-1,-2,-1),
//   VectorInt3(1,1,2),
//   VectorInt3(-1,1,2),
//   VectorInt3(1,-1,2),
//   VectorInt3(1,1,-2),
//   VectorInt3(-1,-1,2),
//   VectorInt3(-1,1,-2),
//   VectorInt3(1,-1,-2),
//   VectorInt3(-1,-1,-2)
//};
//
//const VectorInt3 FeatureNNInteraction::contactShell_2[24]=
//{
//   VectorInt3(2,1,0),
//   VectorInt3(-2,1,0),
//   VectorInt3(-2,-1,0),
//   VectorInt3(2,-1,0),
//   VectorInt3(0,2,1),
//   VectorInt3(0,-2,1),
//   VectorInt3(0,-2,-1),
//   VectorInt3(0,2,-1),
//   VectorInt3(1,0,2),
//   VectorInt3(1,0,-2),
//   VectorInt3(-1,0,-2),
//   VectorInt3(-1,0,2),
//   VectorInt3(2,0,1),
//   VectorInt3(-2,0,1),
//   VectorInt3(2,0,-1),
//   VectorInt3(-2,0,-1),
//   VectorInt3(1,2,0),
//   VectorInt3(-1,2,0),
//   VectorInt3(1,-2,0),
//   VectorInt3(-1,-2,0),
//   VectorInt3(0,1,2),
//   VectorInt3(0,-1,2),
//   VectorInt3(0,1,-2),
//   VectorInt3(0,-1,-2)
//};
//
//const VectorInt3 FeatureNNInteraction::contactShell_4[6]=
//{
//   VectorInt3(2,0,0),
//   VectorInt3(-2,0,0),
//   VectorInt3(0,2,0),
//   VectorInt3(0,-2,0),
//   VectorInt3(0,0,2),
//   VectorInt3(0,0,-2)
//};
/////////// end version with FeatureExcludedVolume /////////////////////////////
