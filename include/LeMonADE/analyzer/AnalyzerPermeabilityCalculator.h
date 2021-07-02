#ifndef AnalyzerPermeabilityCalculator_H
#define AnalyzerPermeabilityCalculator_H

#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>

/**
 * @file
 *
 * @class AnalyzerPermeabilityCalculator
 *
 * @brief Analyzer to calculate permeability across the bilayer which has the normal along z axis on LatticePowerOfTwo.
 *
 * @details This analyzer calculates the average permeability across the bilayer(normal along z axis) 
 * of objects (like copolymers and nano particles (see UpdaterLipidsCreator)) which are setup in as 
 * std::vector<MonomerGroup<molecules_type > > and solvent monomers. The analyzer calculates the solvent
 * and object permeability based on bufferDumpSize choosen by user and dumps it in file with default name
 * Permeability.dat.<br>
 * 
 * \b Important: <br>
 * (i)It is reccomended to create solvent and objects monomer with tags greater than and equal to 5(see UpdaterLipidsCreator).
 * If one needs to change this number from 5 to some other tag, it can be done by changing \SOLVENTAG below.
 * Thus, to use this analyzer, the general scheme for tags should be \b { Lipid monomers(heads and tails)< SOLVENTAG \b }
 * and \b { Non-Lipids monomers including Solvent >= SOLVENTAG \b }.<br>
 * 
 * (ii)In the case of the bilayers with pores, it is essential to find the position of pores are in xy plane.
 * The calculation of permeability is difficult in the region of pore, as there is no interface in this region. 
 * To solve this problem, this analyzer uses method similar to S. Podogin et. al. ACS Nano, 6:10555–10561, 2012.<br>
 * The bilayer plane is divided into bins and permeability is calculated in seperately in these bins. If a bin
 * has less than 30% of average lipid monomer in other bins, it is considered to be a pore state. This threshold 
 * of 30% can be controlled by \b PORE_THRESHOLD below.<br>
 * 
 * (iii) It is reccomended to simulate the system for 100 MCS(simInterval_=100)(or track permeable objects after 100 MCS)
 * before executing this analyzer to get proper results.<br>
 * 
 * @tparam IngredientsType
 *
 * @param ingredients_ the system holding the simulation box with bilayer, solvent and objects.
 * @param groups_ subgroup of molecules_type having groups of objects which translocate through bilayer.
 * 
 * @param binSize_ size of bin in which the bilayer plane is divided and permeability
 * is calculated in each of the bins.
 * @param factorSigma_ this is the factor which is multiplied to spread of lipid monomers along z-axis(standard deviation).
 * This was choosen to be 3 in S. Podogin et. al. ACS Nano, 6:10555–10561, 2012 and 5 in Werner et. al., Soft Matter 
 * 8:11714,2012. Please see these articles for further information.
 * @param outputFile_ the user defined output file name for dumping permeability values on disk.
 * @param simInterval_ the interval of MCS for which system is simulated before the analyzer is executed. It basically
 * tells interval(MCS) after which objects and solvent molecules are tracked. 
 * @param bufferDumpSize_ the buffer size after which to print the permeability values into the file. This value 
 * also automatically sets interval of MCS after which permeability is calculated. For example, for default permeability
 * value is calculated after every simInterval_*bufferDumpSize_=10000 MCS.
 * @param midplaneUpdaterInterval_ the interval of MCS after which midplane values should be updated. For example, for default 
 * values, it would be updated every simInterval_*midplaneUpdaterInterval_=15000MCS. 
 **/

#define SOLVENTAG 5
#define PORE_THRESHOLD 0.3

using namespace std;
template<class IngredientsType>
class AnalyzerPermeabilityCalculator{

protected:
    //!Typedefs for various underlying container holding the monomers
    //!and monomer subgroups.  
    typedef typename IngredientsType::molecules_type molecules_type;
    typedef std::vector< MonomerGroup< molecules_type> > group_type;
    typedef MonomerGroup< molecules_type> group_type_ind;
    const IngredientsType& ingredients;
    const molecules_type& molecules;
    const group_type& groups;
    
    //!More explainantion of these variables could be found below.    
    int32_t binSize,factorSigma,bufferDumpSize,simInterval,midplaneUpdaterInterval;
    std::vector<VectorDouble3> groupsCOM;
    int32_t boxXm_1,boxYm_1,boxZm_1,binBoxX,binBoxY,boxZ;

    //!Maps for storing indices of particle which are inside the boundaries(close to the bilayer).
    std::map<int, int> particlesInsideBoundaries;
    std::map<int, int> particlesInsideBoundariesOld;
    
    //!Vector for storing the indices of monomers which are solvent.
    std::vector < int32_t > solventIndices;

    //!Vectors and variables used by midplaneUpdater() function.
    std::vector<std::vector < float > >  midplane;
    std::vector<std::vector < float > >  counterMidplane;
    std::vector<std::vector < float > >  poreFlag;
    std::vector<VectorInt3> neigbours;
    int32_t counterSolvent,counterObjects,counterExecute;
    float sigma;
    
    //!Variables used by dumpPermeabilityPerMcs() function. 
    std::string outputFile;
    bool isFirstFileDump;
    
public:
    AnalyzerPermeabilityCalculator(const IngredientsType& ingredients_,const group_type& groups_,int32_t binSize_=8, int32_t factorSigma_=3,std::string outputFile_="Permeability.dat", int32_t bufferDumpSize_=100,int32_t midplaneUpdaterInterval_=150,int32_t simInterval_=100);
    
    //!General updater functions
    void execute();
    void initialize(){};
    void cleanup(){}; 
    
    //!Function for updating permeability values for objects and solvent monomers.     
    void permeabilityUpdate();    
    
    //!Function for updating position of midplane and width of bilayer.
    void midplaneUpdater();
    
    //!Function used by midplaneUpdater() to find the nearest neigbouring bin's midplane value.
    void findNearestMidplaneValue(int32_t x, int32_t y);
    
    //!Function to initiate the array with initial maps of particles inside/outside boundaries.
    void initParticleArrays();
    
    //!Function to find distances using minimum image convention.    
    int32_t reduceDistanceInPeriodicSpace(int32_t distance, int32_t period);
    
    //!Function to find the center of mass of subgroup of monomers(generally objects).    
    VectorDouble3 centerOfMass(const group_type_ind&  m);
    
    //!Function for dumping the permeability values in the file.        
    void dumpPermeabilityPerMcs();

};

/**
* @brief Constructor
**/
template<class IngredientsType>
AnalyzerPermeabilityCalculator<IngredientsType>::AnalyzerPermeabilityCalculator(const IngredientsType& ingredients_,const group_type& groups_, int32_t binSize_, int32_t factorSigma_,std::string outputFile_ ,int32_t bufferDumpSize_, int32_t midplaneUpdaterInterval_, int32_t simInterval_)
:ingredients(ingredients_)
,molecules(ingredients_.getMolecules())
,groups(groups_)
,binSize(binSize_)
,factorSigma(factorSigma_)
,outputFile(outputFile_)
,bufferDumpSize(bufferDumpSize_)
,midplaneUpdaterInterval(midplaneUpdaterInterval_)
,simInterval(simInterval_)
{   
    
    //The list of neigbours on 2-dimensional surface used by 
    //function findNearestMidplaneValue() to find a neigbour 
    //of a bin which doesn't have a pore.
    

    neigbours.push_back(VectorInt3(-1,0,0));
    neigbours.push_back(VectorInt3(0,-1,0));
    neigbours.push_back(VectorInt3(1,0,0));
    neigbours.push_back(VectorInt3(0,1,0));
    neigbours.push_back(VectorInt3(-1,1,0));
    neigbours.push_back(VectorInt3(1,-1,0));
    neigbours.push_back(VectorInt3(1,1,0));
    neigbours.push_back(VectorInt3(-1,-1,0));
    
    //The variables for number of bins in x and y direction
    //depending on periodic box length in those directions 
    //respectively.
   
    binBoxX=ingredients.getBoxX()/binSize;
    binBoxY=ingredients.getBoxY()/binSize;
    
    //Counters for various utilities below.
    
    counterExecute=0;
    counterObjects=0;
    counterSolvent=0;
    
    //Default values for box_length and boxXm_1 
    //for folding back and minimum image convention 
    //utilities.
    
    boxZ=ingredients.getBoxZ();
    
    boxXm_1=ingredients.getBoxX()-1;
    boxYm_1=ingredients.getBoxY()-1;
    boxZm_1=ingredients.getBoxZ()-1;
    
    //Initiate arrays of the size as number of bins in x and 
    //y direction.
    
    midplane.resize(binBoxX, std::vector<float>(binBoxY, boxZ/2.0));
    counterMidplane.resize(binBoxX, std::vector<float>(binBoxY, 0));
    poreFlag.resize(binBoxX, std::vector<float>(binBoxY, 0));
    sigma=5.0;
    isFirstFileDump=1;
    
    //Create an array of solvent indices so that it can be used to
    //directly access the position of solvent monomer using 
    //ingredients::molecules_type
    
    for(int32_t i=0;i<ingredients.getMolecules().size();i++)
        if(ingredients.getMolecules()[i].getAttributeTag()==SOLVENTAG)
            solventIndices.push_back(i);
        
   //groupsCOM vector stores the center of mass of objects(groups.size())
   //and solvent monomers.  
    
    groupsCOM.resize(groups.size()+solventIndices.size());
    
    //To start maps to store particles according to inside/outside 
    //of the boundaries.

    midplaneUpdater();
    initParticleArrays();
}

/**
* @details Execution of analyzer where it checks which particle are
* within the boundaries close to midplane updates the 
* particlesInsideBoundaries maps. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPermeabilityCalculator<IngredientsType>::execute(){
    
    //Increase the counter. 
    counterExecute++;

    int32_t x,y,z,distance,f2s,sign;
    
    //Loop to go over all the particles(objects and solvent)
    //to see if they are in the boundaries.
        
    for(size_t i=0;i<groups.size()+solventIndices.size();i++){
        
         //The objects boundaries can be set wider than
         //solvent boundaries.
    
        if(i<groups.size()){
            groupsCOM[i]=centerOfMass(groups[i]);
            f2s=factorSigma;
        }
        
        else{
            groupsCOM[i]=ingredients.getMolecules()[solventIndices[i-groups.size()]];
            f2s=3;
        }
        
         //Folding back coordinates for midplane array.
         //then, find the distance between midplane and particle 
         //using minimum image convention.
            
        x=(int32_t(groupsCOM[i].getX())&boxXm_1)/binSize;
        y=(int32_t(groupsCOM[i].getY())&boxYm_1)/binSize;
        
        distance=int32_t(groupsCOM[i].getZ()) - midplane[x][y];
        distance=reduceDistanceInPeriodicSpace(distance,boxZ);
        
        if(abs(distance)<=f2s*sigma){
             //check if particle is in upper side or
             //lower side.
            
            sign = (distance<0)?-1:1;
            particlesInsideBoundaries[i]=sign;
        }
    }
    
    //Update the permeability values. 
    permeabilityUpdate();
    
    
    //Dump values in the file if buffer size is reached.    
    if(counterExecute%bufferDumpSize==0)
        dumpPermeabilityPerMcs();
    
    //Update midplane values according to set up midplaneUpdaterInterval.    
    if(counterExecute%midplaneUpdaterInterval==0)
        midplaneUpdater();
    
    //Reset map values.    
    for(int32_t i=0;i<particlesInsideBoundariesOld.size();i++){
        particlesInsideBoundariesOld[i]=particlesInsideBoundaries[i];
        particlesInsideBoundaries[i]=0;
    }
    
}

/**
* @details Updates the permeability values by comparing maps 
* particlesInsideBoundariesOld and particlesInsideBoundaries. If a particle
* was inside before(Old) and exit boundary from a direction opposite to it entered in, 
* then it counts as one translocation event. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPermeabilityCalculator<IngredientsType>::permeabilityUpdate(){
     
    int32_t x,y,upDownIndicator;
    
     //Loop to go over all the particles(objects and solvent)
     //which were present in last execution and check if they 
     //still present, if not then, check if they went from 
     //opposite direction they came in and count a translocation 
     //event. If a particle is still present, then store the
     //direction it came inside from.
    
    for(int32_t i=0;i<particlesInsideBoundariesOld.size();i++){
        
        if(abs(particlesInsideBoundariesOld[i])==1&&particlesInsideBoundaries[i]==0){
                        
            x=(int32_t(groupsCOM[i].getX())&boxXm_1)/binSize;
            y=(int32_t(groupsCOM[i].getY())&boxYm_1)/binSize;
            
            upDownIndicator=int32_t(groupsCOM[i].getZ())-midplane[x][y];
            upDownIndicator=reduceDistanceInPeriodicSpace(upDownIndicator,boxZ);
            
            if(particlesInsideBoundariesOld[i]*upDownIndicator<0){

                if(i<groups.size())
                    counterObjects++;
                else
                    counterSolvent++;
            }
        }
        //Store(or remember) the direction of the particles from where it entered in the boundary.
        else if(abs(particlesInsideBoundariesOld[i])==1&&abs(particlesInsideBoundaries[i])==1)
            particlesInsideBoundaries[i]=particlesInsideBoundariesOld[i];
    }
};

/**
* @details Updates the midplane and width(sigma) value of bilayer. Since the
* bilayer is always translating in z axis. The boundaries position and size 
* need to updated regularly. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPermeabilityCalculator<IngredientsType>::midplaneUpdater(){

     //Reset all values to zero
     
     for(int32_t x1=0;x1<binBoxX;x1++)
         for(int32_t y1=0;y1<binBoxY;y1++){
             midplane[x1][y1]=0;
             counterMidplane[x1][y1]=0;
             poreFlag[x1][y1]=0;
         }
             
     int32_t x,y,z;
     int32_t referenceI=0;
     
     //If bilayer is near the boundaries of z axis 
     //of the box(z=0 or z=boxZ-1) even folded back values 
     //fail to give appropriate for the midplane values.
     //This analyzer uses a reference lipid monomer on 
     //the bilayer and use minimum image convention to find the  
     //midplane value close to this monomer.
     
     //Find the index of reference monomer
     for(int32_t i=0;i<ingredients.getMolecules().size();i++){
         
         int32_t type = ingredients.getMolecules()[i].getAttributeTag();
         
         if(type<SOLVENTAG) {referenceI=i; break;}
     }
     
     int32_t z_Ref=ingredients.getMolecules()[referenceI].getZ();
     
     //Find midplane values wrt this reference value.     
     for(int32_t i=0;i<ingredients.getMolecules().size();i++){
         
         int32_t type = ingredients.getMolecules()[i].getAttributeTag();
         
         if(type>=SOLVENTAG) continue;
         
         x=(int32_t(ingredients.getMolecules()[i].getX())&boxXm_1)/binSize;
         y=(int32_t(ingredients.getMolecules()[i].getY())&boxYm_1)/binSize;

         z=ingredients.getMolecules()[i].getZ();
         
         midplane[x][y] += (reduceDistanceInPeriodicSpace(z-z_Ref,boxZ));
         counterMidplane[x][y]++;
     }
     
    float maxValue=0;         
         
     for(int32_t x1=0;x1<binBoxX;x1++)
         for(int32_t y1=0;y1<binBoxY;y1++){
             
             midplane[x1][y1] /= counterMidplane[x1][y1];
             
             //Shift the values to original box coordinates.
             midplane[x1][y1] = (int32_t(midplane[x1][y1])+z_Ref)&boxZm_1;
             
             
             //The maximum value of lipids in a bin(when a pore is present)
             //is close to average number of lipids in bins without pores.
             maxValue = (maxValue>counterMidplane[x1][y1])?maxValue:counterMidplane[x1][y1];    
        }
        
     //Find pore flags.
     for(int32_t x1=0;x1<binBoxX;x1++)
         for(int32_t y1=0;y1<binBoxY;y1++)
             if(counterMidplane[x1][y1]<(PORE_THRESHOLD*maxValue)){poreFlag[x1][y1]=1;}
             
     //If pores are present set the midplane value 
     //of those bins to their neigbour's midplane values.             
     for(int32_t x1=0;x1<binBoxX;x1++)
         for(int32_t y1=0;y1<binBoxY;y1++)
             if(poreFlag[x1][y1]==1)
                 findNearestMidplaneValue(x1,y1);
             
    int32_t counterSigma=0;
    float distance;
    sigma=0;
    
    //Find the width of bilayer.
    for(int32_t i=0;i<ingredients.getMolecules().size();i++){
        int32_t type = ingredients.getMolecules()[i].getAttributeTag();
        
        if(type>=SOLVENTAG) continue;
        
        x=(int32_t(ingredients.getMolecules()[i].getX())&boxXm_1)/binSize;
        y=(int32_t(ingredients.getMolecules()[i].getY())&boxYm_1)/binSize;
        z=ingredients.getMolecules()[i].getZ();
        
        if(!poreFlag[x][y]){
            distance=(z-midplane[x][y]);
            distance=reduceDistanceInPeriodicSpace(distance,boxZ);
            sigma += distance*distance;
            counterSigma++;
        }
    }
    
    sigma=sqrt(sigma/counterSigma);
};

/**
* @details Finds a niegbouring bin without any pore and sets a midplane of bin with given
* coordinates to that neigbouring bin.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*
* @param x1 x coordinate of the bin with a pore.
* @param y1 y coordinate of the bin with a pore.
*/

template<class IngredientsType>
void AnalyzerPermeabilityCalculator<IngredientsType>::findNearestMidplaneValue(int32_t x1, int32_t y1){
    
    int32_t neighX,neighY=0;
    
    for(int nL=0;nL<neigbours.size();nL++){
        //Fold back neighbour list in case pore is at an edge of
        //the box.
        neighX=((neigbours[nL][0]+x1)%binBoxX+binBoxX)%binBoxX;
        neighY=((neigbours[nL][1]+y1)%binBoxY+binBoxY)%binBoxY;

        if(poreFlag[neighX][neighY]==0){
            midplane[x1][y1]=midplane[neighX][neighY];
            break;
         }
    }    
};

/**
* @details Initializes the maps particlesInsideBoundaries and particlesInsideBoundariesOld.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
*/

template<class IngredientsType>
void AnalyzerPermeabilityCalculator<IngredientsType>::initParticleArrays(){
    
    int32_t x,y,distance,sign,f2s;
    
    for(size_t i=0;i<groups.size()+solventIndices.size();i++){
        
        if(i<groups.size()){
            groupsCOM[i]=centerOfMass(groups[i]);
            f2s=factorSigma;
        }
        
        else{
            groupsCOM[i]=ingredients.getMolecules()[solventIndices[i-groups.size()]];
            f2s=3;
        }
    
        x=(int32_t(groupsCOM[i].getX())&boxXm_1)/binSize;
        y=(int32_t(groupsCOM[i].getY())&boxYm_1)/binSize;
        
        distance=(int32_t(groupsCOM[i].getZ())&boxZm_1) - midplane[x][y];
        distance=reduceDistanceInPeriodicSpace(distance,boxZ);
        
        if(abs(distance)<=f2s*sigma){
            sign = (distance<0)?-1:1;
            particlesInsideBoundariesOld.insert(make_pair(i ,sign));
        }
        else
            particlesInsideBoundariesOld.insert(make_pair(i,0));
        
        particlesInsideBoundaries.insert(make_pair(i,0));
 }
 
};

/**
* @details Prints the permeability per mcs into files.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void  AnalyzerPermeabilityCalculator<IngredientsType>::dumpPermeabilityPerMcs(){
    
    std::vector< std::vector< double > > resultsTimeseries(3,std::vector< double >(1));
    resultsTimeseries[0][0] = ingredients.getMolecules().getAge();
    resultsTimeseries[1][0] = counterObjects/(bufferDumpSize*simInterval);
    resultsTimeseries[2][0] = counterSolvent/(bufferDumpSize*simInterval);

    //if it is written for the first time, include comment in the output file
    if(isFirstFileDump){
    
            std::stringstream commentTimeSeries;
            commentTimeSeries<<"Created by AnalyzerPermeabilityCalculator\n";
            commentTimeSeries<<"file contains time series of average permeability averaged per \n";
            commentTimeSeries<<"format: mcs\t P_objects\t P_solvent\n";

            ResultFormattingTools::writeResultFile(
                    outputFile,
                    ingredients,
                    resultsTimeseries,
                    commentTimeSeries.str()
            );

            isFirstFileDump=false;
    }
    //otherwise just append the new data
    else{
            ResultFormattingTools::appendToResultFile(outputFile,
                                                        resultsTimeseries);
    }
    
    counterObjects=0;
    counterSolvent=0;

 }
 
/**
* @details Calculates the center of mass of objects given a subgroup of 
* monomers of the object.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template < class IngredientsType > 
VectorDouble3 AnalyzerPermeabilityCalculator<IngredientsType>::centerOfMass(const group_type_ind& m)
{
    VectorDouble3 CoM_sum;
    
    for ( uint i = 0; i < m.size(); ++i)
    {
      CoM_sum.setX( CoM_sum.getX() + m[i].getX() );
      CoM_sum.setY( CoM_sum.getY() + m[i].getY() );
      CoM_sum.setZ( CoM_sum.getZ() + m[i].getZ() );
    }
    
    float inv_N = 1.0 / double ( m.size() );
    
    VectorDouble3 CoM (
	float ( CoM_sum.getX() ) * inv_N, 
	float ( CoM_sum.getY() ) * inv_N,
	float ( CoM_sum.getZ() ) * inv_N);
      
    return CoM;
 
};

/**
* @details Finds the distances in the minimum image convention give a period. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template <class IngredientsType> 
int32_t AnalyzerPermeabilityCalculator<IngredientsType>::reduceDistanceInPeriodicSpace(int32_t distance, int32_t period)
{
    if(distance > period/2)
    {while(-(distance-period)<period/2) distance -= (period); return distance;}
    
    else if(-distance > period/2) 
    {while((distance+period)<period/2) distance += (period); return distance;}
    
    else
    return distance;
}
  
#endif
