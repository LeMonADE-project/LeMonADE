#ifndef PAIRDISTRIBUTION_H
#define PAIRDISTRIBUTION_H

#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>

using namespace std;

/**
 * @file AnalyzerPoreFinder.h
 *
 * @class VectorInt2
 *
 * @brief This class is created as substitute to VectorInt3 in in class Vector3D.h, as the 
 * most calculations required in class AnalyzerPoreFinder only requires a two dimensional space.
 **/


class VectorInt2{
public:
    VectorInt2(int32_t x_=0,int32_t y_=0):
    x(x_),
    y(y_){}
    
    void setX(int32_t x_){x=x_;};
    void setY(int32_t y_){y=y_;};
    int32_t getX(){return x;};
    int32_t getY(){return y;};
    
    float getLength(){return sqrt(x*x+y*y);};

    void setAllCoordinates(int32_t x_,int32_t y_){x=x_; y=y_;};
    
    VectorInt2 operator+(const VectorInt2& b) {
         VectorInt2 vec;
         vec.x = this->x + b.x;
         vec.y = this->y + b.y;
         return vec;
      };
      
    VectorInt2 operator/(int32_t divisor) {
         VectorInt2 vec;
         vec.x = round(this->x, divisor);
         vec.y = round(this->y, divisor);
         return vec;
      };
      
   friend ostream& operator<<(ostream& out, const VectorInt2& b){
         out<<b.x<<"\t"<<b.y;
         return out;
      };

    int32_t operator[](int32_t i) {
        if(i==0)
            return x;
        else if(i==1)
            return y;
        else
            throw std::runtime_error("Error : value of range");
      };
      
    int32_t round(int32_t i1,int32_t i2){
        float value=float(i1)/i2;
        return(floor(value+0.5));};

private:
    int32_t x,y;
    
};

/**
 * @file AnalyzerPoreFinder.h
 *
 * @class AnalyzerPoreFinder 
 *
 * @brief  Analyzer to find the coordinates of pore on a bilayer, using heirarchial cluster analysis.
 * It also calculates the position of pore points and copolymers center of mass with respect to pore centriod.
 * 
 * @details This analyzer finds the cluster of points in set of all the lattice points(in a periodic 
 * box) on x-y plane which represent the coordinates of pore. This is done by considering all point on
 * x-y plane which are free of tail monomers all along z axis as pore points. The cluster analysis then, 
 * finds out all the cluster of points which are connected with the cluster size more than
 * CLUSTER_SIZE_THRESHOLD. This threshold can be changed with the macro below.
 * Furthermore, in considering the copolymers it is essential to count only the polymer whose COM distance from the 
 * bilayer midplane is below a threshold which can set by POLYMER_DIST_THRESHOLD.<br>
 * 
 * Also important to set all lipid monomers below SOLVENTAG and all non lipid monomers including solvent to SOLVENTAG.<br>
 * Please see AnalyzerPermeabilityCalculator and UpdaterLipidsCreator for further details.<br>  
 * 
 * This class produces a file with time series of pore sizes, averaged position of pore points and copolymers
 * around the centriod of pore in files named as "PoreSize.dat","PoreCoordinates.dat" and "PolymerCoordinates.dat"
 * respectively, with an option to add suffix to file name using an option in the constructor.
 * 
 * @tparam IngredientsType
 *
 * @param groups_ subgroup of molecules_type having groups of objects which translocate through bilayer.
 * @param ingredients_ the system holding the simulation box with bilayer, solvent and objects.
 * 
 * @param filesuffix_ option to add suffix to filename which are output at the end.
 *  
 **/

#define CLUSTER_SIZE_THRESHOLD 20
#define POLYMER_DIST_THRESHOLD 10
#define SOLVENTAG 5
#define TAILTAG 2

template<class IngredientsType>
class AnalyzerPoreFinder: public AbstractAnalyzer
{
protected:
    //!Typedefs for various underlying container holding the monomers
    //!and monomer subgroups.  
    typedef typename IngredientsType::molecules_type molecules_type;
    typedef std::vector< MonomerGroup< molecules_type> > group_type;
    IngredientsType& ingredients;
    const molecules_type& molecules;
    const group_type& groups;
    typedef MonomerGroup< molecules_type> group_type_ind;

    //!Box variables.    
    int32_t boxX,boxY,boxZ,boxXm_1,boxYm_1,boxZm_1;
    
    //!vectors storing VectorInt2 objects
    std::vector<VectorInt2> coordinatesOfPore;
    std::vector<VectorInt2> centriod;
    std::vector<VectorInt2> neigbours;
    std::vector<VectorInt2> clusterList;

    //!vectors various arrays used by functions below.
    
    std::vector<std::vector< bool > > tailOccupied;
    std::vector<std::vector< int > > coordinatesOfPoreStore;
    std::vector<std::vector< int > > coordinatesOfPolymers;
    std::vector<int> sizeOfCluster;
    std::vector< std::vector< int32_t > > poreSizeStore;

    //!string variables for name of output files.    
    std::string outputFilePoreSize,outputFilePore,outputFilePolymer;
    
    int32_t counterForPore,counterForPolymer;
    int32_t binforradius;
    int32_t referenceI;

public:
    AnalyzerPoreFinder(IngredientsType& i, const group_type& g,std::string fileSuffix_="");
    virtual ~AnalyzerPoreFinder(){}
    
    //!Function to fill lattice with tail monomer types.
    void LatticeFiller();
    
    //!Function to store lattice points which are occupied by tail monomers.
    void tailOccupiedFinder();
    
    //!Function for cluster analysis of points free of any tail monomers.
    void clusterAnalysis();
    
    //!Function used by clusterAnalysis() to find a lipid site which is free of lipids and not been checked already.    
    VectorInt2 lookForLipidFreeSites();
    
    //!Function used by clusterAnalysis() addNeigbours of an last unchecked values in the cluster list.
    void addNeigbours(int32_t lastUnChecked);
    
    //!Function used by clusterAnalysis() to store various values to the memory for later dump.
    void storeValues();
    
    //!Function finding the polymer coordinates wrt the centriod of pore.
    void polymerCoordinatesFinder();
    
    //!Function to find the center of mass of subgroup of monomers(generally objects).    
    VectorDouble3 centerOfMass(const group_type_ind& m);
    
    //!Function to find distances using minimum image convention.    
    int32_t reduceDistanceInPeriodicSpace( int32_t distance, int32_t period);
    
    //!General updater functions
    bool execute();
    void cleanup();
    void initialize(){};

};

/**
* @brief Constructor
**/
template<class IngredientsType>
AnalyzerPoreFinder<IngredientsType>::AnalyzerPoreFinder(IngredientsType& i, const group_type& g, std::string fileSuffix_)
:ingredients(i)
,molecules(i.getMolecules())
,groups(g)
{ 
    //Initate the box variables and counter for execute function 
    // and polymers as zero.
    
    counterForPore=0;
    counterForPolymer=0;
    boxX=ingredients.getBoxX();
    boxY=ingredients.getBoxY(); 
    boxZ=ingredients.getBoxZ();
    boxXm_1=ingredients.getBoxX()-1;
    boxYm_1=ingredients.getBoxY()-1; 
    boxZm_1=ingredients.getBoxZ()-1;
    
   //Add suffix if needed to added for output files.
       
    outputFilePoreSize="PoreSize"+fileSuffix_+".dat";
    outputFilePore="PoreCoordinates"+fileSuffix_+".dat";
    outputFilePolymer="PolymerCoordinates"+fileSuffix_+".dat";
    
   //Finding the index of the lipid monomer on the bilayer for reference to
   //calculate value of midplane. This is mainly done for case when bilayer is 
   //at the of periodic box in z direction(z=0 or z=boxZ-1).  
    
    for(int32_t i=0;i<ingredients.getMolecules().size();i++){
        int32_t type = ingredients.getMolecules()[i].getAttributeTag();
        if(type<SOLVENTAG) {referenceI=i; break;}}

   //Initate various array with box size in x-y direction.
    tailOccupied.resize(boxX,std::vector<bool>(boxY));
    coordinatesOfPoreStore.resize(boxX,std::vector<int32_t>(boxY,0));
    coordinatesOfPolymers.resize(boxX,std::vector<int32_t>(boxY,0));

   //Initate pore size array with three rows.
   //for MCS, poresize, and third one radius of gyration. 
    poreSizeStore.resize(3);
    
   //Initate nearest neigbours on 2d surface on simple cubic lattice.
    VectorInt2 V(1,0);
    neigbours.push_back(V);
    V.setAllCoordinates(-1,0);
    neigbours.push_back(V);
    V.setAllCoordinates(0,1);
    neigbours.push_back(V);
    V.setAllCoordinates(0,-1);
    neigbours.push_back(V);
}

/**
* @details Execution of analyzer where it finds the point filled with tail monomers, 
* then analyses the cluster of points free of tail monomers. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
bool AnalyzerPoreFinder<IngredientsType>::execute(){

    //Fill the lattice 
    //and find all the sites filled with tail monomers in x-y plane.
    LatticeFiller();
    tailOccupiedFinder();
    
    //Perform cluster analysis.     
    clusterAnalysis();
    
    //If no pore found, return true as there is no need to find polymer coordinates
    // in closed pore.
    if(coordinatesOfPore.size()==0)
        return true;
    
    //If there is a pore find where centerOfMass of copolymers
    //lie wrt centriod of the pore
    polymerCoordinatesFinder();
    
    return true;
}

/**
* @details Fill lattice if there is a tail monomer with its type or zero otherwise 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::LatticeFiller(){
  
    for(int32_t x=0;x<boxX;x++)
      for(int32_t y=0;y<boxY;y++)
          for(int32_t z=0;z<boxZ;z++)
                ingredients.setLatticeEntry(x,y,z,0);
 
for(int32_t i=0;i<ingredients.getMolecules().size();i++){
  VectorInt3 pos=ingredients.getMolecules()[i];
  int32_t type=ingredients.getMolecules()[i].getAttributeTag();
  
  if(type!=TAILTAG) continue;
  
  VectorInt3 dx(1,0,0);
  VectorInt3 dy(0,1,0);
  VectorInt3 dz(0,0,1);
  
  ingredients.setLatticeEntry(pos,type);
  ingredients.setLatticeEntry(pos+dx,type);
  ingredients.setLatticeEntry(pos+dy,type);
  ingredients.setLatticeEntry(pos+dx+dy,type);
  ingredients.setLatticeEntry(pos+dz,type);
  ingredients.setLatticeEntry(pos+dz+dx,type);
  ingredients.setLatticeEntry(pos+dz+dy,type);
  ingredients.setLatticeEntry(pos+dz+dx+dy,type);
 }
 
} 

/**
* @details Find the lattice sites in x-y plane which are occupied with tail monomers.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::tailOccupiedFinder(){
    
   for(int32_t x=0;x<boxX;x++)
       for(int32_t y=0;y<boxY;y++)
           for(int32_t z=0;z<boxY;z++)
               if(ingredients.getLatticeEntry(x,y,z)==TAILTAG){
                   tailOccupied[x][y]=true;
                   break;
                 }
               else
                   tailOccupied[x][y]=false;
}

/**
* @details Analyse the cluster of points to find pore points on lattice.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::clusterAnalysis(){
    //Restart the arrays with zero size
    coordinatesOfPore.resize(0);
    centriod.resize(0);
    sizeOfCluster.resize(0);
    sizeOfCluster.push_back(0);

    
     // Outside do-while loop looks for points on lattice in x-y plane 
     // which are free of lipids, once one such point is found, the inside 
     // loop finds the cluster associated with this point. At the end, all 
     // the clusters with the size bigger than CLUSTER_SIZE_THRESHOLD are saved 
     // in coordinatesofPore for further processing.   
     
    do{
        //Restart the clusterList
        clusterList.resize(0);
        //initial coordinate of point which is tail free 
        VectorInt2 iniCoordinate=lookForLipidFreeSites();
        
        //If lookForLipidFreeSites returns an impossible coordinate
        //this means all the lattice points have been looked upon.
        // thus break the loop. 
          
        if(iniCoordinate.getX()==(boxX+1)) break;
        
        //Otherwise add this clusterList
        clusterList.push_back(iniCoordinate);
        int32_t lastUnChecked=0;
         
        do{
             //check the for neigbours of the 
             //lastUnChecked cluster points which are tail free
            addNeigbours(lastUnChecked);
            //increase the counter of unchecked lattice sites.
            lastUnChecked++;
           //If the lastUnChecked has yet gone through 
            // all the clusterList, keep do-ing
        }while(lastUnChecked!=clusterList.size());

        //if clusterlist size is larger than threshold then save.
        if(clusterList.size()>CLUSTER_SIZE_THRESHOLD){
            sizeOfCluster.push_back(clusterList.size()+coordinatesOfPore.size());
            coordinatesOfPore.insert(coordinatesOfPore.end(), clusterList.begin(), clusterList.end());
        }
        
   }while(true);
   
   //Once the cluster is found, 
    // store the values.
   storeValues();
}

/**
* @details Adds neigbours of a lattice point which are free of tails.
* This function also sets tail occupied tag to true as an indicator 
* that the lattice site has been checked before.  
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::addNeigbours(int32_t lastUnChecked){
    
    VectorInt2 neigbourToPoint;
    int32_t x,y;
    VectorInt2 pointTobeChecked=clusterList[lastUnChecked];
    
    for(int32_t i=0;i<neigbours.size();i++){
        
        neigbourToPoint=pointTobeChecked+neigbours[i];
        
        //fold back in the case of boundary cases
        x=(neigbourToPoint[0])&boxXm_1;
        y=(neigbourToPoint[1])&boxYm_1;
        
        neigbourToPoint.setAllCoordinates(x,y);
        
         //If point was not checked before 
         //add to list and set the checked 
         //tag to true
         
        if(tailOccupied[x][y]==false){            
            clusterList.push_back(neigbourToPoint);
            tailOccupied[x][y]=true;
        }
    }
    
    tailOccupied[pointTobeChecked.getX()][pointTobeChecked.getY()]=true;
 }
 
/**
* @details Find the and return lattice sites which are tail free and have not been checked before, 
* If no point is found then send coordinates specific boxlength+1 to indicate to clusterAnalysis()
* that search is over.
*   
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
VectorInt2 AnalyzerPoreFinder<IngredientsType>::lookForLipidFreeSites(){
    
    for(int32_t x=0;x<boxX;x++)
        for(int32_t y=0;y<boxY;y++){
            if(tailOccupied[x][y]==false){
                tailOccupied[x][y]=true;
                VectorInt2 vec(x,y);
                return vec;
            }
            tailOccupied[x][y]=true;
        }
    
    VectorInt2 vec(boxX+1,boxY+1);
    return vec;
}

/**
* @details Finds the centriod of all the pores, then stores the values of pore size,
* pore coordinates wrt to centriod, and radius into container to dumped in a file at 
* cleanup. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::storeValues(){
    
    //Centriods calculation
    //Similar to midplane calculation here we need a reference point on 
    //x-y plane to cover for boundary cases.

   for(int32_t i=1;i<sizeOfCluster.size();i++){
       VectorInt2 centriodSum(0,0);
       for(int32_t j=sizeOfCluster[i-1];j<sizeOfCluster[i];j++){
          
           int32_t x = coordinatesOfPore[j][0]-coordinatesOfPore[sizeOfCluster[i-1]][0];
           int32_t y = coordinatesOfPore[j][1]-coordinatesOfPore[sizeOfCluster[i-1]][1];
           x= reduceDistanceInPeriodicSpace(x,boxX);
           y= reduceDistanceInPeriodicSpace(y,boxY);
          
           VectorInt2 Vec(x,y);
           centriodSum = centriodSum + Vec;
    }
   
    centriodSum = centriodSum/(sizeOfCluster[i]-sizeOfCluster[i-1]);
    centriod.push_back(centriodSum);
   }
   //reshift the coordinates according to box coordinates along with the fold back.
   
   for(int32_t i=0;i<centriod.size();i++){
       centriod[i].setX((centriod[i][0]+coordinatesOfPore[sizeOfCluster[i]][0])&boxXm_1);
       centriod[i].setY((centriod[i][1]+coordinatesOfPore[sizeOfCluster[i]][1])&boxYm_1);
   }
   
   //find the radius of gyration and store coordinates of pore wrt Centriods.
   double radiusofGyration=0;
   for(int32_t i=1;i<sizeOfCluster.size();i++){
       double rG_temp=0;
       for(int32_t j=sizeOfCluster[i-1];j<sizeOfCluster[i];j++){
           int32_t x= reduceDistanceInPeriodicSpace(coordinatesOfPore[j][0]-centriod[i-1][0],boxX);
           int32_t y= reduceDistanceInPeriodicSpace(coordinatesOfPore[j][1]-centriod[i-1][1],boxY);
           VectorInt2 Vec(x,y);
           //to take care of the negative values
           coordinatesOfPoreStore[x+boxX/2][y+boxY/2]++;
           rG_temp += Vec.getLength()*Vec.getLength();
       }
       counterForPore++;
       rG_temp /= (sizeOfCluster[i]-sizeOfCluster[i-1]);
       rG_temp = sqrt(rG_temp);
       radiusofGyration += rG_temp;
    }
    if(sizeOfCluster.size())
        radiusofGyration /= (sizeOfCluster.size()-1);
   
   //finds the weighted average of pores size(where weight is pore size 
   //itself) to be stored in container which dumps in a file.
   float squarePoreSize=0,poreSize=0;
   for(int32_t i=1;i<sizeOfCluster.size();i++){
       squarePoreSize += (sizeOfCluster[i]-sizeOfCluster[i-1])*(sizeOfCluster[i]-sizeOfCluster[i-1]);
       poreSize += (sizeOfCluster[i]-sizeOfCluster[i-1]);
   }

   if(poreSize)
       squarePoreSize /= poreSize;
   
   poreSizeStore[0].push_back(ingredients.getMolecules().getAge());
   poreSizeStore[1].push_back(squarePoreSize);
   poreSizeStore[2].push_back(radiusofGyration);
   
}

/**
* @details Finds coordinates of centerOfMass of copolymers 
* with respect to the centriod of pore, only considering the polymer 
* which POLYMER_DIST_THRESHOLD away from the midplane.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::polymerCoordinatesFinder(){
    /**find the midplane */
    int32_t midplane=0,lipidcounter=0;
    for(int32_t i=0;i<ingredients.getMolecules().size();i++){
        int32_t type = ingredients.getMolecules()[i].getAttributeTag();
        if(type>=SOLVENTAG) continue;
        midplane += (reduceDistanceInPeriodicSpace(ingredients.getMolecules()[i].getZ()-ingredients.getMolecules()[referenceI].getZ(),boxZ));
        lipidcounter++;
    }
    
    midplane /= lipidcounter;
    midplane = (midplane+ingredients.getMolecules()[referenceI].getZ())&boxZm_1;
    //store polymer COM wrt the centriod.
    for(int32_t CN=0;CN<centriod.size();CN++){
        for(int32_t i=0;i<groups.size();i++){
            VectorDouble3 groupCOM=centerOfMass(groups[i]);
            int32_t distance=reduceDistanceInPeriodicSpace(midplane-groupCOM[2],boxZ);
            
            if(fabs(distance)<POLYMER_DIST_THRESHOLD){
                int32_t x=reduceDistanceInPeriodicSpace(groupCOM[0]-centriod[CN][0],boxX);
                int32_t y=reduceDistanceInPeriodicSpace(groupCOM[1]-centriod[CN][1],boxY);
                //to take care of the negative values
                coordinatesOfPolymers[x+boxX/2][y+boxY/2]++;
            }
        }
        counterForPolymer++;
    }
    
}

/**
* @details Calculates the center of mass of objects given a subgroup of 
* monomers of the object.
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template < class IngredientsType > 
VectorDouble3 AnalyzerPoreFinder<IngredientsType>::centerOfMass(const group_type_ind& m)
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
int32_t AnalyzerPoreFinder<IngredientsType>::reduceDistanceInPeriodicSpace(int32_t distance, int32_t period)
{
    
	if(distance > period/2){
            while(-(distance-period)<period/2) 
                distance -= (period);
            
            return distance;}
	
	else if(-distance > period/2){ 
            while((distance+period)<period/2) 
                distance += (period);
            
            return distance;}
	
	else 
            return distance;
}
/**
* @details dump all containers into files. 
* 
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void AnalyzerPoreFinder<IngredientsType>::cleanup(){
    
    std::stringstream comments;
    comments<<"Created by AnalyzerPoreFinder\n";
    comments<<"file contains time series of weighted average pore size and radius of gyration of pore(second cumulant of average pore size) \n";
    comments<<"format: mcs\t PoreSize\t Radius_Gyration\n";

    ResultFormattingTools::writeResultFile(
            outputFilePoreSize,
            ingredients,
            poreSizeStore,
            comments.str()
    );
    
    comments.str("");
    comments<<"Created by AnalyzerPoreFinder\n";
    comments<<"file contains averaged pore points from centriod of a pore\n";
    comments<<"format: x\t y\t Average_value\n";
    std::vector< std::vector< float > > coordinatesOfPoreFormat;
    coordinatesOfPoreFormat.resize(3);
    
    for(int32_t x=0;x<coordinatesOfPoreStore.size();x++)
        for(int32_t y=0;y<coordinatesOfPoreStore[x].size();y++){
            coordinatesOfPoreFormat[0].push_back(x);
            coordinatesOfPoreFormat[1].push_back(y);
            coordinatesOfPoreFormat[2].push_back((float)coordinatesOfPoreStore[x][y]/counterForPore);
        }
    ResultFormattingTools::writeResultFile(
            outputFilePore,
            ingredients,
            coordinatesOfPoreFormat,
            comments.str()
    );
    
    comments.str("");
    comments<<"Created by AnalyzerPoreFinder\n";
    comments<<"file contains averaged position of copolymers from centriod of a pore\n";
    comments<<"format: x\t y\t Average_value\n";
    
    if(counterForPolymer){
        std::vector< std::vector< float > > coordinatesOfPolyFormat;
        coordinatesOfPolyFormat.resize(3);
        for(int32_t x=0;x<coordinatesOfPolymers.size();x++)
            for(int32_t y=0;y<coordinatesOfPolymers[x].size();y++){
                coordinatesOfPolyFormat[0].push_back(x);
                coordinatesOfPolyFormat[1].push_back(y);
                coordinatesOfPolyFormat[2].push_back((float)coordinatesOfPolymers[x][y]/counterForPolymer);
            }
            
        ResultFormattingTools::writeResultFile(
            outputFilePolymer,
            ingredients,
            coordinatesOfPolyFormat,
            comments.str()
    );
    }
}
#endif
