#ifndef UpdaterLipidsCreator_H
#define UpdaterLipidsCreator_H

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

/**
 * @file
 *
 * @class UpdaterLipidsCreator
 *
 * @brief Updater to create bilayer in the middle of the box along z axis with solvent and/or nano particles 
 * and facially amphiphillic copolymers at random positions.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients.
 * This updater requires FeatureAttributes. There are choices of two types of head structure (i) Ring shaped
 * (ii) linear polymer with 3 monomers as was used in Werner et. al., Soft Matter 8:11714,2012. 
 * There is a possibility to create two types with different tail lengths specified by n1_ and n2_.<br>
 * \b Attributes: <br>
 * \b Lipids_shorter_length: Head monomers Tag 1 and Tail monomers Tag 2.<br>
 * \b Lipids_longer_length: Head monomers Tag 3 and Tail monomers Tag 4.<br>
 * \b Solvent_monomers: Tag 5.<br>
 * \b Nano_particles: Tag 7.<br>
 * \b Facially_Amphiphlic_Copolymers: backbone monomers Tag 7 and side group monomers Tag 8.<br>
 * 
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding an empty simulation box for system setup.
 * 
 * @param nlipids_ number of lipids required in bilayer.
 * @param Bilayer_Half_ if needed a half membranes(covers the 75% of xy plane).
 * @param ringHead if needed a ring shaped head instead of normal three monomers.
 * @param MeanLatticeOccupany_ Mean Lattice Occupany of simple cubic lattice(default 0.5).
 * @param n1_ is length of tails of lipid type 1.
 * @param n2_ is length of tails of lipid type 2.
 * @param n_fraction_ fraction of lipid type 1.
 * @param NanoType_ if Nano particles are to be included type 0 then tetrahedral nano particle 
 * would be created otherwise triangular shaped.
 * @param NanoNum_ number of nano particles in the system.
 * @param polylength_ length of facially amphiphillic copolymers to be created.
 * @param numpoly_ number of facially amphiphillic copolymers to be created.
 **/

using namespace std;

template<class IngredientsType>
class UpdaterLipidsCreator:public AbstractUpdater
{
public:
    UpdaterLipidsCreator(IngredientsType& ingredients_,int32_t nlipids_=300,double MeanLatticeOccupany_=0.5, int32_t n1_=5, int32_t n2_=0, double n_fraction_=0,bool ringHead_=false,bool Bilayer_Half_=false,int32_t NanoType_=1,int32_t NanoNum_=0,int32_t polylength_=5,int32_t numpoly_=0);

    //!General updater functions
    bool execute();
    void initialize(){execute();};
    void cleanup(){};   
    
    
    //!Function used by function lipidsBilayerCreator to find a vacant position where lipid has to be created.     
    VectorInt3 vacantLipidPosFinder(MoveAddMonomerSc<int32_t>& move,int32_t n1,int32_t upperleaflet_);
    
    //!Function used by function lipidsBilayerCreator to create a single lipid.     
    void singleLipidCreator(MoveAddMonomerSc<int32_t>& move, VectorInt3 pos,int32_t TotalMonLipids,int32_t upperleaflet);


    //!Function for creating lipid bilayer in the middle of the box.
    void lipidsBilayerCreator(MoveAddMonomerSc<int32_t>& move);
    
    //!Function used by the functions faciallyAmphCopolymers, nanoTriangle and nanoTetrahedron
    //!to find shapes vacant for objects to be created.     
    VectorInt3 vacantShapeFinder(MoveAddMonomerSc<int32_t>& move, std::vector<VectorInt3> shape);

    //!Function for creating Facially Amphiphlic Copolymers.
    void faciallyAmphCopolymers(MoveAddMonomerSc<int32_t>& move,int32_t polylength, int32_t numpoly);
    
    //!Function for creating Nano particles in the triangular shape.    
    void nanoTriangle(MoveAddMonomerSc<int32_t>& move,int32_t TotNanoNum_);

    //!Function for creating Nano particles in the tetrahedral shape.
    void nanoTetrahedron(MoveAddMonomerSc<int32_t>& move, int32_t TotNanoN_);
    
    //!Function for creating ring shaped head on the lipid molecules.
    void ringHeadCreator(MoveAddMonomerSc<int32_t>& move,VectorInt3 pos, int32_t upperleaflet_, int32_t type);
    
    //!Function for creating solvent when all other particles have been created.
    void solventCreator(MoveAddMonomerSc<int32_t>& move, double solventdens_);
    
    //!Function to get upperleaflet lipids.     
    int32_t getUpperleafletLipids(){return leafletcounter;}

protected:
    //!Ingredients holding the simulation box
    IngredientsType& ingredients;

    //!Number of lipids to be created
    int32_t nlipids;
    
    //!Bool variable for whether half bilayer is needed or ring head is needed in the lipids.
    bool Bilayer_Half,ringHead;
    
    //!Mean lattice Occupancy of the lattice.
    double MeanLatticeOccupany;

    //!Mean lattice Occupancy of the lattice.
    std::vector<VectorInt3> upperleafletBondvecs,lowerleafletBondvecs;
    
    RandomNumberGenerators randomNumbers;

    int32_t leafletcounter;
    int32_t n1,n2,shorterlipid_length,longerlipid_length,shorterlipid_number,longerlipid_number;
    int32_t NanoType,NanoNumber,polylength,numpoly;
    double n_fraction;
};

/**
* @brief Constructor
**/

template<class IngredientsType>
UpdaterLipidsCreator<IngredientsType>::UpdaterLipidsCreator(IngredientsType& ingredients_,int32_t nlipids_,double MeanLatticeOccupany_, int32_t n1_, int32_t n2_, double n_fraction_, bool ringHead_,bool Bilayer_Half_, int32_t NanoType_,int32_t NanoNum_,int32_t polylength_,int32_t numpoly_)
    :ingredients(ingredients_),
    nlipids(nlipids_),
    MeanLatticeOccupany(MeanLatticeOccupany_),
    Bilayer_Half(Bilayer_Half_),
    ringHead(ringHead_),
    n1(n1_),
    n2(n2_),
    n_fraction(n_fraction_),
    NanoType(NanoType_),
    NanoNumber(NanoNum_),
    polylength(polylength_),
    numpoly(numpoly_)
{

    upperleafletBondvecs.push_back(ingredients.getBondset().getBondVector(21));
    upperleafletBondvecs.push_back(ingredients.getBondset().getBondVector(33));
    upperleafletBondvecs.push_back(ingredients.getBondset().getBondVector(27));

    lowerleafletBondvecs.push_back(ingredients.getBondset().getBondVector(18));
    lowerleafletBondvecs.push_back(ingredients.getBondset().getBondVector(30));
    lowerleafletBondvecs.push_back(ingredients.getBondset().getBondVector(24));
    
    randomNumbers.seedAll();

    /** The choice of these bond vectors implies two conditions:
    * 
    * (i) the height of lipids in z axis would be lipid_length+2.
    * 
    * (ii) all the lipids are created 4 lattice units away from each other 
    * to avoid the excluded volume conflicts. This restricts the number of 
    * lipids that can be filled, but as it turns out maximum number of lipids 
    * which are filled this way already voilates the "almost flat condition", 
    * which is generally required by the bilayer related simulations, therefore, 
    * these bond vectors suffice the simulations of bilayer.
    */
    
    /********************************************************************************/
    
    /**If both n1 and n2 are present decide what are the number and total length of shorter 
     * and longer lipids (assuming three head monomers). There are actually
     * 4 head monomers in the case of ring head. But this case is separately handled
     * function singleLipidCreator.
     */
    if(n1&&n2){
        shorterlipid_length=(n1<n2)?2*n1+3:2*n2+3;
        longerlipid_length=(n1>n2)?2*n1+3:2*n2+3;
        shorterlipid_number=(shorterlipid_length==2*n1+3)?(nlipids*n_fraction):((1-n_fraction)*nlipids);
        longerlipid_number=(longerlipid_length==2*n1+3)?(nlipids*n_fraction):((1-n_fraction)*nlipids);
    }
    
    //if one of them is zero there is only one type of lipid,
    //assign its important properties to shorterlipid_length 
    else{
        shorterlipid_length=(n1)?2*n1+3:2*n2+3;
        longerlipid_length=shorterlipid_length;
        shorterlipid_number=nlipids+1;
        longerlipid_number=-1;}
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template<class IngredientsType>
bool UpdaterLipidsCreator<IngredientsType>::execute()
{
   
  MoveAddMonomerSc<int32_t> move;
  
  //First create the Lipids
  lipidsBilayerCreator(move);
  
  //Create facially Amphiphlic Copolymers
  faciallyAmphCopolymers(move,polylength, numpoly);
  
  //Create Nano particles
  if(NanoType)
   nanoTriangle(move, NanoNumber);
  else
   nanoTetrahedron(move, NanoNumber);
  
  //Lastly create Solvent
  solventCreator(move, MeanLatticeOccupany);
  return true;
}

/**
* @brief Creates a bilayer with midplane at the middle of z axis of the box. 
* It distributes equal number of lipids in both leaflets of bilayer.
* Also, it creates n_fraction of n1 lipids.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::lipidsBilayerCreator(MoveAddMonomerSc<int32_t>& move)
{
    //Initialize the variables to calculate number of lipids in leaflet and number
    //of types of lipids. 
    int32_t typeofLipid;
    int32_t shorterlipidcounter=0,longerlipidcounter=0;
    int32_t leaflet;
    leafletcounter=0;
    bool Oneleaffull=false,OneofTypeComplete=false;
    VectorInt3 InitialPos;
   
    for(int32_t lipidnumber =0;lipidnumber<nlipids;lipidnumber++){
        
         //Randomly choose a leaflet(upper=1 and lower=0) randomly 
         //until it's filled by half of lipids.
        
        if(Oneleaffull==false)
            leaflet=rand()%2;
        
        //Randomly choose a type until one of them is completely created. 
         
        if(OneofTypeComplete==false)
            typeofLipid=(rand()%2)?shorterlipid_length:longerlipid_length;

        //Find the vacant position and then fill it with lipid. 
         
        InitialPos=vacantLipidPosFinder(move,typeofLipid,leaflet); 
        singleLipidCreator(move,InitialPos,typeofLipid,leaflet);
        
        //Count the both type of lipids.
         
        if(typeofLipid==shorterlipid_length)
            shorterlipidcounter++;
        else
            longerlipidcounter++;  
        
        //Count the number of lipids in the upper leaflet.
         
        if(leaflet) 
            leafletcounter++;
        
        //Tells the system that one of the leaflet is full.
         
        if(leafletcounter==nlipids/2||(lipidnumber-leafletcounter+1==nlipids/2)){
            Oneleaffull=true;
            leaflet=(leafletcounter==nlipids/2)?0:1;
            }
        
        //Tells the system that creation one of the type of lipids is finished 
         
        if(shorterlipidcounter==(shorterlipid_number)||longerlipidcounter==(longerlipid_number)){
            OneofTypeComplete=true;
            typeofLipid=(shorterlipidcounter==shorterlipid_number)?longerlipid_length:shorterlipid_length;
            }
        }        
};	

/**
* @brief Finds a vacant position for a lipids given the type of lipid(shorter or longer) and leaflet(upper or lower)s. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to check if the lattice sites are empty adding monomers later.
* @param lipid_type shorter or longer lipid.
* @param leaflet upper or lower leaflet. 
*/

template<class IngredientsType>
VectorInt3 UpdaterLipidsCreator<IngredientsType>::vacantLipidPosFinder(MoveAddMonomerSc<int32_t>& move,int32_t lipid_type,int32_t leaflet){
   
    //Initializing a move to check for vacant position in lattice.
    move.init(ingredients);
    VectorInt3 position;
    int32_t x,y,z;
    
    // \grow decides if lipid grows upwards or downwards (lower or upper leaflets respectively).
    int32_t grow=(leaflet)?1:-1;
    
    //Here, we choose shorterlipid_length as if position is vacant for shorter lipids,
    //it is also vacant for longer lipids.  
      
    int32_t height_lipids=shorterlipid_length+2;

    //do-while loop performs till a vacant position is found
      
       do{ 
           
         //In the case of bilayer half only 3/4 of x-axis is covered
         //the factor 1/4 comes due to the choice of bond vectors.
         // See the explanation in constructor.
        
         if(Bilayer_Half)
             x=rand()%(3*(ingredients.getBoxX())/16);
         else
             x=rand()%((ingredients.getBoxY())/4);  
         
         x=x*4;

         //In the case of ringHead, head monomers covers more lattice
         //units on y-axis. 
         
         if(ringHead){ 
             y=rand()%((ingredients.getBoxY())/6);
             y=y*6;}
         else{
            y=rand()%((ingredients.getBoxY())/4);
            y=y*4;}
            
        //The point on the z axis where lipid starts growing from.
        
         z=(ingredients.getBoxZ()/2)+height_lipids*grow;
         
         position.setAllCoordinates(x,y,z);
         move.setPosition(position);
         }while(move.check(ingredients)==false);
        
     //Change z position if longer lipids have to be grown.
        
    if(lipid_type==shorterlipid_length)
        return position;
    
    else{
        z=(ingredients.getBoxZ()/2)+(longerlipid_length+2)*grow;
        position.setAllCoordinates(x,y,z);   
        return position;
    }
}	

/**
* @brief Creates a single lipid on given any intial position, total monomers present in lipid and leaflet(upper or lower). 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to apply or create a single lipid.
* @param Initpos intial position of lipid from where lipid starts growing.
* @param TotalMonLipids Total monomers in the lipid to be grown.
* @param leaflet upper or lower leaflet. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::singleLipidCreator(MoveAddMonomerSc<int32_t>& move,VectorInt3 initPos,int32_t totalMonLipids,int32_t leaflet){
       
       std::vector<VectorInt3> bondVecSet;
       
       if(leaflet)
           bondVecSet=upperleafletBondvecs;
       else
           bondVecSet=lowerleafletBondvecs;            
           
        // Initialize tail length variable which length of tails which can calculated by
        // assuming three head monomers. 
       
       int32_t tail_length=(totalMonLipids-3)/2,lastMonomerAdded;
 
       //type_chr handles the tag to be given to the head and tail monomers
       //of different kind of lipids.
       
       int32_t type_chr=(totalMonLipids==shorterlipid_length)?1:3;
       
       VectorInt3 position;
       position=initPos;
       VectorInt3 bondVector(0,0,0);
       
       for(int32_t monNum=0;monNum<totalMonLipids;monNum++){
   
           //type assigns different tag to head monomers(monNum<3)
           //and tail monomers

           int32_t type=(monNum<3)? type_chr : type_chr+1;         
           move.init(ingredients);
           move.setTag(type);
           move.setPosition(position+bondVector);
           
           //The if-else statement is safety measure to make sure that initPos is selected
           //correctly.
           
           if(move.check(ingredients)==true)
           {               
                position+=bondVector;
                //if ring head is required
                if(ringHead&&monNum==0){
                ringHeadCreator(move, initPos,leaflet,type);
                monNum=monNum+3;
                bondVector=bondVecSet[1];
                lastMonomerAdded = ingredients.getMolecules().size()-1;
                position=ingredients.getMolecules()[lastMonomerAdded];
                position+=bondVector;
                move.setPosition(position);
                move.setTag(type_chr+1);
                }
                
                move.apply(ingredients);
                lastMonomerAdded = ingredients.getMolecules().size()-1;
                bondVector=bondVecSet[0];
                
                //Always connect the last monomer to second last monomer unless
                //it is tail ends.  
                 
                if ( monNum > 0&& monNum!=(totalMonLipids-tail_length))      
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);                  
            
                // Choosing the different bondvector if it is the last head monomer.
                
                if(monNum==2){                                                      
                bondVector=bondVecSet[1];
                }

                //Connection to second tail and last head monomer.
                 
                if(monNum==totalMonLipids-tail_length)                          
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-tail_length-1);
                
                //Change the position and bondvector at the end first tail.
                 
                if(monNum==totalMonLipids-tail_length-1){
                position=ingredients.getMolecules()[lastMonomerAdded-tail_length];
                bondVector=bondVecSet[2]; 
                 }
                }
                
            else
                throw std::runtime_error("Excluded volume conflict: could not add a single lipid");
	}
      }

/**
* @brief Creates a ring shape head in lipids, given the position, leaflet(upper or lower) and type. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to apply or create a ring head.
* @param Initpos intial position of lipid; from where lipid starts growing.
* @param leaflet upper or lower leaflet. 
* @param type tag assigned to the monomers. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::ringHeadCreator(MoveAddMonomerSc<int32_t>& move, VectorInt3 Initpos,int32_t leaflet,int32_t type){
    
    move.init(ingredients);
    int32_t Factor=(leaflet)?-1:1;
    VectorInt3 v1(1,0,Factor*2);
    VectorInt3 v2(-1,0,Factor*2);
    VectorInt3 v3(0,2,Factor*2);
    VectorInt3 v4(0,0,Factor*4);
    move.setPosition(Initpos);
    move.setTag(type);
    move.apply(ingredients);
    
    move.setPosition(Initpos+v1);
    move.setTag(type);
    move.apply(ingredients);
    int32_t lastMonomerAdded = ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
    
    move.setPosition(Initpos+v2);
    move.setTag(type);
    move.apply(ingredients);
    lastMonomerAdded = ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-2);

    move.setPosition(Initpos+v4);
    move.setTag(type);
    move.apply(ingredients);
    lastMonomerAdded = ingredients.getMolecules().size()-1;
    ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
    
    ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-2);
}

/**
* @brief Finds the vacant shape on avoiding excluded volume conflicts. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to find a vacant shape on the lattice.
* @param shape vector of 3d coordinates of shape which has to be found on the lattice. 
*/


template<class IngredientsType>
VectorInt3 UpdaterLipidsCreator<IngredientsType>::vacantShapeFinder(MoveAddMonomerSc<int32_t>& move, std::vector<VectorInt3> shape)
{
     //x,y,z are the initial coordinates of shape which would be returned later. 
     //Midplane and ShiftValue are found so that they do not have any conflicts with the lipids.
     
    int32_t x,y,z,counteroutside=0;
    VectorInt3 position;

    bool occupied;
            
      //Do till an unoccupied position is found.
     
    do{
        occupied =false;
       
        //Random initial position
        x=rand()%(ingredients.getBoxX());
        y=rand()%(ingredients.getBoxY());
        z=rand()%(ingredients.getBoxZ());
        
        move.init(ingredients);
        position.setAllCoordinates(x,y,z);
        move.setPosition(position);
        
        //Find if the object can inserted at this place
        for(int32_t monNum=0;monNum<shape.size();monNum++){
            move.init(ingredients);
            move.setPosition(position+shape[monNum]);
            if(move.check(ingredients)==false) {occupied= true; break;}
            }
            
    counteroutside++;
    //Throw error if too many objects are inserted.
    if(counteroutside==10000)
        throw std::runtime_error("Could not find vacant position for object; too many objects");
    
    }while(occupied==true);
    
    return position;
}

/**
* @brief Creates Facially Amphiphlic Copolymers given length and number of polymers. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to find a vacant shape on the lattice.
* @param polylength length of copolymers. 
* @param numpoly number of copolymers. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::faciallyAmphCopolymers(MoveAddMonomerSc<int32_t>& move,int32_t polylength, int32_t numpoly)
{
    //Create a vector of 3d coordinates to store shape of copolymer.

    std::vector<VectorInt3> shapevector;
    VectorInt3 vec;
    
    for(int32_t polyMono=0;polyMono<polylength;polyMono++){
        vec.setAllCoordinates(2*polyMono,0,0);
        shapevector.push_back(vec);
        vec.setAllCoordinates(2*polyMono,2,0);
        shapevector.push_back(vec);
    }
    
    VectorInt3 Initpos,position;
    for(int32_t i=0;i<numpoly;i++){
        //Find a vacant space
        Initpos=vacantShapeFinder(move,shapevector);
        int32_t x=Initpos.getX();
        int32_t y=Initpos.getY();
        int32_t z=Initpos.getZ();
        
        //Create the copolymer based on InitialPos.
        int32_t lastMonomerAdded;
        for(int32_t polyMono=0;polyMono<polylength;polyMono++){
            
            position.setAllCoordinates(x+2*polyMono,y,z);
            move.setPosition(position);
            move.setTag(7);
            move.apply(ingredients);

            position.setAllCoordinates(x+2*polyMono,y+2,z);
            move.setPosition(position);
            move.setTag(8);
            move.apply(ingredients);
            
            lastMonomerAdded= ingredients.getMolecules().size()-1;
            if(polyMono){
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
                ingredients.modifyMolecules().connect(lastMonomerAdded-1,lastMonomerAdded-3);
            }
            else
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
        }
    }
}

/**
* @brief Creates nano particles in a triangular shape given number of nano particles. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used create the nano particles.
* @param totNanoNum number of nano particles. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::nanoTriangle(MoveAddMonomerSc<int32_t>& move,int32_t totNanoNum)
{   
    // Create a vector of 3d coordinates to store shape of nano particle.
    
    VectorInt3 v1(2,1,0);
    VectorInt3 v2(1,0,2);    

    std::vector<VectorInt3> shapevector;
    VectorInt3 vec(0,0,0);
    shapevector.push_back(vec);
    shapevector.push_back(vec+v1);
    shapevector.push_back(vec+v2);
    
    VectorInt3 Initpos;
    
    for(int32_t nanoNumber=0;nanoNumber<totNanoNum;nanoNumber++){

        //Find a vacant space.
        Initpos=vacantShapeFinder(move,shapevector);
        int32_t x=Initpos.getX();
        int32_t y=Initpos.getY();
        int32_t z=Initpos.getZ();
        
        //create the particle.
        move.setPosition(Initpos);
        move.setTag(7);
        move.apply(ingredients);
        
        move.setPosition(Initpos+v1);
        move.setTag(7);
        move.apply(ingredients);
        
        int32_t lastMonomerAdded = ingredients.getMolecules().size()-1;
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
    
        move.setPosition(Initpos+v2);
        move.setTag(7);
        move.apply(ingredients);
        
        lastMonomerAdded = ingredients.getMolecules().size()-1;
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-2);
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
        
  }
}

/**
* @brief Creates nano particles in a tetrahedral shape given number of nano particles. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used create the nano particles.
* @param totNanoNum number of nano particles. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::nanoTetrahedron(MoveAddMonomerSc<int32_t>& move,int32_t totNanoNum)
{   
    VectorInt3 v1(1,0,2);
    VectorInt3 v2(-1,0,2);
    VectorInt3 v3(0,2,1);

    std::vector<VectorInt3> shapevector;
    VectorInt3 vec(0,0,0);
    shapevector.push_back(vec);
    shapevector.push_back(vec+v1);
    shapevector.push_back(vec+v2);
    shapevector.push_back(vec+v3);
    
    VectorInt3 Initpos;
    
   for(int32_t nanoNumber=0;nanoNumber<totNanoNum;nanoNumber++){

        //Find a vacant space.
        Initpos=vacantShapeFinder(move,shapevector);
        int32_t x=Initpos.getX();
        int32_t y=Initpos.getY();
        int32_t z=Initpos.getZ();
        
        //create the particle.
        move.setPosition(Initpos);
        move.setTag(7);
        move.apply(ingredients);
        
        move.setPosition(Initpos+v1);
        move.setTag(7);
        move.apply(ingredients);
        
        int32_t lastMonomerAdded = ingredients.getMolecules().size()-1;
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
    
        move.setPosition(Initpos+v2);
        move.setTag(7);
        move.apply(ingredients);
        
        lastMonomerAdded = ingredients.getMolecules().size()-1;
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-2);
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
        
        move.setPosition(Initpos+v3);
        move.setTag(7);
        move.apply(ingredients);
        
        lastMonomerAdded = ingredients.getMolecules().size()-1;
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-2);
        ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-3);
        
  }
}   
      
/**
* @brief Creates explicit solvent particles given number mean lattice occupancy. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used create the solvent particles.
* @param meanLatticeOccupany mean lattice occupancy of lattice. 
*/

template<class IngredientsType>
void UpdaterLipidsCreator<IngredientsType>::solventCreator(MoveAddMonomerSc<int32_t>& move,double meanLatticeOccupany)
{

    //Initialize the variables required to create solvent.
    //Calculate the number of solvent particles needed to be created in order to maintain the given
    // mean lattice occupancy.
    
    int32_t x,y,z;
    int32_t total_monomers_added = ingredients.getMolecules().size();
    size_t numsol = (((ingredients.getBoxX())*(ingredients.getBoxY())*(ingredients.getBoxZ())*MeanLatticeOccupany)-(total_monomers_added*8))/8;
    
    int32_t initIndex=total_monomers_added-1;

    size_t i=0;
    int32_t Midplane=ingredients.getBoxZ()/2;
    
    while(i<numsol){
        x=rand()%(ingredients.getBoxX());
        y=rand()%(ingredients.getBoxY());
        z=rand()%(ingredients.getBoxZ());

        //It was found that creating solvent far away from the midplane and lipids 
        //prevents system to go to meta-stable states before reaching equilibrium state.
        
        if(Bilayer_Half){
        z=rand()%(ingredients.getBoxZ());
        if(x<3*ingredients.getBoxX()/4){
            if(i<numsol/2)
                z=rand()%(Midplane-10);
            if(i>numsol/2)
                z=rand()%(Midplane-10) + (Midplane+10);
            } 
        }
        
        else{
            if(i<numsol/2)
                z=rand()%(Midplane-10);
            if(i>numsol/2)
                z=rand()%(Midplane-10) + (Midplane+10);}

        VectorInt3 position(x,y,z);
        move.setTag(5);
        move.setPosition(position);
        if(move.check(ingredients)==true){
            move.apply(ingredients);
            i++;
        }   
        int32_t lastMonomerAdded=ingredients.getMolecules().size()-1;
        ingredients.setCompressedOutputIndices(initIndex+1,lastMonomerAdded); 
        }
};
	
#endif	 
