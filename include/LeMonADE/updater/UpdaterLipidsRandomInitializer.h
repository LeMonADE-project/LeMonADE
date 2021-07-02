#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterLipidsCreator.h>
/**
 * @file
 *
 * @class UpdaterLipidsRandomInitializer
 *
 * @brief Updater to create lipids at random position and random orientations with choice of two types of lipid
 * solvent monomers.
 *
 * @details This is upgrade to UpdaterLipidsCreator where instead of creating the lipids in form of bilayer, 
 * lipids are generated at random postions and orientations.There is a possibility to create two types 
 * with different tail lengths specified by n1_ and n2_. <br> 
 * \b Attributes: <br>
 * \b Lipids_shorter_length: Head monomers Tag 1 and Tail monomers Tag 2.<br>
 * \b Lipids_longer_length: Head monomers Tag 3 and Tail monomers Tag 4.<br>
 * \b Solvent_monomers: Tag 5.<br>
 * 
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding an empty simulation box for system setup.
 * 
 * @param nlipids_ number of lipids required in bilayer.
 * @param MeanLatticeOccupany_ Mean Lattice Occupany of simple cubic lattice(default 0.5).
 * @param n1_ is length of tails of lipid type 1.
 * @param n2_ is length of tails of lipid type 2.
 * @param n_fraction_ fraction of lipid type 1.
 * 
 **/

template<class IngredientsType>
class UpdaterLipidsRandomInitializer:public UpdaterLipidsCreator<IngredientsType>
{
    typedef UpdaterLipidsCreator<IngredientsType> BaseClass;
    typedef UpdaterSimpleSimulator<IngredientsType, MoveLocalSc> SimpSimu;
    
public:
    UpdaterLipidsRandomInitializer();
    //!Allowing access to constructors of UpdaterLipidsCreator.       
    using BaseClass::UpdaterLipidsCreator;
    //!Allowing access to other function of this class.       
    using BaseClass::solventCreator;
    using BaseClass::ingredients;
    //general updater functions
    void initialize();
    bool execute();
    void cleanup(){};
    
    //!Function for initiating lipid molecules at random position in the box.   
    void randomInitializer(MoveAddMonomerSc<int32_t>& move);
    
    //!Function for creating a single lipid in random orientation.   
    void singleLipidCreator(MoveAddMonomerSc<int32_t>& move,int32_t totalMonLipids,SimpSimu& simulater);
    
    //!Function draw a random bond vector which are allowed by BFM creatorian.   
    VectorInt3 drawRandomDirection() const;
    
private:
    
    //!Redecleration of the variables from the base class to be used in this class.       
    int32_t longerlipid_length,shorterlipid_length,nlipids,n_fraction,shorterlipid_number,longerlipid_number;
    float MeanLatticeOccupany;
    bool rI;
    
};

/**
* @brief Default constructor
**/
template<class IngredientsType>
UpdaterLipidsRandomInitializer<IngredientsType>::UpdaterLipidsRandomInitializer(){
    std::cout<<"hello";
    exit(0);
};

template<class IngredientsType>
void UpdaterLipidsRandomInitializer<IngredientsType>::initialize()
{
    //Use the base class variables in this class.
    this->nlipids = BaseClass::nlipids;
    this->shorterlipid_length = BaseClass::shorterlipid_length;
    this->longerlipid_length = BaseClass::longerlipid_length;
    this->shorterlipid_number = BaseClass::shorterlipid_number;
    this->longerlipid_number = BaseClass::longerlipid_number;
    this->n_fraction = BaseClass::n_fraction;
    this->MeanLatticeOccupany = BaseClass::MeanLatticeOccupany;
    
    execute();
};

template<class IngredientsType>
bool UpdaterLipidsRandomInitializer<IngredientsType>::execute()
{
    //Move to add monomers
     MoveAddMonomerSc<int32_t> move;
     
    //Initialize the lipids at random position and random orientation.
    // and then create solvent.
     randomInitializer(move);
     solventCreator(move, rI=true);
};

/**
* @brief Creates lipids randomly everywher in the box. 
* It creates n_fraction of n1 lipids.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template<class IngredientsType>
void UpdaterLipidsRandomInitializer<IngredientsType>::randomInitializer(MoveAddMonomerSc<int32_t>& move)
{
    //Initialize the variables to calculate number of lipids number
    //of lipids type of lipid. Also Instantiate an object of simple simulater 
    //as the it might be required to move the system. 
    
    int32_t typeofLipid;
    SimpSimu sim(ingredients,1000);
    int32_t shorterlipidcounter=0,longerlipidcounter=0;
    bool OneofTypeComplete=false;
    VectorInt3 InitialPos;
   
    for(int32_t lipidnumber =0;lipidnumber<nlipids;lipidnumber++){
                
        //Randomly choose a type until one of them is completely created. 
         
        if(OneofTypeComplete==false)
            typeofLipid=(rand()%2)?shorterlipid_length:longerlipid_length;

        //Create a lipid. 
    
        singleLipidCreator(move,typeofLipid,sim);
        
        //Count the both type of lipids.
         
        if(typeofLipid==shorterlipid_length)
            shorterlipidcounter++;
        else
            longerlipidcounter++;  
        
        //Tells the system that creation one of the type of lipids is finished 
         
        if(shorterlipidcounter==(shorterlipid_number)||longerlipidcounter==(longerlipid_number)){
            OneofTypeComplete=true;
            typeofLipid=(shorterlipidcounter==shorterlipid_number)?longerlipid_length:shorterlipid_length;
            }
        }        
};

/**
* @brief Creates a single lipid given total monomer in the lipid. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
* 
* @param move MoveAddMonomerSc type move used to add monomers of lipid.
* @param totalMonLipids total number of monomers in lipid to be created .
* @param sim Simple simulater to move the system if no possible bond vector
* can be used to add a lipid monomer. 
*/

template<class IngredientsType>
void UpdaterLipidsRandomInitializer<IngredientsType>::singleLipidCreator(MoveAddMonomerSc<int32_t>& move,int32_t totalMonLipids,SimpSimu& sim){

     //Initialize the variables to calculate rejected trials to create monomers
     //and other variables to be used later. 
      
    int32_t rejectedMonoTrials,x,y,z,type,lastMonomerAdded,monNum;
    VectorInt3 position;
    
    int32_t tail_length=(totalMonLipids-3)/2;
    int32_t type_chr=(totalMonLipids==shorterlipid_length)?1:3;
    
    VectorInt3 bondVector(0,0,0);
    
    //Random initial position.    
    x=rand()%(ingredients.getBoxX());
    y=rand()%(ingredients.getBoxY());
    z=rand()%(ingredients.getBoxZ());
    position.setAllCoordinates(x,y,z);
    rejectedMonoTrials=0;
    monNum=0;
    while(monNum<totalMonLipids){
        
        type=(monNum<3)? type_chr : type_chr+1;
        
        move.init(ingredients);
        move.setTag(type);
        move.setPosition(position+bondVector);
        //Check if position is available otherwise
        //try a new bond vector.        
        
        if(move.check(ingredients)==true){
            
            move.apply(ingredients);
            rejectedMonoTrials=0;
            position+=bondVector;
            
            lastMonomerAdded = ingredients.getMolecules().size()-1;
            
            if(monNum >0 && monNum!=(totalMonLipids-tail_length))
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-1);
            
            if(monNum==totalMonLipids-tail_length)
                ingredients.modifyMolecules().connect(lastMonomerAdded,lastMonomerAdded-tail_length-1);
                
            
            if(monNum==totalMonLipids-tail_length-1)
                position=ingredients.getMolecules()[lastMonomerAdded-tail_length];
                
            monNum++;
            }
        
        else{
            rejectedMonoTrials++;
            
             //If on the average all the bond vectors are used,
             //move the system by simply simulating for 1000 MCS.
            
            if(rejectedMonoTrials==108){
                sim.initialize();
                sim.execute();
                lastMonomerAdded = ingredients.getMolecules().size()-1;
                position=ingredients.getMolecules()[lastMonomerAdded];
                
                if(monNum==totalMonLipids-tail_length-1)
                    position=ingredients.getMolecules()[lastMonomerAdded-tail_length];
                
                rejectedMonoTrials=0;
            }
        }
        
        bondVector=drawRandomDirection();
        } 
};

/**
* @brief Returns a random vector from the set of bond vectorss. 
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template<class IngredientsType>
VectorInt3 UpdaterLipidsRandomInitializer<IngredientsType>::drawRandomDirection() const
{
    VectorInt3 bondVector;
    size_t randVectorId=rand()%(ingredients.getBondset().size());
    randVectorId=randVectorId%46;
    randVectorId += (randVectorId<17) ? 17 : 0;
    bondVector=ingredients.getBondset().getBondVector(randVectorId);
    return bondVector;
};
 
