/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project: topological effects 
 *****************************************************************************/

#ifndef LEMONADE_ANALYZER_ABSTRACT_MSD_H
#define LEMONADE_ANALYZER_ABSTRACT_MSD_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>
/*************************************************************************
 * definition of AnalyzerAbstractMSD class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerAbstractMSD
 *
 * @brief Analyzer for evaluating ...
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details 
 */
template < class IngredientsType > class AnalyzerAbstractMSD : public AbstractAnalyzer
{

private:
    	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;
	//! Position is calculated for the groups in this vector
	std::vector<MonomerGroup<molecules_type> > groups;
	//! vector of mcs times for writing the time series
	std::vector<double> MCSTimes;
	//! sum of position and square position and the number of measurements 
	std::vector<std::vector<double> > Fluctuations;
	//! time series of the msd 
	std::vector< std::vector<VectorDouble3> > PositionTimeSeries;
	//! write out the results
	void writeResults();
	//! analyze after equilibration time
	uint32_t equilibrationTime;
	//! name of the output file
	std::string outputFilename;
	//! 
	bool firstTime; 
	//!
	uint32_t startMCS;
	//! calculate the center of mass of a monomer group 
	template <class groupType>
	VectorDouble3 COMGroup(groupType group)
	{
	  VectorDouble3 COM;
	  for(size_t i = 0; i < group.size(); i++ )
	      COM+=group[i];
	  return COM/(double)group.size();
	}
	
	//!add values to Fluctuations
	void AddFluctuations(uint32_t ID, VectorDouble3 vec)
	{
	  Fluctuations[0][ID]+=(vec*vec);
	  Fluctuations[1][ID]+=vec.getX();
	  Fluctuations[2][ID]+=vec.getY();
	  Fluctuations[3][ID]+=vec.getZ();
	  Fluctuations[4][ID]=((double)Fluctuations[0][ID])/((double)numberOfMeasurements)
	  -((double)(
	    Fluctuations[1][ID]*Fluctuations[1][ID]
	   +Fluctuations[2][ID]*Fluctuations[2][ID]
	   +Fluctuations[3][ID]*Fluctuations[3][ID]))/((double)(numberOfMeasurements*numberOfMeasurements));
	  
	}
	//!
	uint32_t numberOfMeasurements;
	//convinience function 
	void CollectPositions();
	
	bool SystemCOMIsReference;
	//! used to determine the reference position
	MonomerGroup<molecules_type>  ReferenceGroup;
	
protected:
	//! reference to the complete system
	const IngredientsType& ingredients;
  	//! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
	void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
	//! get the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
	std::vector<MonomerGroup<molecules_type> > getMonomerGroups( ){return groups;}
	
	void setSystemCOMAsReference(bool SystemCOMIsReference_){ SystemCOMIsReference=SystemCOMIsReference_;}
public:

	//! constructor
	AnalyzerAbstractMSD(const IngredientsType& ingredients_, uint32_t equilibrationTime_=0);

	//! destructor. does nothing
	virtual ~AnalyzerAbstractMSD(){}
	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();
	//! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();
	
	//! Set equilibration time 
	void setEquilibrationTime(uint32_t time){equilibrationTime=time;}
	//! Get equilibration time 
	uint32_t getEquilibrationTime(){return equilibrationTime;} 
	//! set output filename
	void setOutputFilename(std::string outputFilename_){outputFilename=outputFilename_;}
	//! get output filename
	std::string getOutputFilename(){return outputFilename;}

};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * */
template<class IngredientsType>
AnalyzerAbstractMSD<IngredientsType>::AnalyzerAbstractMSD(
	const IngredientsType& ingredients_, uint32_t equilibrationTime_)
:ingredients(ingredients_),
equilibrationTime(equilibrationTime_),
firstTime(true),
SystemCOMIsReference(false),
ReferenceGroup(MonomerGroup<molecules_type>(ingredients.getMolecules()))
{
}
template<class IngredientsType>
void AnalyzerAbstractMSD<IngredientsType>::CollectPositions()
{
	{
	    if( (ingredients.getMolecules().getAge()-startMCS) >= equilibrationTime)
	    {
	      VectorDouble3 ReferencePosition(0,0,0);
	      if (SystemCOMIsReference)
	      {
		  ReferencePosition=COMGroup(ReferenceGroup) ;
	      }
	      if(firstTime)
	      {
		Fluctuations.resize(0);
		Fluctuations.resize(5, std::vector<double>(groups.size(),0)); 
		numberOfMeasurements++;
		//calculate the COM of each group in groups 
		for(uint32_t i=0;i< groups.size();i++)		
		{
		  VectorDouble3 COM(COMGroup(groups[i])-ReferencePosition);
		  std::vector<VectorDouble3> vec;
		  vec.push_back(COM);
		  PositionTimeSeries.push_back(vec);
		  AddFluctuations(i,COM);
		}
		firstTime=false;
	      }else
	      {
		numberOfMeasurements++;
		//calculate the COM of each group in groups     
		for(uint32_t i=0;i< groups.size();i++)		
		{
		  VectorDouble3 COM(COMGroup(groups[i])-ReferencePosition);
		  PositionTimeSeries[i].push_back(COM);
		  AddFluctuations(i,COM);
		}
	      }
	      MCSTimes.push_back(ingredients.getMolecules().getAge()-startMCS);
	    }
	}
}


/**
 * @details 
 * */
template< class IngredientsType >
void AnalyzerAbstractMSD<IngredientsType>::initialize()
{
  startMCS=ingredients.getMolecules().getAge();
  MonomerGroup<molecules_type> mygroup(ingredients.getMolecules());
  
  if(SystemCOMIsReference)
  {
//     RefernceGroup.push_back(MonomerGroup<typename IngredientsType::molecules_type>(ingredients.getMolecules()));
    for ( size_t i = 0; i < ingredients.getMolecules().size(); ++i )
	mygroup.push_back(i);
    ReferenceGroup = mygroup;
  }

  CollectPositions();
}

/**
 * @details 
 * */
template< class IngredientsType >
bool AnalyzerAbstractMSD<IngredientsType>::execute()
{
  //takes into account a certain equilibration time
  //after which the system is analyzed
  CollectPositions();
  return true;
}


template<class IngredientsType>
void AnalyzerAbstractMSD<IngredientsType>::cleanup()
{ 
  //calculate the mean square displacements for each group
    std::vector<std::vector<double> > MSD; 
    for (uint32_t i=0; i < PositionTimeSeries.size(); i++)
    {
      std::vector<double> MSDTimeSeriesPerGroup; 
      for (uint32_t Delta=0; Delta < PositionTimeSeries[i].size(); Delta++)
      {
	double MSDPerTimeDifference(0);
	for (uint32_t j=0; j < (PositionTimeSeries[i].size() - Delta) ; j++)
	{
	  VectorDouble3 Diff( PositionTimeSeries[i][j+Delta]-PositionTimeSeries[i][j]);
	  MSDPerTimeDifference+=(Diff.getX()*Diff.getX()+Diff.getY()*Diff.getY()+Diff.getZ()*Diff.getZ());
	}
	MSDPerTimeDifference/=((double)(PositionTimeSeries[i].size() - Delta));
	MSDTimeSeriesPerGroup.push_back(MSDPerTimeDifference);
      }
      MSD.push_back(MSDTimeSeriesPerGroup);
    } 

    //write out the time series MSD
    MSD.insert(MSD.begin(), MCSTimes);
    std::stringstream commentTimeSeriesMSD;
    commentTimeSeriesMSD<<"Created by AnalyzerAbstractMSD\n";
    commentTimeSeriesMSD<<"file contains time series of MSD of all groups\n";
    commentTimeSeriesMSD<<"format: Delta time  \t  Group_1 \t ... \t Group_N  \n";

    ResultFormattingTools::writeResultFile( outputFilename, ingredients, MSD, commentTimeSeriesMSD.str() );

    std::vector<std::vector<double> > AverageFluctuations;
    std::vector<double> GroupID;
    for (uint32_t i=0; i<groups.size();i++){GroupID.push_back(i+1);}
    AverageFluctuations.push_back(GroupID);
    AverageFluctuations.push_back(Fluctuations[4]);
    
    //write out the End Fluctuations
    
    std::stringstream commentTimeSeriesAverageFluctuations;
    commentTimeSeriesAverageFluctuations<<"Created by AnalyzerAbstractMSD\n";
    commentTimeSeriesAverageFluctuations<<"file contains average fluctuations for each group\n";
    commentTimeSeriesAverageFluctuations<<"IDs start at 1\n";
    commentTimeSeriesAverageFluctuations<<"format: Group ID  \t  <R2>-<R>2   \n";
    std::stringstream AverageFluctuationsFilename;
    AverageFluctuationsFilename << "AverageFluctuations" << outputFilename;
    ResultFormattingTools::writeResultFile( AverageFluctuationsFilename.str(), ingredients, AverageFluctuations, commentTimeSeriesAverageFluctuations.str() );
    
}


#endif /*LEMONADE_ANALYZER_ABSTRACT_MSD_H*/


