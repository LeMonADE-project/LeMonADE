/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project: topological effects 
 *****************************************************************************/

#ifndef LEMONADE_ANALYZER_ABSTRACT_DUMP_H
#define LEMONADE_ANALYZER_ABSTRACT_DUMP_H

#include <string>
#include <vector>
#include <sstream>

#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/*************************************************************************
 * definition of AnalyzerAbstractDump class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerAbstractDump
 *
 * @brief Abstract analyzer which provides a dump function to write out data after a certain buffer is reached
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details Only works if the number of measurements per time step stays constant. 
 */
template < class IngredientsType , class T > class AnalyzerAbstractDump : public AbstractAnalyzer
{

private:
	//! analyze after equilibration time
	uint32_t equilibrationTime;
	//! name of the output file
	std::string outputFilename;
	/**/
	uint32_t NColumns;

protected:
	//! reference to the complete system
	const IngredientsType& ingredients;
	//! vector of mcs times for writing the time series
	std::vector<T> MCSTimes;
	//! time series of the data 
	std::vector< std::vector<T> > Data;
public:

	//! constructor
	AnalyzerAbstractDump(const IngredientsType& ingredients_, std::string outputFilename_, uint32_t equilibrationTime_=0);

	//! destructor. does nothing
	virtual ~AnalyzerAbstractDump(){}
	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();
	//! Called by TaskManager::execute()
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();
private:
	/*stringstream to be written out to the file*/
	std::string commentTimeSeries;
	/*threeshold for executing the function dumpTimeSeries()*/
	uint32_t bufferSize;
	/*the header for the output file has to be written only once*/
	bool isFirstFileDump;
	/*control bool for the last writing*/
	bool isCleanup;
protected:
	//!dump data and cleans buffer 
	void dumpTimeSeries();
private:
	/*resize the data container*/
	inline void resizeData(){ Data.resize(0); Data.resize(NColumns,std::vector<T>(0)); }
public:
	//! Set equilibration time 
	void setEquilibrationTime(uint32_t time){equilibrationTime=time;}
	//! Get equilibration time 
	uint32_t getEquilibrationTime(){return equilibrationTime;} 
	//! set output filename
	void setOutputFilename(std::string outputFilename_){outputFilename=outputFilename_;}
	//! get output filename
	std::string getOutputFilename(){return outputFilename;}
	//! set the comment stringstream
	void setComment(std::string commentTimeSeries_){commentTimeSeries=commentTimeSeries_;}
	//! get the comment stringstream
// 	std::stringstream getComment(){return commentTimeSeries;}
	//! set the number of columns 
	void setNumberOfColumns(uint32_t NColumns_){NColumns=NColumns_;resizeData();}
	//! get the number of colums
	uint32_t getNumberOfColumns(){ return NColumns; }

};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * */
template<class IngredientsType, class T >
AnalyzerAbstractDump<IngredientsType, T >::AnalyzerAbstractDump(
	const IngredientsType& ingredients_, const std::string outputFilename_, const uint32_t equilibrationTime_ )
:ingredients(ingredients_),
outputFilename(outputFilename_),
equilibrationTime(equilibrationTime_),
bufferSize(100),
NColumns(0),
Data(NColumns,std::vector<T>(0)),
isFirstFileDump(true),
isCleanup(false)
{
  std::stringstream comment; 
  comment<<"Created by AnalyzerAbstractDump\n";
  comment<<"format: mcs \t data \n";
  commentTimeSeries=comment.str();
}


/**
 * @details 
 * */
template< class IngredientsType, class T >
void AnalyzerAbstractDump<IngredientsType, T >::initialize()
{
}

/**
 * @details 
 * */
template< class IngredientsType, class T >
bool AnalyzerAbstractDump<IngredientsType, T >::execute()
{
  
  MCSTimes.push_back( ingredients.getMolecules().getAge() );
  dumpTimeSeries();
  return true;
}


template<class IngredientsType, class T >
void AnalyzerAbstractDump<IngredientsType, T >::cleanup()
{ 
  isCleanup=true;
  dumpTimeSeries();
}
/*****************************************************************************/
/**
 * @details Saves the current content of DistanceTimeSeries to the file outputfile
 * . The output format is:
 * mcs 
 * */
template<class IngredientsType, class T >
void AnalyzerAbstractDump<IngredientsType, T >::dumpTimeSeries()
{
    if (Data[0].size() > bufferSize || isCleanup ){
	if ( NColumns != Data.size() || NColumns == 0 )
	    throw std::runtime_error("[AnalyzerAbstractDump]::dumpTimeSeries() Did not set up NColumns!");
	for (size_t i = 0; i < Data.size(); ++i)
	  std::cout << i << " " << Data[i].size() << " " << MCSTimes.size() << "\n"; 
	//fist make a single vector<vector<T> > for writing the results
	std::vector<std::vector<T> > resultsTimeseries(Data);
	resultsTimeseries.insert(resultsTimeseries.begin(),MCSTimes);
	//if it is written for the first time, include comment in the output file
	if(isFirstFileDump){
		ResultFormattingTools::writeResultFile(
			outputFilename,
			ingredients,
			resultsTimeseries,
			commentTimeSeries
		);

		isFirstFileDump=false;
	}
	//otherwise just append the new data
	else{
		ResultFormattingTools::appendToResultFile(outputFilename,
							  resultsTimeseries);
	}
	//set all time series vectors back to zero size
	MCSTimes.resize(0);
	resizeData();
    }
}

#endif /*LEMONADE_ANALYZER_ABSTRACT_DUMP_H*/


