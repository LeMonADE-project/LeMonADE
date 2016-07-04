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

#ifndef LEMONADE_ANALYZER_RG_ANALYZER_H
#define LEMONADE_ANALYZER_RG_ANALYZER_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>

/*************************************************************************
 * definition of AnalyzerRadiusOfGyration class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerRadiusOfGyration
 *
 * @brief Analyzer for evaluating the radius of gyration RG2
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * 
 * @details For the MonomerGroups provided as argument to the constructor, this 
 * Analyzer calculates the average Radius of gyration. It produces two output files:
 * - prefix_averages.dat: Contains the averages Rg2_x, Rg2_y, Rg2_z,Rg2_tot for all groups
 * - prefix_timeseries.dat: Contains the time series of Rg2_x, Rg2_y, Rg2_z,Rg2_tot for all groups
 * .
 * "prefix" in the filenames stands for a string which can be set in the constructor
 * or using the access function setOutputFilePrefix(string)m, and defaults to "Rg2".
 *
 */
template < class IngredientsType > class AnalyzerRadiusOfGyration : public AbstractAnalyzer
{
	
private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;
	//! reference to the complete system
	const IngredientsType& ingredients;
	//! Rg2 is calculated for the groups in this vector
	const std::vector<MonomerGroup<molecules_type> >& groups;
	//! these vectors keep track of the Rg^2 averages of all groups to be analyzed
	std::vector<double> Rg2SumX,Rg2SumY,Rg2SumZ,Rg2SumTotal;
	//! timeseries of the Rg^2 of all groups. Values are saved to disk in intervals of saveInterval 
	std::vector< std::vector<double> > Rg2TimeSeriesX,Rg2TimeSeriesY,Rg2TimeSeriesZ,Rg2TimeSeriesTotal;
	//! vector of mcs times for writing the time series
	std::vector<double> MCSTimes;
	//! number of samples for the averages Rg2SumX, etc.
	uint32_t nValues;
	//! mcs below this number are not considered in averages
	uint64_t lowerLimitForAnalysis;
	//! interval for saving Rg2TimeSeries to file
	uint32_t saveInterval;
	//! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
	std::string outputFilePrefix;
	//! flag used in dumping time series output
	bool timeSeriesFirstDump;
	//! save the current values in Rg2TimeSeriesX, etc., to disk
	void dumpTimeSeries();
	//! calculate the Rg squared of the monomer group
	VectorDouble3 calculateRg2Components(const MonomerGroup<molecules_type>& group) const;
	
public:
  
	//! constructor
	AnalyzerRadiusOfGyration(const IngredientsType& ing,
				 const std::vector<MonomerGroup<molecules_type> >& groupVector,
			  std::string filenamePrefix=std::string("Rg2"),
			  uint64_t lowerMcsLimit=0);
	
	//! destructor. does nothing
	virtual ~AnalyzerRadiusOfGyration(){}
	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();
	//! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();
	//! Set the prefix for the output files
	void setOutputFilePrefix(std::string prefix){outputFilePrefix=prefix;}
	//! Set the number of values, after which the time series is saved to disk
	void setSaveInterval(uint32_t interval){saveInterval=interval;}
	
};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * @param groupVector vector of MonomerGroup to be analyzed
 * @param filenamePrefix prefix for the two output files. defaults to "Rg2".
 * @param lowerMcsLimit mcs below this are ignored for evaluation of averages
 * */
template<class IngredientsType>
AnalyzerRadiusOfGyration<IngredientsType>::AnalyzerRadiusOfGyration(
	const IngredientsType& ing, 
	const std::vector< MonomerGroup< molecules_type > >& groupVector,
	std::string filenamePrefix,
	uint64_t lowerMcsLimit)
:ingredients(ing)
,groups(groupVector)
,lowerLimitForAnalysis(lowerMcsLimit)
,saveInterval(100)
,outputFilePrefix(filenamePrefix)
,timeSeriesFirstDump(true)
,nValues(0)
{
}

/**
 * @throw std::runtime_error if there are no monomers in the groups
 * */
template< class IngredientsType >
void AnalyzerRadiusOfGyration<IngredientsType>::initialize()
{
	//test if the groups contain monomers and exit otherwise
	if(groups.size()==0)
		throw std::runtime_error("AnalyzerRadiusOfGyration::initialize(): no monomers to analyze");
	//give a warning if a single group is empty
	for(size_t n=0;n<groups.size();n++){
		if(groups[n].size()==0){
			std::cerr<<"AnalyzerRadiusOfGyration::initialize():\n";
			std::cerr<<"Warning: group "<<n<<" of "<<groups.size()<<" is empty\n";
		}
	}
	
	//resize the vectors to the number of groups
	Rg2SumX.resize(groups.size(),0.0);
	Rg2SumY.resize(groups.size(),0.0);
	Rg2SumZ.resize(groups.size(),0.0);
	Rg2SumTotal.resize(groups.size(),0.0);
	Rg2TimeSeriesX.resize(groups.size());
	Rg2TimeSeriesY.resize(groups.size());
	Rg2TimeSeriesZ.resize(groups.size());
	Rg2TimeSeriesTotal.resize(groups.size());
}

/**
 * @details Calculates the current Rg2, adds it to the average, saves it in the 
 * time series, and saves the time series to disk in regular intervals. 
 * */
template< class IngredientsType >
bool AnalyzerRadiusOfGyration<IngredientsType>::execute()
{
	VectorDouble3 Rg2Components;
	
	for(size_t n=0;n<groups.size();n++)
	{
		//this vector will contain (Rg^2_x, Rg^2_y, Rg^2_z), i.e. the squared components!
		Rg2Components=calculateRg2Components(groups[n]);
		//add to averages only above certain mcs
		if (ingredients.getMolecules().getAge()>=lowerLimitForAnalysis){
			//add values to groupwise averages
			Rg2SumX[n]+=Rg2Components.getX();
			Rg2SumY[n]+=Rg2Components.getY();
			Rg2SumZ[n]+=Rg2Components.getZ();
			Rg2SumTotal[n]+=(Rg2Components.getX()+Rg2Components.getY()+Rg2Components.getZ());
			nValues++;
		}
		//save values in time series
		Rg2TimeSeriesX[n].push_back(Rg2Components.getX());
		Rg2TimeSeriesY[n].push_back(Rg2Components.getY());
		Rg2TimeSeriesZ[n].push_back(Rg2Components.getZ());
		Rg2TimeSeriesTotal[n].push_back(Rg2Components.getX()+Rg2Components.getY()+Rg2Components.getZ());	
		
	}
	MCSTimes.push_back(ingredients.getMolecules().getAge());
	//save to disk in regular intervals
	if(MCSTimes.size()>=saveInterval)
		dumpTimeSeries();
	
	return true;
}


template<class IngredientsType>
void AnalyzerRadiusOfGyration<IngredientsType>::cleanup()
{
	std::cout<<"AnalyzerRadiusOfGyration::cleanup()\n";
	std::cout<<"printing radius of gyration results to file...";
	
	//first write the remaining data from the time series
	dumpTimeSeries();
	
	//now write the averages into a separate file
	
	//open a file
	std::ofstream fileAverages;
	std::string fname(outputFilePrefix);
	fname+=std::string("_averages.dat");
	fileAverages.open(fname.c_str());
	
	if(!fileAverages.is_open())
		throw std::runtime_error("AnalyzerRadiusOfGyration::cleanup(): error opening output file");
	//add comments and meta-data to file
	ingredients.printMetaData(fileAverages);
	fileAverages<<"# Created by AnalyzerRadiusOfGyration\n";
	fileAverages<<"# file contains average Rg_squared (Rg2) for every group G1...Gn\n";
	fileAverages<<"# format: Rg2X_G1\t Rg2Y_G1\t Rg2Z_G1\t Rg2Total_G1\t ...";
	fileAverages<<"Rg2X_Gn\t Rg2Y_Gn\t Rg2Z_Gn\t Rg2Total_Gn\n";
	//save the values
	for(size_t n=0;n<groups.size();n++)
	{
		fileAverages<<Rg2SumX[n]/double(nValues)<<"\t";
		fileAverages<<Rg2SumY[n]/double(nValues)<<"\t";
		fileAverages<<Rg2SumZ[n]/double(nValues)<<"\t";
		fileAverages<<Rg2SumTotal[n]/double(nValues)<<"\t";
	}
	
	fileAverages.close();
	std::cout<<"done\n";
}

/**
 * @details Saves the current content of the time series of the Rg2 of all 
 * groups to the file filenamePrefix_timeseries.dat. The output format
 * is: mcs Rg2_x_g1 Rg2_y_g1 Rg2_z_g1 Rg2_tot_g1 ...Rg2_x_gn Rg2_y_gn Rg2_z_gn Rg2_tot_gn 
 * where the index g1...gn refers to the n groups in the vector groups.
 * */
template<class IngredientsType>
void AnalyzerRadiusOfGyration<IngredientsType>::dumpTimeSeries()
{
	//fist make a single vector<vector<double> > for writing the results
	std::vector<std::vector<double> > resultsTimeseries;
	resultsTimeseries.push_back(MCSTimes);
	for(size_t n=0;n<groups.size();n++)
	{
		resultsTimeseries.push_back(Rg2TimeSeriesX[n]);
		resultsTimeseries.push_back(Rg2TimeSeriesY[n]);
		resultsTimeseries.push_back(Rg2TimeSeriesZ[n]);
		resultsTimeseries.push_back(Rg2TimeSeriesTotal[n]);		
	}
	
	//if it is written for the first time, include comment in the output file
	if(timeSeriesFirstDump){
		std::stringstream commentTimeSeries;
		commentTimeSeries<<"Created by AnalyzerRadiusOfGyration\n";
		commentTimeSeries<<"file contains time series of Rg_squared (Rg2) for every group G1...Gn\n";
		commentTimeSeries<<"format: mcs\t Rg2X_G1\t Rg2Y_G1\t Rg2Z_G1\t Rg2Total_G1\t";
		commentTimeSeries<<"Rg2X_Gn\t Rg2Y_Gn\t Rg2Z_Gn\t Rg2Total_Gn\n";
	
		ResultFormattingTools::writeResultFile(
			outputFilePrefix+std::string("_timeseries.dat"),
						       ingredients,
					 resultsTimeseries,
					 commentTimeSeries.str());
		
		timeSeriesFirstDump=false;
	}
	//otherwise just append the new data
	else{
		ResultFormattingTools::appendToResultFile(
			outputFilePrefix+std::string("_timeseries.dat"),
					 resultsTimeseries);
	}
	//set all time series vectors back to zero size
	MCSTimes.resize(0);
	for(size_t n=0;n<groups.size();n++)
	{
		Rg2TimeSeriesX[n].resize(0);
		Rg2TimeSeriesY[n].resize(0);
		Rg2TimeSeriesZ[n].resize(0);
		Rg2TimeSeriesTotal[n].resize(0);
	}	
}

/**
 * @details calculates the three components Rg^2_x, Rg^2_y,Rg^2_z and returns 
 * them in a vector.
 * @return VectorDouble3 containing the components Rg^2_x, Rg^2_y,Rg^2_z, or (0.0,0.0,0.0) if group is empty)
 * @param group the monomer group of which the Rg2 is calculated
 * */
template<class IngredientsType>
VectorDouble3 AnalyzerRadiusOfGyration<IngredientsType>::calculateRg2Components(
	const MonomerGroup<molecules_type>& group) const
{
	//if group is empty, return zero vector
	if(group.size()==0){
		return VectorDouble3(0.0,0.0,0.0);
	}
	
	
	VectorDouble3 sum_sqr; VectorDouble3 CoM_sum;
	//first calculate the center of mass
	for ( size_t n = 0; n < group.size(); ++n)
	{
		CoM_sum.setX( CoM_sum.getX() + group[n].getX() );
		CoM_sum.setY( CoM_sum.getY() + group[n].getY() );
		CoM_sum.setZ( CoM_sum.getZ() + group[n].getZ() );
	}
	double inv_N = 1.0 / double ( group.size() );
	
	VectorDouble3 CoM (double ( CoM_sum.getX() ) * inv_N, 
			   double ( CoM_sum.getY() ) * inv_N,
			   double ( CoM_sum.getZ() ) * inv_N);
	
	//now calculate the Rg2 using the center of mass
	for ( uint32_t n = 0; n < group.size(); ++n)
	{
		double diffX,diffY,diffZ;
		diffX = double(group[n].getX()) - CoM.getX();	
		diffY = double(group[n].getY()) - CoM.getY();
		diffZ = double(group[n].getZ()) - CoM.getZ();
		sum_sqr +=VectorDouble3(diffX*diffX,diffY*diffY,diffZ*diffZ);
	}
	return sum_sqr / double ( group.size() );
		
}
#endif


