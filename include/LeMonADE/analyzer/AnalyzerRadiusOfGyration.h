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

#ifndef LEMONADE_ANALYZER_RADIUS_OF_GYRATION_H
#define LEMONADE_ANALYZER_RADIUS_OF_GYRATION_H

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
 * @brief Analyzer for evaluating the squared radius of gyration Rg^2
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details This analyzer calculates the squared radius of gyration of either the
 * complete system (if no subgroup of monomers was specified), or the average of
 * the Rg^2 of a set of subgroups (if subgroups were specified using the provided
 * function void setMonomerGroups(std::vector<MonomerGroup<molecules_type> >) )
 * The analyzer calculates Rg^2 in every timestep and saves the time series
 * of this data to disk in the format
 * mcs <Rg2_x> <Rg2_y> <Rg2_z>  <Rg2_tot>
 * The default output filename is Rg2TimeSeries.dat and it can be changed by argument
 * to the constructor or by using the setter function provided.
 * If a more sophisticated grouping of monomers into groups is required, one can
 * also write a new analyzer, inheriting from AnalyzerRadiusOfGyration, and
 * overwriting the initialize function.
 */
template < class IngredientsType > class AnalyzerRadiusOfGyration : public AbstractAnalyzer
{

private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;
	//! reference to the complete system
	const IngredientsType& ingredients;
	//! Rg2 is calculated for the groups in this vector
	std::vector<MonomerGroup<molecules_type> > groups;
	//! timeseries of the Rg^2. Components: [0]-> Rg^2_x, [1]->Rg^2_y, [2]->Rg^2_z, [3]:Rg^2_tot
	std::vector< std::vector<double> > Rg2TimeSeries;
	//! vector of mcs times for writing the time series
	std::vector<double> MCSTimes;
	//! max length of internal buffer for each coordinate, before saving to disk
	uint32_t bufferSize;
	//! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
	std::string outputFile;
	//! flag used in dumping time series output
	bool isFirstFileDump;
	//! save the current values in Rg2TimeSeriesX, etc., to disk
	void dumpTimeSeries();
	//! calculate the Rg squared of the monomer group
	VectorDouble3 calculateRg2Components(const MonomerGroup<molecules_type>& group) const;
protected:
	//! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
	void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

	//! constructor
	AnalyzerRadiusOfGyration(const IngredientsType& ing,
			  std::string filename="Rg2TimeSeries.dat");

	//! destructor. does nothing
	virtual ~AnalyzerRadiusOfGyration(){}
	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();
	//! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();
	//! Set the number of values, after which the time series is saved to disk
	void setBufferSize(uint32_t size){bufferSize=size;}
	//! Change the output file name
	void setOutputFile(std::string filename){outputFile=filename;isFirstFileDump=true;}

};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * @param filename output file name. defaults to "Rg2TimeSeries.dat".
 * */
template<class IngredientsType>
AnalyzerRadiusOfGyration<IngredientsType>::AnalyzerRadiusOfGyration(
	const IngredientsType& ing,
	std::string filename)
:ingredients(ing)
,Rg2TimeSeries(4,std::vector<double>(0))
,bufferSize(100)
,outputFile(filename)
,isFirstFileDump(true)
{
}


/**
 * @details fills all monomers into the groups vector as a single group,
 * if the group size is zero. Normally this should always be the case,
 * unless the groups were already set explicitly by using the setter function
 * setMonomerGroups
 * */
template< class IngredientsType >
void AnalyzerRadiusOfGyration<IngredientsType>::initialize()
{
	//if no groups are set, use the complete system by default
	//groups can be set using the provided access function
	if(groups.size()==0){
		groups.push_back(MonomerGroup<molecules_type>(ingredients.getMolecules()));
		for(size_t n=0;n<ingredients.getMolecules().size();n++)
			groups[0].push_back(n);
	}

}

/**
 * @details Calculates the current Rg2, saves it in the
 * time series, and saves the time series to disk in regular intervals.
 * */
template< class IngredientsType >
bool AnalyzerRadiusOfGyration<IngredientsType>::execute()
{
	VectorDouble3 Rg2Components(0.0,0.0,0.0);

	for(size_t n=0;n<groups.size();n++)
	{
		//this vector will contain (Rg^2_x, Rg^2_y, Rg^2_z), i.e. the squared components!
		Rg2Components+=calculateRg2Components(groups[n])/double(groups.size());
	}

	Rg2TimeSeries[0].push_back(Rg2Components.getX());
	Rg2TimeSeries[1].push_back(Rg2Components.getY());
	Rg2TimeSeries[2].push_back(Rg2Components.getZ());
	Rg2TimeSeries[3].push_back(Rg2Components.getX()+Rg2Components.getY()+Rg2Components.getZ());
	MCSTimes.push_back(ingredients.getMolecules().getAge());
	//save to disk in regular intervals
	if(MCSTimes.size()>=bufferSize)
		dumpTimeSeries();

	return true;
}


template<class IngredientsType>
void AnalyzerRadiusOfGyration<IngredientsType>::cleanup()
{
	std::cout<<"AnalyzerRadiusOfGyration::cleanup()...";
	//write the remaining data from the time series
	dumpTimeSeries();
	std::cout<<"done\n";
}

/**
 * @details Saves the current content of Rg2TimeSeries to the file
 * Rg2TimeSeries.dat. The output format is:
 * mcs Rg2_x Rg2_y Rg2_z Rg2_tot
 * */
template<class IngredientsType>
void AnalyzerRadiusOfGyration<IngredientsType>::dumpTimeSeries()
{
	//fist make a single vector<vector<double> > for writing the results
	std::vector<std::vector<double> > resultsTimeseries=Rg2TimeSeries;
	resultsTimeseries.insert(resultsTimeseries.begin(),MCSTimes);

	//if it is written for the first time, include comment in the output file
	if(isFirstFileDump){
		std::stringstream commentTimeSeries;
		commentTimeSeries<<"Created by AnalyzerRadiusOfGyration\n";
		commentTimeSeries<<"file contains time series of average Rg_squared (Rg2) over all analyzed groups\n";
		commentTimeSeries<<"format: mcs\t Rg2X\t Rg2Y\t Rg2Z\t Rg2Total\n";

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
	//set all time series vectors back to zero size
	MCSTimes.resize(0);
	Rg2TimeSeries.resize(0);
	Rg2TimeSeries.resize(4,std::vector<double>(0));
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


