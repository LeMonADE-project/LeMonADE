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

#include <LeMonADE/LeMonADE.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/NumericTools.h>
#include <LeMonADE/utility/NumericFunctions.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/Vector3D.h>


/************************************************************************
 * this analyzer calculates the average Rg^2 of groups of particles.
 * it produces two ouput files:
 * - RgSquaredAverages: the average Rg^2 (component wise) for every 
 *   group given to analyze (separately)
 * -RgSquaredTimeseries: time series of the total Rg^2 (componentwise)
 *  of all monomers given to analyze. If a vector of monomer groups is
 *  handed to the constructor (e.g. to define polymer chains), then the 
 *  time series is the Rg^2 averaged over all these groups (i.e. the average
 *  chain Rg^2) 
 * **********************************************************************/

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
 *
 *
 * @deprecated
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 *
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test!!!
 */
template < class IngredientsType > class AnalyzerRadiusOfGyration : public AbstractAnalyzer
{
	typedef typename IngredientsType::molecules_type molecules_type;
	
	const IngredientsType& ingredients;
	
	std::vector<MonomerGroup<molecules_type> > groups;
	
	TimeSeries<VectorDouble3>  Rg2TimeSeries;
	std::vector<Moments<double> > Rg_squaredX,Rg_squaredY,Rg_squaredZ,Rg_squared;
	
	//! mcs below this number are not considered in analysis
	uint64_t lowerLimit;
	
	//! prefix for output filename
	std::string outputFilePrefix;
public:
  
	/*
	 * there are three different constructors: 
	 * - the first calculates the Rg^2 of the complete system
	 * - the second takes a monomer group as additional argument and 
	 *   calculates the average Rg^2 of the monomers in this group
	 * - the third takes a std::vector of monomer groups and calculates
	 *   the Rg^2 of all molecules in this group. The time series output is
	 *   then the average Rg^2, averaged over both time and groups.
	 * */
	AnalyzerRadiusOfGyration(const IngredientsType& ing,std::string filenamePrefix="",uint64_t lowerMcsLimit=0);
	AnalyzerRadiusOfGyration(const IngredientsType& ing,
				 const MonomerGroup<molecules_type>& g,
			  std::string filenamePrefix="",
			  uint64_t lowerMcsLimit=0);
	AnalyzerRadiusOfGyration(const IngredientsType& ing,
				 const std::vector<MonomerGroup<molecules_type> >& g,
			  std::string filenamePrefix="",
			  uint64_t lowerMcsLimit=0);
	
	virtual ~AnalyzerRadiusOfGyration(){}
	
	virtual bool execute();
	virtual void cleanup();
	
	  /**
	   * @brief This function is called \a once in the beginning of the TaskManager.
	   *
	   * @details It´s a virtual function for inheritance.
	   * Use this function for initializing tasks (e.g. init SDL)
	   *
	   **/
	  virtual void initialize(){};


};

/*************************************************************************
 * implementation of memeber execute()
 * ***********************************************************************/

template<class IngredientsType>
AnalyzerRadiusOfGyration<IngredientsType>::AnalyzerRadiusOfGyration(const IngredientsType& ing,
								    std::string filenamePrefix,
								    uint64_t lowerMcsLimit)
:ingredients(ing)
,Rg_squaredX(1)
,Rg_squaredY(1)
,Rg_squaredZ(1)
,Rg_squared(1)
,lowerLimit(lowerMcsLimit)
,outputFilePrefix(filenamePrefix)
{
	groups.push_back(MonomerGroup<molecules_type>(&(ing.getMolecules())));
	
	//add all monomers to the group
	for(size_t n=0; n<ing.getMolecules().size();++n)
	{
		groups[0].push_back(n);
	}
}

template<class IngredientsType>
AnalyzerRadiusOfGyration<IngredientsType>::AnalyzerRadiusOfGyration(const IngredientsType& ing, 
								    const MonomerGroup< molecules_type >& g
								    ,std::string filenamePrefix
								    ,uint64_t lowerMcsLimit)
:ingredients(ing)
,Rg_squaredX(1)
,Rg_squaredY(1)
,Rg_squaredZ(1)
,Rg_squared(1)
,lowerLimit(lowerMcsLimit)
,outputFilePrefix(filenamePrefix)
{
	//add the group given as argument to groups. groups will only consist of this single group
	groups.push_back(g);
}


template<class IngredientsType>
AnalyzerRadiusOfGyration<IngredientsType>::AnalyzerRadiusOfGyration(const IngredientsType& ing, 
								    const std::vector< MonomerGroup< molecules_type > >& g,
								    std::string filenamePrefix,
								    uint64_t lowerMcsLimit)
:ingredients(ing)
,groups(g)
,Rg_squaredX(g.size())
,Rg_squaredY(g.size())
,Rg_squaredZ(g.size())
,Rg_squared(g.size())
,lowerLimit(lowerMcsLimit)
,outputFilePrefix(filenamePrefix)
{
}


//measure the squared Rg for all groups 
template< class IngredientsType >
bool AnalyzerRadiusOfGyration<IngredientsType>::execute()
{
	if (ingredients.getMolecules().getAge()<lowerLimit) return true;
	
	VectorDouble3 Rg2Components;
	VectorDouble3 Rg2ComponentsSum;
	for(size_t n=0;n<groups.size();n++)
	{
		//this vector contains (Rg^2_x, Rg^2_y, Rg^2_z), i.e. the squared components!
		Rg2Components=squaredRadiusOfGyrationComponents(groups[n]);
		//save into time series
		Rg2ComponentsSum+=Rg2Components;
		
		//add values to groupwise averages
		Rg_squaredX[n].add(Rg2Components.getX());
		Rg_squaredY[n].add(Rg2Components.getY());
		Rg_squaredZ[n].add(Rg2Components.getZ());
		Rg_squared[n].add(Rg2Components.getX()+Rg2Components.getY()+Rg2Components.getZ());
		
	}
	Rg2ComponentsSum/=double(groups.size());
	Rg2TimeSeries.add(ingredients.getMolecules().getAge(),Rg2ComponentsSum);
	return true;
}

template<class IngredientsType>
void AnalyzerRadiusOfGyration<IngredientsType>::cleanup()
{
	std::cout<<"printing radius of gyration results to file...";
	//first print the time series
	std::vector<std::vector<double> > resultsTimeseries;
	resultsTimeseries.resize(5);
	//copy results of first group including times
	std::map<uint64_t,VectorDouble3> groupResults;
	std::map<uint64_t,VectorDouble3>::const_iterator it;
	
	groupResults=Rg2TimeSeries.get();
	for(it=groupResults.begin();it!=groupResults.end();++it)
	{
		resultsTimeseries[0].push_back(double(it->first));
		resultsTimeseries[1].push_back(it->second.getX());
		resultsTimeseries[2].push_back(it->second.getY());
		resultsTimeseries[3].push_back(it->second.getZ());
		resultsTimeseries[4].push_back(it->second.getX()+it->second.getY()+it->second.getZ());
	}
	
	
	//now write to file (including a comment)
	std::stringstream filenameTimeSeries;
	filenameTimeSeries<<outputFilePrefix<<"_RgSquaredTimeseries.dat";
	
	std::stringstream commentTimeSeries;
	commentTimeSeries<<"Created by RgAnalyzer\n";
	commentTimeSeries<<"sample size: "<<groups.size()<<" groups\n";
	
	commentTimeSeries<<"format: mcs\t <Rg^2_x>\t <Rg^2_y>\t <Rg^2_z> <Rg^2>\n";
	
	
	ResultFormattingTools::writeResultFile(filenameTimeSeries.str(),ingredients,resultsTimeseries,commentTimeSeries.str());
	
	///////////////////////////////////////////////////////////////////////
	//now write the averages into a separate file
	std::vector<std::vector<double> > averages;
	averages.resize(5);
	for(size_t n=0;n<groups.size();n++)
	{
		averages[0].push_back(groups[n].size());
		averages[1].push_back(Rg_squaredX[n].m_1());
		averages[2].push_back(Rg_squaredY[n].m_1());
		averages[3].push_back(Rg_squaredZ[n].m_1());
		averages[4].push_back(Rg_squared[n].m_1());
	}
	
	std::stringstream commentAverages;
	commentAverages<<"Created by RgAnalyzer\n";
	commentAverages<<"sample size: "<<groups.size()<<" groups\n";
	commentAverages<<"file contains average Rg_squared for every group\n";
	commentAverages<<"format: N(size of group in monomers)\t Rg_squaredX\t Rg_squaredY\t Rg_squaredZ\t Rg_squared_total\n";
	
	std::stringstream filenameAverages;
	filenameAverages<<outputFilePrefix<<"_RgSquaredAverages.dat";
	ResultFormattingTools::writeResultFile(filenameAverages.str(),ingredients,averages,commentAverages.str());
	
	std::cout<<"done\n";
	
	
}




#endif


