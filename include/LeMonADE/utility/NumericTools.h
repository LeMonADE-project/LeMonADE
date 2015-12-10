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

#ifndef LEMONADE_UTILITY_NUMERICTOOLS_H
#define LEMONADE_UTILITY_NUMERICTOOLS_H

/**
 * @file
 * @brief Some numeric tools that may be useful for analyzing data
 * */

#include <vector>


/**
 * @class Moments
 *
 * @brief calculates first and second momoments, as well as first and second central moments of scalar types
 *
 * @deprecated Untested!
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test!
 **/
template <class ScalarType> class Moments
{
public:
  Moments():sum(0),squareSum(0),nValues(0){};
  
  //! Add a new value to the collection
  void add(ScalarType val){sum+=val;squareSum+=val*val;++nValues;}

  //! Reset collected values to 0
  void reset(){sum=0;squareSum=0;nValues=0;}

  //! Merge another object of the same type into this one, such that the samples of both are combined
  void merge(Moments<ScalarType>& rhs){sum+=rhs.sum;squareSum+=rhs.squareSum;nValues+=rhs.nValues;}
  
  //! Zeroth moment(number of samples)
  uint32_t m_0() const {return nValues;}

  //! First moment (average)
  double m_1() const{ return double(sum)/double(nValues);}

  //! Second moment (average square)
  double m_2() const{return double(squareSum)/double(nValues);}

  //! First central moment (always 0)
  double mu_1() const{return 0.0;}

  //! Second central moment (variance)
  double mu_2() const{return m_2()-m_1()*m_1();}
 
private:
  ScalarType sum;
  ScalarType squareSum;
  uint32_t nValues;
};

/**
 * @class Histogram
 *
 * @brief Simple class for collecting a histogram of double values in integer bins
 *
 * @deprecated Untested!
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test!
 **/
class Histogram
{
public:

	//! Add a value to the histogram
	void addValue(int32_t bin, double statisticalWeight=1.0)
	{
		histogram[bin]+=statisticalWeight;
		sum+=statisticalWeight;
	}

	/**
	 * @brief Get the resulting histogram as a table (can be used with ResultFormattingTools::writeResultFile)
	 * the values of the bins run from lowerLimit to upperLimit
	 **/
	std::vector< std::vector < double > > getResultTable(int32_t lowerLimit, int32_t upperLimit)
	{
		size_t vectorLength(upperLimit-lowerLimit+1);
		
		std::vector< std::vector <double> > results;
		results.resize(2);
		results[0].resize(vectorLength);
		results[1].resize(vectorLength);
		
		int32_t mapIndex; 
		size_t vectorIndex;
		
		for(mapIndex=lowerLimit,vectorIndex=0; mapIndex<=upperLimit; ++mapIndex,++vectorIndex)
		{
			results[0][vectorIndex]=double(mapIndex);
			results[1][vectorIndex]=histogram[mapIndex];
		}
		
		return results;
		
	}
	
	//! As getResultTable(int32_t lowerLimit, int32_t upperLimit), but returns the complete range of values
	std::vector< std::vector < double > > getResultTable()
	{
		//get upper and lower limit for the bins
		int32_t lowerLimit=histogram.begin()->first;
		int32_t upperLimit=histogram.rbegin()->first;
		
		return getResultTable(lowerLimit,upperLimit);
		
	}
	
	
	//! As getResultTable(), but the histogram is normalized to 1
	std::vector< std::vector < double > > getNormalizedResultTable()
	{
		//get upper and lower limit for the bins
		int32_t lowerLimit=histogram.begin()->first;
		int32_t upperLimit=histogram.rbegin()->first;
		
		return getNormalizedResultTable(lowerLimit,upperLimit);
		
	}
	
	/**
	 * @brief As getNormalizedResultTable(), but returns only the part of the histogram with
	 * bins between lower limit and upper limit
	 **/
	std::vector< std::vector < double > > getNormalizedResultTable(int32_t lowerLimit,int32_t upperLimit)
	{
	
		std::vector< std::vector < double > > result=getResultTable(lowerLimit,upperLimit);
		if(sum>0.0)
		{
			for (size_t n=0;n<result[0].size();n++)
			{
				result[1][n]/=sum;
			}
		}
		
		return result;
		
	}
	
	
private:
	//! Map of bins to collected values
	std::map<int32_t,double> histogram;

	//! Number of samples used( including statisticalWeight)
	double sum;
};


/////////////////////////////////////////////////////////////////////////////
/**
 * @class SignOfSlope
 *
 * @brief Simple class for holding the slope of preceding value.
 *
 * @deprecated Untested!
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test or delete it!
 **/
struct SignOfSlope
{
  
    double old;
    bool sign;
    
    SignOfSlope():old(0),sign(true){}
    
    void add(double val){sign = val > old ? true : false; }
    bool get() const { return sign;}
    
};

/******************************************************************************
 * THE FOLLOWING CLASSES CAN BE USED ALSO WITH GenericAnalyzer
 * ***************************************************************************/

/*****************************************************************************/
/**
 * @struct Average
 *
 * @brief Collects scalar values and calculates average.
 *
 * @deprecated Untested!
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test or delete it!
 **/
struct Average
{
  
  double sum;
  uint64_t N;
  
  Average():sum(0),N(0){}
  
  void add(double val){sum+=val;++N;}
  
  /**
   * @brief Adds a value to the Average.
   *
   * @details This does essentially the same as the operator+=\n
   * It is still useful though, if one wants to use Average together with
   * GenericAnalyzer<IngredientsType,IngredientsType::molecules_type,Accumulator,PrimitiveAnalyzer>
   * In this case this interface allows for the PrimitiveAnalyzer to return an Average of some
   * changing subgroups and the ensemble average is still calculated correctly.
   **/
  void add(Average val){sum+=val.sum;N+=val.N;}
  
  Average& operator+=(Average& rhs){
	  sum+=rhs.sum;N+=rhs.N;
// 	  std::cout << "added average.\n";
	  return *this;
	  
	}
  
  double get() const { return sum/double(N);}
  
  void reset(){N=0;sum=0.0;}
  
  static std::string getName() { return "Average" ;}
  
  void write ( std::ofstream& out )
  {
	  out << this->get();
  }
  
};

/*****************************************************************************/
/**
 * @class TimeSeries
 *
 * @brief Accumulates a time series of data of type DataType.
 *
 * @details The class requires new values to be added as std::pair<int64_t,DataType>, 
 * where the first value in the pair is an integer time, and the second is the data
 *
 * @deprecated Untested!
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Write a test or delete it!
 **/
template<class DataType> class TimeSeries
{
public:
  TimeSeries(){}
  
  /*****************************************************************************/
  /**
   * @brief Combine two time series into on with averaged data (different lengths are allowed)
   **/
  TimeSeries<DataType>& operator+=(TimeSeries<DataType>& rhs);
  
  /*****************************************************************************/
  /**
   * @brief Add a new data point to the time series. 
   **/
  void add(std::pair<uint64_t,DataType> data);
  void add(uint64_t time, DataType data);
  
  /*****************************************************************************/
  /**
   * @brief Returns the complete, averaged time series as a map of time to value
   **/
  std::map<uint64_t,DataType> get();
  
  /*****************************************************************************/
  /**
   * @brief Returns a name string, e.g. for putting together filenames
   **/
  static std::string getName() { return "TimeSeries" ;}
  
  /*****************************************************************************/
  /**
   * @brief writes formatted result into file stream
   **/
  void write(std::ofstream& out);

  /*****************************************************************************/
  /**
   * @brief clears the data
   **/
  void clear(){data.clear();nValues.clear();}

  /*****************************************************************************/
  /**
   * @brief clears the data
   **/
  size_t size(){return data.size();}
  
  
private:
  /*****************************************************************************/
  /**
   * @var data 
   * @brief Map of timestamp to accumulated data
   * 
   * @var nValues
   * @brief Map of timestamp to number of data points. Final result is data[i]/nValues[i]
   **/
  std::map<uint64_t,DataType>  data;
  std::map<uint64_t,uint64_t>  nValues;
};

/*****************************************************************************/
/******************************************************************************
 * implementation of class TimeSeries
 * ****************************************************************************/

/*****************************************************************************/
//operator+=
/*****************************************************************************/
template<class DataType>
TimeSeries<DataType >& TimeSeries<DataType>::operator+=(TimeSeries< DataType >& rhs)
{
  
    typedef typename std::map<uint64_t,DataType> ::iterator iterator;
    //add data and nValues from rhs. If a timestamp exists in rhs, but not here,
    //it is created due to behaviour of operator[] of std::map. thus, all data survive.
    for(iterator it=rhs.data.begin();it!=rhs.data.end();++it){
      data[it->first]+=it->second;
      nValues[it->first]+=rhs.nValues[it->first];
    }
    //return reference to self
    return *this;
}

/*****************************************************************************/
//add(...). Adding of new values is only allowed for increasing timesteps
/*****************************************************************************/
template<class DataType>
void TimeSeries<DataType>::add(std::pair<uint64_t,DataType> dataPoint)
{
  //add new value to the end of map, if timestep is increasing. throw exception otherwise.
  if(data.size()!=0 && dataPoint.first > data.end()->first){
    data.insert(data.end(),dataPoint);
    nValues.insert(nValues.end(),std::pair<uint64_t,uint64_t>(dataPoint.first,1));
  }
  //if no values there yet, add the new data in any case
  else if(data.size()==0){
    data.insert(data.end(),dataPoint);
    nValues.insert(std::pair<uint64_t,uint64_t>(dataPoint.first,1));
  }
  //throw exception of timestep is not growing
  else
    throw std::runtime_error("TimeSeries::add  Time series is not continuous.");
}

template<class DataType>
void TimeSeries<DataType>::add(uint64_t time, DataType data)
{
	add(std::make_pair(time,data));
}

/*****************************************************************************/
//get()
/*****************************************************************************/
template<class DataType>
std::map<uint64_t, DataType > TimeSeries<DataType>::get()
{
    //calculate average data[i]/nValues[i], store series in map, and return it.
    std::map<uint64_t,DataType> tmp;
  
    typedef typename std::map <uint64_t,DataType>::iterator iterator;
    for (iterator it=data.begin();it!=data.end();++it){
      tmp.insert(std::pair<uint64_t,DataType>(it->first,(it->second)/nValues[it->first]));
    }
    
    return tmp;
}

/*****************************************************************************/
//write(...)
/*****************************************************************************/
template<class DataType>
void TimeSeries<DataType>::write(std::ofstream& out)
{
    //get correctly averaged result
    std::map<uint64_t,DataType> tmp=get();
    typedef typename std::map <uint64_t,DataType>::iterator iterator;
    
    //write complete series to file stream
    for (iterator it=tmp.begin();it!=tmp.end();++it){
      out<<it->first<<"\t"<<it->second<<std::endl;
    }
}

#endif /* LEMONADE_UTILITY_NUMERICTOOLS_H */
