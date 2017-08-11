/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Hauke Rabbel)
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

/* ****************************************************************************
 * This example shows how to write a simple analyzer. This analyzer calculates
 * the center of mass of all the particles in the system and writes an output
 * file at the end. You can use this example as a starting point for writing
 * any analyzers you need later.
 * ***************************************************************************/

#include <string>
#include <ostream>            //for writing the output file

#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/Vector3D.h>

/* ****************************************************************************
 * Here comes the declaration of the class. Analyzers generally inherit
 * from the class AbstractAnalyzer and have to implement the three functions
 * initialize(), execute() and cleanup().
 * Furthermore, they normally hold a const reference to either ingredients
 * or molecules, because they need to access the particle data. The reference
 * should always be const, because the Analyzer is not supposed to change the
 * system, only analyze the contents.
 * Below the class declaration we implement the functions
 * ***************************************************************************/


template<class IngredientsType> class Ex5Analyzer:public AbstractAnalyzer
{
public:

  //constructor, get a reference to the ingredients for accessing the data later
  Ex5Analyzer(const IngredientsType& ing, std::string filename);

  //we have to implement these three methods
  void initialize();
  bool execute();
  void cleanup();

private:

  const IngredientsType& ingredients;

  std::string outputFilename;
  //this will count the number of samples, so we can calculate the average
  uint32_t nValues;
  //this holds the number of particles in the system. set by initialize()
  double nParticles;
  //this will be the result at the end
  VectorDouble3 centerOfMass;

};


/* ****************************************************************************
 * Here come the implementations
 * ***************************************************************************/

//constructor
//we initialize the reference to ingredients and set nValues to zero
template<class IngredientsType>
 Ex5Analyzer<IngredientsType>::Ex5Analyzer(const IngredientsType& ing,
			 std::string filename)
 :ingredients(ing)
{
  outputFilename=filename;
  nValues=0;
}

/* ****************************************************************************
 * initialize
 * in initialize we save the number of particles to a variable
 * ***************************************************************************/
template<class IngredientsType>
void Ex5Analyzer<IngredientsType>::initialize()
{
  nParticles=double(ingredients.getMolecules().size());
}



/* ****************************************************************************
 * execute...this is where the analysis happens
 * we calculate the center of mass of the system and add the value to the
 * vector centerOfMass. In cleanup, we will divide the collected value
 * by the number of samples.
 * ***************************************************************************/
template<class IngredientsType>
bool Ex5Analyzer<IngredientsType>::execute()
{
  VectorDouble3 currentCenterOfMass;
  //loop over all particles and calculate the current center of mass
  for(size_t n=0;n<nParticles;n++){
    VectorDouble3 particleCoordinate=ingredients.getMolecules()[n];
    currentCenterOfMass+= (particleCoordinate/nParticles);
  }

  //add the current center of mass to the accumulated value
  centerOfMass+=currentCenterOfMass;
  //increase the counter for the number of samples
  nValues++;

  return true;
}

/* ****************************************************************************
 * cleanup()... here we calculate the average and write an output file
 * ***************************************************************************/
template<class IngredientsType>
void Ex5Analyzer<IngredientsType>::cleanup()
{
  //calculate the average
  centerOfMass/=double(nValues);

  //open an output file and write the center of mass to it
  //we use the filename given to the constructor
  std::ofstream outputFile;
  outputFile.open(outputFilename.c_str());
  outputFile<<centerOfMass.getX()<<"\t"
            <<centerOfMass.getY()<<"\t"
            <<centerOfMass.getZ()<<std::endl;

  outputFile.close();

}
