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
 * This example shows how to write a simple updater. In this case, the updater
 * will add a number of single monomers to the system.
 * Basically, this is similar to what has been done in the previous
 * example, only that it is encapsulated in an updater - and thus it is reusable
 * in other programs.
 * ***************************************************************************/

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

/* ****************************************************************************
 * Here comes the declaration of the class. Updaters generally inherit
 * from the class AbstractUpdater and have to implement the three functions
 * initialize(), execute() and cleanup().
 * Furthermore, they normally hold a reference to either ingredients
 * or molecules, because they need to access the particle data.
 * Below the class declaration we implement the functions
 * ***************************************************************************/

template<class IngredientsType> class Ex5Updater:public AbstractUpdater
{
public:

  //constructor, get a reference to the ingredients for accessing the data later
  //the second argument specifies how many particles to construct
  Ex5Updater(IngredientsType& ing, uint32_t nParticles);

  //we have to implement these three methods
  void initialize();
  bool execute();
  void cleanup();

private:

  IngredientsType& ingredients;
  uint32_t numberOfParticles;

  uint32_t boxX,boxY,boxZ;
};


/* ****************************************************************************
 * Here come the implementations
 * ***************************************************************************/

//constructor
//we initialize the reference to ingredients and set numberOfParticles
template<class IngredientsType>
 Ex5Updater<IngredientsType>::Ex5Updater(IngredientsType& ing,
			uint32_t nParticles)
 :ingredients(ing)
{
  numberOfParticles=nParticles;
}

/* ****************************************************************************
 * initialize...in this case we save the box size to three separate coordinates.
 * note that for this to work we have to use the FeatureBox in the main program
 * ***************************************************************************/
template<class IngredientsType>
void Ex5Updater<IngredientsType>::initialize()
{
  boxX=ingredients.getBoxX();
  boxY=ingredients.getBoxY();
  boxZ=ingredients.getBoxZ();
}

/* ****************************************************************************
 * execute...this is where the particle creation happens
 * ***************************************************************************/
template<class IngredientsType>
bool Ex5Updater<IngredientsType>::execute()
{

  //first we need random numbers, because we place the particles at random
  //positions. therefore, we get an instance of the RandomNumberGenerators
  //class. there is no need to seed the generator here, as it should be done
  //in the main() function of the program
  RandomNumberGenerators randomNumbers;

  //now place the particles in the box. for this example we neglect excluded
  //volume. for each particle we get three random coordinates
  for(uint32_t n=0;n<numberOfParticles;n++){
    uint32_t x=randomNumbers.r250_rand32() % boxX;
    uint32_t y=randomNumbers.r250_rand32() % boxY;
    uint32_t z=randomNumbers.r250_rand32() % boxZ;

    //now we add the particle
    //Note: in the example about monte carlo moves, there will be a better way
    //introduced to add particles, but this also works and is easier for now.
    ingredients.modifyMolecules().addMonomer(x,y,z);

  }

  return true;
}

/* ****************************************************************************
 * cleanup()... no need to do anything here
 * ***************************************************************************/
template<class IngredientsType>
void Ex5Updater<IngredientsType>::cleanup()
{
}
