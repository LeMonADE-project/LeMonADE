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

#ifndef LEMONADE_UPDATER_UPDATERSIMPLESIMULATOR_H
#define LEMONADE_UPDATER_UPDATERSIMPLESIMULATOR_H

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/updater/moves/MoveLocalBase.h>

//simple simulation updater. takes the type of move as template 
//argument and the number of mcs to be executed as argument for the 
//constructor

/**
 * @file
 *
 * @class UpdaterSimpleSimulator
 *
 * @brief Simple simulation updater for general purpose.
 *
 * @details It takes the type of move as template argument MoveType
 * and the number of mcs to be executed as argument for the constructor
 *
 * @todo remove IngredientsUpdater and use AbstractUpdater
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam <MoveType> name of the specialized move.
 */
template<class IngredientsType,class MoveType>
class UpdaterSimpleSimulator:public IngredientsUpdater<IngredientsType>
{

public:
  /**
   * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
   *
   * @param ing a reference to the IngredientsType - mainly the system
   * @param steps MCS per cycle to performed by execute()
   */
  UpdaterSimpleSimulator(IngredientsType& ing,uint32_t steps)
  :IngredientsUpdater<IngredientsType>(ing),nsteps(steps)
  {}
  
  /**
   * @brief This checks all used Feature and applies all Feature if all conditions are met.
   *
   * @details This function runs over \a steps MCS and performs the moves.
   * It setting the age of the system and prints a simple simple simulation speed
   * in the number of attempted monomer moves per s (tried and performed monomer moves/s).
   *
   * @return True if function are done.
   */
  bool execute()
  {
	time_t startTimer = time(NULL); //in seconds
	std::cout<<"mcs "<<this->getIngredients().getMolecules().getAge() << " passed time " << ((difftime(time(NULL), startTimer)) ) <<std::endl;

    for(int n=0;n<nsteps;n++){
      
	for(int m=0;m<this->getIngredients().getMolecules().size();m++)
	{
		move.init(this->getIngredients());
		
		if(move.check(this->getIngredients())==true)
		{
			move.apply(this->getIngredients());
		}
	}
      
    }
    
    this->getIngredients().modifyMolecules().setAge(this->getIngredients().modifyMolecules().getAge()+nsteps);
    
    std::cout<<"mcs "<<this->getIngredients().getMolecules().getAge() << " with " << (((1.0*nsteps)*this->getIngredients().getMolecules().size())/(difftime(time(NULL), startTimer)) ) << " [attempted moves/s]" <<std::endl;
    std::cout<<"mcs "<<this->getIngredients().getMolecules().getAge() << " passed time " << ((difftime(time(NULL), startTimer)) ) << " with " << nsteps << " MCS "<<std::endl;

    return true;
  }
  
  /**
   * @brief This function is called \a once in the beginning of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for initializing tasks (e.g. init SDL)
   *
   **/
  virtual void initialize(){};

  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup(){};

private:
  //! Specialized move to be used
  MoveType move;

  //! Number of mcs to be executed
  uint32_t nsteps;
};

#endif
