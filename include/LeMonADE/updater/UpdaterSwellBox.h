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

#ifndef LEMONADE_UPDATER_UPDATERSWELLBOX_H
#define LEMONADE_UPDATER_UPDATERSWELLBOX_H

#include <string>
#include <vector>
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 
#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/MonomerGroup.h>
/**
 * @file
 *
 * @class UpdaterSwellBox
 *
 * @brief Increase the box size for a swelling simulation.
 *
 * @details First the largest cluster is detected. If monomers of this cluster are close to the box boundary, 
 * the box size is increased by the step size.
 *
 * @tparam ingredients_ Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * @tparam maxBoxSize_  the maximum box size. Exceeding leads to a runtime_error 
 * @tparam step_ inteval to increase the box size
 */
template<class IngredientsType>
class UpdaterSwellBox:public AbstractUpdater
{
  
public:
  /**
   * @brief Standard Constructor initialized with ref to Ingredients and MCS per cycle
   *
   * @tparam ingredients_ Ingredients class storing all system information( e.g. monomers, bonds, etc).
   * @tparam maxBoxSize_  the maximum box size. Exceeding leads to a runtime_error 
   * @tparam step_ inteval to increase the box size
  */
  UpdaterSwellBox(IngredientsType& ingredients_, const uint32_t maxBoxSize_,  uint32_t step_,  const uint32_t boundarySize_=2 )
  :ingredients(ingredients_), maxBoxSize(maxBoxSize_), step(step_), boundarySize(boundarySize_),nMolecules(0),isExecuted(false)
  {}
 
  /**
   * @brief This checks all used Feature and applies all Feature if all conditions are met.
   *
   * @details 
   *
   * @return True if function are done and should run further steps. Otherwise return false.
   */
  bool execute()
  {
    if(isExecuted==false)
    {
      std::cout << "UpdaterSwellBox::execute() " << std::endl;    
      const typename IngredientsType::molecules_type& getMolies = ingredients.getMolecules(); 
      typename IngredientsType::molecules_type& setMolies = ingredients.modifyMolecules(); 
      //search largest cluster: use the attribute tag for that
      //at first store the initial attributes to reset them after the search
      //then color the whole graph. connected parts have the same color
      auto nMonomers(getMolies.size());
      //store initial attributes
      std::vector<int> tmpAttributes;
      for (auto i =0 ; i < nMonomers;i++)
      {
          tmpAttributes.push_back(getMolies[i].getAttributeTag());
          setMolies[i].setAttributeTag(0);
      }
      //color graph
      int freeColor(1);
      for (auto i=0 ; i < nMonomers;i++)
      {
        std::vector<size_t> rememberBranch;
        rememberBranch.push_back(i);
        if(getMolies[i].getAttributeTag() == 0) 
        {
          while ( rememberBranch.size() > 0 )
          {
            auto ID(rememberBranch.back());
            auto color(getMolies[ID].getAttributeTag());
            rememberBranch.pop_back();
            
            if ( color == 0)
            {
                setMolies[ID].setAttributeTag(freeColor);
                auto nLinks(getMolies.getNumLinks(ID));
                if(nLinks != 0 )
                {
                  for(auto j=0; j < nLinks ; j++)
                  {
                    auto Neighbor(getMolies.getNeighborIdx(ID,j));
                    if( getMolies[Neighbor].getAttributeTag() == 0 ) rememberBranch.push_back(Neighbor);
                  }
                }
            }
          }
          freeColor++;
        }
      } 
      freeColor--;
      nMolecules=freeColor;
      std::vector<std::vector<uint32_t> > ColoredGraphIDs(freeColor+1,std::vector<uint32_t>(0));
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is "<< ColoredGraphIDs.size() << std::endl; 
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is Mol0 "<< ColoredGraphIDs[0].size() << std::endl;
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is Mol1 "<< ColoredGraphIDs[1].size() << std::endl;
      for (auto i=0 ; i < nMonomers;i++)
      {
        uint32_t atti(getMolies[i].getAttributeTag());
//         std::cout << "UpdaterSwellBox::execute(): ID="<<i << " color="  <<atti << std::endl; 
        if (atti==0 )//this should never happen because all monomers should be covered above and colored 
        {
          std::stringstream error_message;
          error_message << "UpdaterSwellBox::execute(): Found monomer which is not colored! ID is " << i << "\n";
          throw std::runtime_error(error_message.str());
        }
        ColoredGraphIDs[atti].push_back(i);
      }
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is "<< ColoredGraphIDs.size() << std::endl; 
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is Mol0 "<< ColoredGraphIDs[0].size() << std::endl;
//       std::cout << "UpdaterSwellBox::execute(): Cluste Size is Mol1 "<< ColoredGraphIDs[1].size() << std::endl;
      auto biggestClusterID(0);
      auto biggestClusterSize(ColoredGraphIDs[biggestClusterID].size());
      for (auto i=0; i < freeColor+1 ; i++)
      {
        std::cout << "ClusterID  "  << i << " cluster size  " << ColoredGraphIDs[i].size() << std::endl; 
        if ( biggestClusterSize < ColoredGraphIDs[i].size() )
        { 
          biggestClusterSize=ColoredGraphIDs[i].size();
          biggestClusterID=i;
        }
      }
      std::cout << "biggest ClusterID  "  << biggestClusterID << " biggestClusterSize " << biggestClusterSize << std::endl; 
      LargestCluster=ColoredGraphIDs[biggestClusterID];
      //reset attributes
      for (auto i =0 ; i < nMonomers;i++)
      {
          tmpAttributes.push_back(getMolies[i].getAttributeTag());
          setMolies[i].setAttributeTag(tmpAttributes[i]);
      }
      //check monomer position 
      auto BoxX(ingredients.getBoxX());
      auto BoxY(ingredients.getBoxY());
      auto BoxZ(ingredients.getBoxZ());
      auto thresholdX(BoxX-boundarySize);
      auto thresholdY(BoxY-boundarySize);
      auto thresholdZ(BoxZ-boundarySize);
      bool increaseBoxSize(false);
      for (auto i =0 ; i < LargestCluster.size();i++)
      {
        auto pos(getMolies[LargestCluster[i]].getVector3D());
        if(pos.getX()>BoxX-boundarySize || pos.getX() < boundarySize ){increaseBoxSize=true;break; }
        if(pos.getY()>BoxY-boundarySize || pos.getY() < boundarySize ){increaseBoxSize=true;break; }
        if(pos.getZ()>BoxZ-boundarySize || pos.getZ() < boundarySize ){increaseBoxSize=true;break; }
      }
      if (increaseBoxSize)
      {
//         std::cout << "UpdaterSwellBox::execute() Set the box to the size:  "<<  BoxX+step << std::endl;
        ingredients.setBoxX(BoxX+step);
        ingredients.setBoxY(BoxY+step);
        ingredients.setBoxZ(BoxZ+step);
        ingredients.synchronize();
        auto newBoxX(ingredients.getBoxX());
        auto newBoxY(ingredients.getBoxY());
        auto newBoxZ(ingredients.getBoxZ());
        for (auto i =0; i < getMolies.size();i++)
        {
          auto oldPos(getMolies[i].getVector3D());
          setMolies[i].modifyVector3D().setAllCoordinates( translateCoordinate(oldPos.getX(), BoxX, newBoxX)
                                                          ,translateCoordinate(oldPos.getY(), BoxY, newBoxY)
                                                          ,translateCoordinate(oldPos.getZ(), BoxZ, newBoxZ)
                                                        );
        }
        std::cout << "UpdaterSwellBox::initialize-> execute:  increase box size to \n"
                  << " x-Box: " << ingredients.getBoxX() <<"\n"
                  << " y-Box: " << ingredients.getBoxY() <<"\n"
                  << " z-Box: " << ingredients.getBoxZ() <<"\n"
                  << "\n";
        ingredients.synchronize();
      }
      isExecuted=true;
      return true;  
    }else 
      return false;
  }

  /**
   * @brief This function is called \a once in the beginning of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for initializing tasks (e.g. init SDL)
   *
   **/
  virtual void initialize(){
    std::cout << "UpdaterSwellBox::initialize() " << std::endl;
    
    if ( ingredients.isPeriodicX() || ingredients.isPeriodicY() || ingredients.isPeriodicZ())
    {
      std::stringstream error_message;
      error_message << "UpdaterSwellBox::initialize: System is periodic in at least on direction ("
                    << ingredients.isPeriodicX() <<","<< ingredients.isPeriodicY() <<","<< ingredients.isPeriodicZ()  <<")\n";
      throw std::runtime_error(error_message.str()); 
    }
    execute();
  }

  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup(){};

  //!
  inline size_t getNMolecules(){return nMolecules;};
  inline size_t getLargestClusterSize(){return LargestCluster.size();};
private:
  //! A reference to the IngredientsType - mainly the system
  IngredientsType& ingredients;
  //!maximum box size
  const uint32_t maxBoxSize;
  //!step width to increase the box 
   uint32_t step;  
  //!distance of monomer to boundary which is checked
  const uint32_t boundarySize; 
  //! number of molecules_type
  uint32_t nMolecules;
  //!execute only in the beginning of the simulation to avoid problems in the bfm file writing
  bool isExecuted;
  //! group of monomers which belong together
//   std::vector<MonomerGroup<molecules_type> > groups;
  //! largest cluster which is associated with the network
  std::vector<uint32_t> LargestCluster;
  //! translation of the coordinate relative to the mid point 
  template <class T >
  inline T translateCoordinate(T pos, uint32_t oldBox, uint32_t newBox)
  {
    return (newBox+oldBox)/2-pos;
  }
};

#endif
