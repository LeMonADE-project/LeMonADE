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

#ifndef LEMONADE_ATTRIBUTE_TAG_H
#define LEMONADE_ATTRIBUTE_TAG_H

#include <stdexcept>
#include <iostream>
#include <sstream>


/*****************************************************************************/
/**
 * @file AttributeTag.h
 *
 * @brief Provides simple classes for different tags. 
 *
 * @details 
 *
 * @todo 
 **/


typedef  int32_t Integer;


/*****************************************************************************/
/**
 * @file AttributeTag.h
 *
 * @brief Class for labels on chains
 *
 * @details The class stores the chainID it belongs to and the number of Labels
 * 	    sitting on the monomer. 
 *
 * @todo 
 **/
class MonomerLabel
{
public:
  //Default contructor 
  MonomerLabel():ChainID(0),NLabels(0),LabelID(0){};
  //@todo:copy constructor
  //
  MonomerLabel(uint32_t ChainID_, uint32_t NLabels_, uint32_t LabelID_):ChainID(ChainID_),NLabels(NLabels_),LabelID(LabelID_){};
  
  ~MonomerLabel(){};
  //assign operator 
  MonomerLabel& operator = ( const MonomerLabel &Label )
  {  };
  //equality operator 
  bool operator==(MonomerLabel& Label )
  {  };
  //! set the chain ID the Label belongs to 
  void     setChainID(uint32_t ChainID_){ChainID=ChainID_ ;}
  //! get the chain ID the Label belongs to
  uint32_t getChainID(		       ){ return ChainID  ;}
  
  //! set the number of labels sitting on the monomer 
  void     setNLabels(uint32_t NLabels_){NLabels=NLabels_ ;}
  //! get the number of labels sitting on the monomer 
  uint32_t getNLabels(		       ){ return NLabels  ;}
  
  //! set the ID of the label(s) on the monomer
  void     setLabelID(uint32_t LabelID_){LabelID=LabelID_ ;}
  //! get the ID of the label(s) on the monomer
  uint32_t getLabelID(		       ){ return LabelID  ;}
  
  
private:
  uint32_t ChainID;
  uint32_t NLabels;
  uint32_t LabelID;
  
};
#endif