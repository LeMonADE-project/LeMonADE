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
#include <vector>


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
  /**
  * @brief standard constructor for a monomer Label 
  */
  MonomerLabel( const int32_t ChainID_, 
	        const int32_t NLabels_, 
	        const int32_t LabelID_):
  ChainID(ChainID_),NLabels(NLabels_),LabelID(LabelID_)
  { if (ChainID_ < 0 || NLabels_ < 0 || LabelID_ < 0  )
      throw std::runtime_error("ChainID, NLabels and LabelID should not be smaller 0! "); 
  }
  /**
  * @brief default constructor
  */
  MonomerLabel():ChainID(0),NLabels(0),LabelID(0){}
  /**
  * @brief copy constructor
  */
  MonomerLabel(const MonomerLabel &Label )
  {
    
    ChainID = Label.ChainID;
    NLabels = Label.NLabels;
    LabelID = Label.LabelID;    
  }

  /**
  * @brief destructor
  */
  ~MonomerLabel(){};
  
  /**
  * @brief assignment operator
  */
  MonomerLabel& operator= ( const MonomerLabel &Label )
  { 
    ChainID = Label.ChainID;
    NLabels = Label.NLabels;
    LabelID = Label.LabelID;
  }
  
  /**
  * @brief equality operator
  */
  bool operator == (const MonomerLabel& Label )
  { 
    if( Label.ChainID != ChainID ) return false; 
    if( Label.LabelID != LabelID ) return false; 
    if( Label.NLabels != NLabels ) return false; 
    return true;
  };
  
    /**
  * @brief unequality operator
  */
  bool operator != (const MonomerLabel& Label )
  { 
    return !(*this == Label );
  };

  
  //! set the chain ID the Label belongs to 
  void     setChainID(const uint32_t ChainID_){ChainID=ChainID_ ;}
  //! get the chain ID the Label belongs to
  const uint32_t& getChainID(		       ) const { return ChainID  ;}
  
  //! set the number of labels sitting on the monomer 
  void     setNLabels(const uint32_t NLabels_){NLabels=NLabels_ ;}
  //! get the number of labels sitting on the monomer 
  const uint32_t& getNLabels(		       ) const { return NLabels  ;}
  
  //! set the ID of the label(s) on the monomer
  void     setLabelID(const uint32_t LabelID_){LabelID=LabelID_ ;}
  //! get the ID of the label(s) on the monomer
  const uint32_t& getLabelID(		       ) const { return LabelID  ;}
  
  //! set the chain ID, number of labels and label id 
  void     setAll(const uint32_t ChainID_,const uint32_t NLabels_,const uint32_t LabelID_)
  {
    ChainID=ChainID_ ;
    NLabels=NLabels_ ;
    LabelID=LabelID_ ;
  }
  void     setAll(const std::vector<uint32_t> in_)
  {
    ChainID=in_[0] ;
    NLabels=in_[1] ;
    LabelID=in_[2] ;
  }
  //! get the chain ID, number of labels and label id 
  const std::vector<uint32_t>& getAll() const { 
    std::vector<uint32_t> Labels(3,0);
    Labels[0]=(ChainID);
    Labels[1]=(NLabels);
    Labels[2]=(LabelID);
    return Labels;
  }
  
private:
  uint32_t ChainID;
  uint32_t NLabels;
  uint32_t LabelID;
  
};
  
  /**
  * @brief \b Stream \b Out \b operator of the MonomerLabel
  *
  * @details Streams out the elements chain ID, number of labels and the label 
  * 	     ID  of the label separated by space
  *
  * @param stream output-stream
  * @param label object of class MonomerLabel
  * @return output-stream
  **/
  std::ostream& operator<< (std::ostream& stream, const MonomerLabel & label)
  {
	  stream 	  
	  << label.getChainID() << "/" 
	  << label.getNLabels() << "/" 
	  << label.getLabelID();
	  return stream;
  };

  /**
  * @brief \b Stream \b In \b operator of the MonomerLabel
  *
  * @details Streams in the elements. And setting at first chainID, then NLabels, at last LabelID.
  *
  * @param stream input-stream
  * @param label object of class MonomerLabel
  * @return input-stream
  **/
  std::istream& operator >> (std::istream& stream, MonomerLabel & label)
  {
	  int temp;
	  stream >> temp; label.setChainID(temp); stream.ignore(1);
	  stream >> temp; label.setNLabels(temp); stream.ignore(1);
	  stream >> temp; label.setLabelID(temp);
	  return stream;
  }

#endif