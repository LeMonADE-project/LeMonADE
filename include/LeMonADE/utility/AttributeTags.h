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

/**
 * @brief this is just a short form 
 */
typedef  int32_t Integer;


/*****************************************************************************/
/**
 * @file AttributeTag.h
 *
 * @brief Class for labels on chains
 *
 * @details The class is a intended to be used as a monomer extension which stores information 
 *	    about the chain ID, the number of labels sitting on the monomer and the label ID.
 *
 * @param chainID the monomer is intended to be in a chain of a certain ID
 * @param nLabels stores the number of labels sitting on the monomer
 * @param labelID stores the label ID (like a color). It is implicitly assumed 
 * 		  that all labels sitting on the same monomer have the same label ID!
 * @todo 
 **/
class MonomerLabel
{
public:
  /**
  * @brief standard constructor for a monomer Label 
  */
  MonomerLabel( const int32_t chainID_, 
	        const int32_t nLabels_, 
	        const int32_t labelID_):
  chainID(chainID_),nLabels(nLabels_),labelID(labelID_)
  { if (chainID_ < 0 || nLabels_ < 0 || labelID_ < 0  )
      throw std::runtime_error("chainID, nLabels and labelID should not be smaller 0! "); 
  }
  /**
  * @brief default constructor
  */
  MonomerLabel():chainID(0),nLabels(0),labelID(0){}
  /**
  * @brief copy constructor
  */
  MonomerLabel(const MonomerLabel &Label )
  {
    
    chainID = Label.chainID;
    nLabels = Label.nLabels;
    labelID = Label.labelID;    
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
    chainID = Label.chainID;
    nLabels = Label.nLabels;
    labelID = Label.labelID;
  }
  
  /**
  * @brief equality operator
  */
  bool operator == (const MonomerLabel& Label )
  { 
    if( Label.chainID != chainID ) return false; 
    if( Label.labelID != labelID ) return false; 
    if( Label.nLabels != nLabels ) return false; 
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
  void     setchainID(const uint32_t chainID_){chainID=chainID_ ;}
  //! get the chain ID the Label belongs to
  const uint32_t& getchainID(		       ) const { return chainID  ;}
  
  //! set the number of labels sitting on the monomer 
  void     setnLabels(const uint32_t nLabels_){nLabels=nLabels_ ;}
  //! get the number of labels sitting on the monomer 
  const uint32_t& getnLabels(		       ) const { return nLabels  ;}
  
  //! set the ID of the label(s) on the monomer
  void     setlabelID(const uint32_t labelID_){labelID=labelID_ ;}
  //! get the ID of the label(s) on the monomer
  const uint32_t& getlabelID(		       ) const { return labelID  ;}
  
  //! set the chain ID, number of labels and label id 
  void     setAll(const uint32_t chainID_,const uint32_t nLabels_,const uint32_t labelID_)
  {
    chainID=chainID_ ;
    nLabels=nLabels_ ;
    labelID=labelID_ ;
  }
  //!set all information 
  void     setAll(const std::vector<uint32_t> in_)
  {
    chainID=in_[0] ;
    nLabels=in_[1] ;
    labelID=in_[2] ;
  }
  //! get the chain ID, number of labels and label id 
  const std::vector<uint32_t> getAll() const { 
    std::vector<uint32_t> label(3,0);
    label[0]=(chainID);
    label[1]=(nLabels);
    label[2]=(labelID);
    return label;
  }
  //! reset all info to the default value
  void reset()
  {
    chainID=0;
    nLabels=0;
    labelID=0;
  }
  /**
  * @brief \b Stream \b Out \b operator of the MonomerLabel
  *
  * @details Streams out the elements chain ID, number of labels and the label 
  * 	     ID  of the label separated by slash 
  *
  * @param stream output-stream
  * @param label object of class MonomerLabel
  * @return output-stream
  **/
  friend std::ostream& operator<< (std::ostream& stream, const MonomerLabel & label)
  {
	  stream 	  
	  << label.getchainID() << "/" 
	  << label.getnLabels() << "/" 
	  << label.getlabelID();
	  return stream;    
  }
  /**
  * @brief \b Stream \b In \b operator of the MonomerLabel
  *
  * @details Streams in the elements. And setting at first chainID, then nLabels, at last labelID,
  * 	     which are assumed to be separated by slashes. Otherwise throw a runtime error.
  *
  * @param stream input-stream
  * @param label object of class MonomerLabel
  * @return input-stream
  **/
  friend std::istream& operator >> (std::istream& stream, MonomerLabel & label)
  {
    	  int temp;
	  stream >> temp; label.setchainID(temp);
	  if (stream.peek()=='/') 
	      stream.ignore(1);
	  else 
	    throw std::runtime_error("MonomerLabel::operator>> Wrong separator!");
	  stream >> temp; label.setnLabels(temp); 
	  if (stream.peek()=='/') 
	      stream.ignore(1);
	  else 
	    throw std::runtime_error("MonomerLabel::operator>> Wrong separator!");
	  stream >> temp; label.setlabelID(temp);
	  return stream;
  }
private:
  //! is intended to store the chain ID of the monomer where the label is sitting on 
  uint32_t chainID;
  //! the number of labels sitting on the monomer 
  uint32_t nLabels;
  //! the label ID, because there could be more than only one (like a color)
  uint32_t labelID;
  
};


#endif
