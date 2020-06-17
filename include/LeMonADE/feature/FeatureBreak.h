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

#ifndef LEMONADE_FEATURE_FEATUREBREAK_H
#define LEMONADE_FEATURE_FEATUREBREAK_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
/*****************************************************************************/
/**
 * @file
 * @date   2020/06/06
 * @author Toni
 *
 * @class FeatureBreak
 * @brief This Feature add new bonds between monomers.
 *
 * @details Works only in combination with an excluded volume feature
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureBreak : public Feature {
  
public:
	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;

	//constructor
	FeatureBreak():breakingEnergy(0.0),prop(1.0){};
	
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template<class IngredientsType>
	void exportRead(FileImport<IngredientsType>& fileReader);

	//! Export the relevant functionality for writing bfm-files to the responsible writer object
	template<class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const{ return true;};

	//! check bas connect move - always true 
	template<class IngredientsType, class SpecializedMove> 
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const{ return true};
	
        //! check bas connect move - always true 
	template<class IngredientsType, class SpecializedMove> 
	bool checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const;
        
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalSc& move) const {return true;};
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalScDiag& move) const {return true;};
	
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const {throw std::runtime_error("*****FeatureBreak::check MoveLocalBcc: wrong lattice type ... \n"); return false;};
	
	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const {throw std::runtime_error("*****FeatureBreak::check MoveLocalBcc: wrong lattice type ... \n"); return false;};

	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const {return true;};
	
	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move){};
	
	//!apply move for the scBFM connection move for connection
	template<class IngredientsType, class SpecializedMove> 
	void applyMove(IngredientsType& ing, const MoveConnectBase<SpecializedMove>& move);
	
        //!apply move for the scBFM connection move for connection
	template<class IngredientsType, class SpecializedMove> 
	void applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move);
        
	//!apply move for the scBFM local move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);	
	
	//!apply move for the scBFM local diagonal move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalScDiag& move);	

	//!
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);	
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
        
        double getBreakingEnergy(){return breakingEnergy;}
        void setBreakingEnergy(double breakingEnergy_) {breakingEnergy=breakingEnergy_;prop=exp(breakingEnergy);}
protected:
  
private:
        template<class IngredientsType>
        void fillReactivePartner(IngredientsType& ingredients);
        
        std::map <uint32_t, uint32_t > ReactiveBonds;
        
        void removeReactiveBond(uint32_t Mon1, uint32_t Mon2);
        void addReactiveBond(uint32_t Mon1, uint32_t Mon2);
        
        double breakingEnergy;
        double prop;
        
};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/******************************************************************************/
/**
 * @fn void FeatureBreak ::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveConnectBase connects two monomers which later can be broken again. 
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureBreak::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
{
  if ( ! ingredients.getMolecules().areConnected(move.getIndex(),move.getPartner())) return false; 
  move.multiplyProbability();
}

/******************************************************************************/
/**
 * @fn void FeatureBreak ::applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveConnectBase connects two monomers which later can be broken again. 
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
void FeatureBreak::applyMove(IngredientsType& ing, const MoveBreakBase<SpecializedMove>& move)
{
  removeReactiveBond(move.getIndex(),move.getPartner());
}
/******************************************************************************/
/**
 * @fn void FeatureBreak ::applyMove(IngredientsType& ing, const MoveConnectSc& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveConnectBase connects two monomers which later can be broken again. 
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
void FeatureBreak::applyMove(IngredientsType& ing,const MoveConnectBase<SpecializedMove>& move)
{
  auto Mon1(move.getIndex());
  auto Mon2(move.getPartner());
//   if (ing.getMolecules()[Mon1].isReactive() && ing.getMolecules()[Mon2].isReactive() )
    addReactiveBond(Mon1,Mon2);
}
/******************************************************************************/
/**
 * @fn void FeatureBreak ::fillReactivePartner(IngredientsType& ingredients)
 * @brief fillReactivePartner the bond list of the reactive monomers. 
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureBreak::fillReactivePartner(IngredientsType& ingredients)
{
  for (auto i=0; i < ingredients.getMolecules().size();i++)
  {
    if (ingredients.getMolecules[i].isReactive() == true )
    {
      for(auto j=0; j < ingredients.getMolecules().getNumLinks(i);j++)
      {
        auto neighbor(ingredients.getMolecules().getNeighborIdx(i,j));
          if ( neighbor > i && (ingredients.getMolecules[neighbor].isReactive() == true) )
            addReactiveBond(i,neighbor);  
      }
    }
  }
}
/******************************************************************************/
/**
 * @fn void FeatureBreak ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the bond list of the reactive monomers. 
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureBreak  ::synchronize(IngredientsType& ingredients)
{
	std::cout << "FeatureBreak::synchronizing reactive partner list ...\n";
	fillReactivePartner(ingredients);
	std::cout << "done\n";
}

/*****************************************************************/
/**
 * @class ReadBreakingEnergy
 *
 * @brief Handles BFM-File-Reads \b !breaking_energy
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType>
class ReadBreakingEnergy: public ReadToDestination<IngredientsType>
{
public:
  ReadBreakingEnergy(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadBreakingEnergy(){}
  virtual void execute();
};


/*****************************************************************/
/**
 * @class WriteBreakingEnergy
 *
 * @brief Handles BFM-File-Write \b !breaking_energy
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBreakingEnergy:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !breaking_energy into the header of the bfm-file.
  WriteBreakingEnergy(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
  virtual ~WriteBreakingEnergy(){}
  virtual void writeStream(std::ostream& strm);
};


/**
 * @brief Executes the reading routine to extract \b !breaking_energy.
 *
 * @throw <std::runtime_error> attributes and identifier could not be read.
 **/
template < class IngredientsType >
void ReadBreakingEnergy<IngredientsType>::execute()
{
	std::cout<<"reading breaking_energy...";
	double breaking_energy = 0.0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();
	std::string line;
	getline(source,line);
	breaking_energy = atof(line.c_str());
	std::cout << "#!breaking_energy=" << breaking_energy << std::endl;
	ingredients.setBreakingEnergy(breaking_energy);
}


//! Executes the routine to write \b !breaking_energy.
template < class IngredientsType>
void WriteBreakingEnergy<IngredientsType>::writeStream(std::ostream& strm)
{
  strm<<"#!breaking_energy=" << (this->getSource().getBreakingEnergy()) << std::endl<< std::endl;
}
/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * !breaking_energy
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureBreak::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("#!breaking_energy",new ReadBreakingEnergy<IngredientsType>(fileReader.getDestination()));
}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * !breaking_energy
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureBreak::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("#!breaking_energy",new WriteBreakingEnergy<IngredientsType>(fileWriter.getIngredients_()));
}
#endif
