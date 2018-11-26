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

#ifndef LEMONADE_FEATURE_FEATUREMOLECULESIOUNSAVECHECK_H
#define LEMONADE_FEATURE_FEATUREMOLECULESIOUNSAVECHECK_H

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/core/MoleculesWrite.h>
#include <LeMonADE/core/MoleculesRead.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondsetUnsaveCheck.h>
#include <LeMonADE/io/FileImport.h>

/**
 * @file class FeatureMoleculesIOUnsaveCheck
 * @author Toni
 * */

/**
 * @class FeatureMoleculesIOUnsaveCheck
 * @brief Feature providing standard read/write commands for molecules
 * @details Feature providing standard read/write commands for molecules. The commands
 * provided are !number_of_monomers,!bonds,!add_bonds,!mcs. Uses FeatureBondsetUnsaveCheck
 * which allows bonds between refolded monomers. 
 * */
class FeatureMoleculesIOUnsaveCheck:public Feature
{
public:
	typedef LOKI_TYPELIST_2(FeatureBox,FeatureBondsetUnsaveCheck< >) required_features_back;
	
	FeatureMoleculesIOUnsaveCheck():numberOfMonomers(0),maxConnectivity(0){}
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template<class IngredientsType> 
	void exportRead(FileImport<IngredientsType>& fileReader);
	
	//! Export the relevant functionality for writing bfm-files to the responsible writer object
	template<class IngredientsType> 
	void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

	template<class IngredientsType> 
	void synchronize(IngredientsType& ingredients);
	
	/**
	* @brief Overloaded function to stream all metadata to an output stream
	* @details Gives the number of vertices (monomers), and maximum connectivity
	* @param stream output stream
	*/
	void printMetaData(std::ostream& streamOut) const
	{
		streamOut << "\tNumber of monomers: " << numberOfMonomers << std::endl;
		streamOut << "\tmax connectivity: " << maxConnectivity << std::endl;
	}
	
	//! add a range of particle indices to be written in compressed (solvent) format
	void setCompressedOutputIndices(size_t startIdx, size_t stopIdx)
	{
		solventIndices[startIdx]=stopIdx;
	}
	
	//! add a range of particle indices to be written in compressed (solvent) format
		void clearCompressedOutputIndices()
		{
			solventIndices.clear();//[startIdx]=stopIdx;
		}

	/**
	 * @brief get a map containing ranges of indices which are written in compressed (solvent) format
	 * @details The key and value provide ranges of particle indices to be compressed. values are inclusive,
	 * i.e. a pair of key=5, value=7 means that particles 5,6,7 are compressed
	 * */
	const std::map<size_t,size_t>& getCompressedOutputIndices() const {return solventIndices;}

private:
	//metadata
	uint numberOfMonomers,maxConnectivity;
	//monomers to be compressed
	std::map<size_t,size_t> solventIndices;
	
};






/**
* The function is called by the Ingredients class when an object of type Ingredients 
* is associated with an object of type FileImport. The export of the Reads is thus
* taken care automatically when it becomes necessary.\n
* Registered Read-In Commands:
* * !number_of_monomers
* * !bonds
* * !add_bonds
* * !remove_bonds
* * !mcs
*
* @param fileReader File importer for the bfm-file
* @tparam IngredientsType Features used in the system. See Ingredients.
*/

template < class IngredientsType >
void FeatureMoleculesIOUnsaveCheck::exportRead(FileImport<IngredientsType>& fileReader)
{
	IngredientsType& destination=fileReader.getDestination();
	fileReader.registerRead("!number_of_monomers", new  ReadNrOfMonomers< IngredientsType > (destination) );
	fileReader.registerRead("!bonds", new  ReadBonds< IngredientsType > (destination) );
	fileReader.registerRead("!add_bonds", new  ReadBonds< IngredientsType > (destination) );
	fileReader.registerRead("!remove_bonds", new  ReadRemoveBonds< IngredientsType > (destination) );
	fileReader.registerRead("!mcs", new  ReadMcs< IngredientsType > (destination) );
}


/**
* The function is called by the Ingredients class when an object of type Ingredients 
* is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
* taken care automatically when it becomes necessary.\n
* Registered Write-Out Commands:
* * !number_of_monomers
* * !bonds
* * !add_bonds
* * !remove_bonds
* * !mcs
*
* @param fileWriter File writer for the bfm-file.
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void FeatureMoleculesIOUnsaveCheck::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
	const IngredientsType& source=fileWriter.getIngredients_();
	fileWriter.registerWrite("!number_of_monomers", new WriteNrOfMonomers <IngredientsType> (source));
	fileWriter.registerWrite("!bonds", new WriteBonds <IngredientsType> (source));
	fileWriter.registerWrite("!add_bonds", new WriteAddBonds<IngredientsType> (source,fileWriter.getCommandType()));
	fileWriter.registerWrite("!remove_bonds", new WriteRemoveBonds<IngredientsType> (source,fileWriter.getCommandType()));
	fileWriter.registerWrite("!mcs", new WriteMcs <IngredientsType> (source));
}





template<class IngredientsType> 
void FeatureMoleculesIOUnsaveCheck::synchronize(IngredientsType& ingredients)
{
	numberOfMonomers=ingredients.getMolecules().size();
	maxConnectivity=ingredients.getMolecules().getMaxConnectivity();
}

#endif /* LEMONADE_FEATURE_FEATUREMOLECULESIO_H */
