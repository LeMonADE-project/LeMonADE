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

#ifndef FEATURE_SHEARFORCE_H
#define FEATURE_SHEARFORCE_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves//MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>


/**
 * @file class FeatureShearForce
 * */

/**
 * @class FeatureShearForce
 * @brief Feature providing a shear force 
 * @details Feature providing a shear force on the top and bottom of z in x direction 
 * and in the opposite direction in the middle of the box
 * */
class FeatureShearForce:public Feature
{
public:

	FeatureShearForce(): ForceOn(true),Amplitude_Shear_Force(0.0){};
	virtual ~FeatureShearForce(){};
	
	//We need 3 Features to handle the oscillatory forces
	typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;
	
	//! For all unknown moves: this does nothing
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const;
	
	//! Overloaded for moves of type MoveScMonomer to check for sinusoidal movement
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;
	
	//! Overloaded for moves of type MoveScMonomer to check for sinusoidal movement
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalScDiag& move) const;
	
	void setAmplitudeShearForce(double amplitudeShearForce)
	{
		Amplitude_Shear_Force = amplitudeShearForce;
	}

	void setForceOn(bool forceOn)
	{
		ForceOn = forceOn;
	}

	double getAmplitudeShearForce() const
	{
		return Amplitude_Shear_Force;
	}

	bool isForceOn() const
	{
		return ForceOn;
	}

	double getTotalShearForce() const
	{
		return Total_Shear_Force;
	}

	void setTotalShearForce(double totalShearForce)
	{
		Total_Shear_Force = totalShearForce;
	}

	//template<class IngredientsType>
	//void synchronize(IngredientsType& ingredients);
	
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template <class IngredientsType>
	void exportRead(FileImport <IngredientsType>& fileReader);
	
	 //! Export the relevant functionality for writing bfm-files to the responsible writer object
	template <class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;



private:
	
	bool ForceOn; //!< Force is On (True) or Off (False)
	double Total_Shear_Force; //!< The total deformation force in MCS on charge
	double Amplitude_Shear_Force; //!< The amplitude of the sinusoidal force on charge

};

template<class IngredientsType>
bool FeatureShearForce::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}

template<class IngredientsType>
bool FeatureShearForce::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
	if(ForceOn){       
        if(move.getDir()==VectorInt3(-1,0,0)||move.getDir()==VectorInt3(1,0,0)){
		      //position perpendicular to force direction
		      uint32_t posPerp(ingredients.getMolecules()[move.getIndex()].getZ());
		      double total_Shear_Force = Amplitude_Shear_Force*((double)(move.getDir().getX()));
		      //Metropolis: zeta = exp (-dV)
		      //dV=-q*E*dr
		      
		      //shearing works like:
		      // >>>>>>>>>>>>>>>>>>> force to the rigth
		      //
		      //  nothing happens here
		      //
		      // <<<<<<<<<<<<<<<<<<< force to the left
		      //
		      //  nothing happens here
		      //
		      // >>>>>>>>>>>>>>>>>>> force to the rigth
		      //In this way the periodic boundaries stabilize the system and it does not rotate.
		      
		      //Box Size in Z direction
		      uint32_t BoxSize(ingredients.getBoxZ());
		      //foldback into Box
		      posPerp=posPerp%BoxSize;
		      //upper force field acts to the right
		      if(posPerp>BoxSize-3){
				      double prob=exp(total_Shear_Force);
// 				      std::cout<<"PosPerUp "<<posPerp<<"prob "<<prob<<" "<<move.getDir().getX()<<std::endl;
				      move.multiplyProbability(prob);
		      }
		      //center force field acts to the left
		      if(posPerp<=((double)BoxSize/2.+1) && posPerp>((double)BoxSize/2.-3.)){
				      double prob=exp(-total_Shear_Force);
// 				      std::cout<<"PosPerMid "<<posPerp<<"prob "<<prob<<" "<<move.getDir().getX()<<std::endl;
				      move.multiplyProbability(prob);
		      }
		      
		      //lower force field acts to the rigth
		      if(posPerp<2){
				      double prob=exp(total_Shear_Force);
// 				      std::cout<<"PosPerDown "<<posPerp<<"prob "<<prob<<" "<<move.getDir().getX()<<std::endl;
				      move.multiplyProbability(prob);
		      }
		}

	}

	return true;
}

template<class IngredientsType>
bool FeatureShearForce::checkMove(const IngredientsType& ingredients, MoveLocalScDiag& move) const
{
	if(ForceOn){       
        if(move.getDir().getX()==-1||move.getDir().getX()==1){
		      //position perpendicular to force direction
		      uint32_t posPerp(ingredients.getMolecules()[move.getIndex()].getZ());
		      double total_Shear_Force = Amplitude_Shear_Force*((double)(move.getDir().getX()));
		      //Metropolis: zeta = exp (-dV)
		      //dV=-q*E*dr
		      
		      //shearing works like:
		      // >>>>>>>>>>>>>>>>>>> force to the rigth
		      //
		      //  nothing happens here
		      //
		      // <<<<<<<<<<<<<<<<<<< force to the left
		      //
		      //  nothing happens here
		      //
		      // >>>>>>>>>>>>>>>>>>> force to the rigth
		      //In this way the periodic boundaries stabilize the system and it does not rotate.
		      
		      //Box Size in Z direction
		      uint32_t BoxSize(ingredients.getBoxZ());
		      //foldback into Box
		      posPerp=posPerp%BoxSize;
		      //upper force field acts to the right
		      if(posPerp>BoxSize-3){
				      double prob=exp(total_Shear_Force);
				      move.multiplyProbability(prob);
		      }
		      //center force field acts to the left
		      if(posPerp<=((double)BoxSize/2.+1) && posPerp>((double)BoxSize/2.-3.)){
				      double prob=exp(-total_Shear_Force);
				      move.multiplyProbability(prob);
		      }
		      
		      //lower force field acts to the rigth
		      if(posPerp<2){
				      double prob=exp(total_Shear_Force);
				      move.multiplyProbability(prob);
		      }
		}

	}

	return true;
}

/*****************************************************************/
/**
 * @class ReadForceFieldOn
 *
 * @brief Handles BFM-File-Reads \b #!force_field_on
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadForceFieldOn: public ReadToDestination<IngredientsType>
{
public:
  ReadForceFieldOn(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadForceFieldOn(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadForceFieldOn<IngredientsType>::execute()
{
	std::cout<<"reading FieldIsOn...";

	bool fieldIsOn = false;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	fieldIsOn = atoi(line.c_str());
	std::cout << "#!force_field_on = " << (fieldIsOn? "True" : " False" ) << std::endl;

	ingredients.setForceOn(fieldIsOn);
}

/*****************************************************************/
/**
 * @class WriteForceFieldOn
 *
 * @brief Handles BFM-File-Write \b #!force_field_on
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteForceFieldOn:public AbstractWrite<IngredientsType>
{
public:
	WriteForceFieldOn(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteForceFieldOn(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteForceFieldOn<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!force_field_on=" << (this->getSource().isForceOn() ? "1" : "0") << std::endl<< std::endl;
}


/*****************************************************************/
/**
 * @class ReadAmplitudeShearForce
 *
 * @brief Handles BFM-File-Reads \b #!shear_force_amplitude
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadAmplitudeShearForce: public ReadToDestination<IngredientsType>
{
public:
  ReadAmplitudeShearForce(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadAmplitudeShearForce(){}
  virtual void execute();
};

template<class IngredientsType>
void ReadAmplitudeShearForce<IngredientsType>::execute()
{
	std::cout<<"reading AmplitudeShearForce...";

	double amplitudeForce = 0.0;
	IngredientsType& ingredients=this->getDestination();
	std::istream& source=this->getInputStream();

	std::string line;
	getline(source,line);
	amplitudeForce = atof(line.c_str());
	std::cout << "#!shear_force_amplitude = " << (amplitudeForce) << std::endl;

	ingredients.setAmplitudeShearForce(amplitudeForce);
}


/*****************************************************************/
/**
 * @class WriteAmplitudeShearForce
 *
 * @brief Handles BFM-File-Write \b #!shear_force_amplitude
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteAmplitudeShearForce:public AbstractWrite<IngredientsType>
{
public:
	WriteAmplitudeShearForce(const IngredientsType& i)
	:AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

	virtual ~WriteAmplitudeShearForce(){}

	virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteAmplitudeShearForce<IngredientsType>::writeStream(std::ostream& stream)
{
	stream<<"#!shear_force_amplitude=" << (this->getSource().getAmplitudeShearForce()) << std::endl<< std::endl;
}



/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!shear_field_on
 * * #!shear_force
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureShearForce::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!force_field_on", new ReadForceFieldOn<FeatureShearForce>(*this));
    fileReader.registerRead("#!shear_force_amplitude", new ReadAmplitudeShearForce<FeatureShearForce>(*this));

}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * #!shear_field_on
 * * #!shear_force
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureShearForce::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
    fileWriter.registerWrite("#!force_field_on",new WriteForceFieldOn<FeatureShearForce>(*this));
    fileWriter.registerWrite("#!shear_force_amplitude", new WriteAmplitudeShearForce<FeatureShearForce>(*this));

}

#endif /*FEATURE_SHEARFORCE_H*/
