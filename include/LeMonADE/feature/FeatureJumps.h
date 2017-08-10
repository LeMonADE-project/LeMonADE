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

#ifndef LEMONADE_FEATURE_FEATUREJUMPS_H
#define LEMONADE_FEATURE_FEATUREJUMPS_H

#include <vector>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureBox.h>

#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/FileImport.h>

/**
 * @class ReadJumpsMcs
 *
 * @brief Handles BFM-File-Reads \b !mcs when jumps are present (deprecated bfm-version).
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType >
class ReadJumpsMcs: public ReadToDestination < IngredientsType >
{
  //necessary because only on the first call the connections have to be set up
  //and the get pointer has to be set back to the first mcs
  //! Holds if first !mcs is read-in, already.
  bool first_call;

  /**
   * @brief Return if first !mcs is read-in, already.
   *
   * @return \a True if first call/occurance of !mcs - \a False otherwise.
   */
  bool isFirstCall() const {return first_call;}


  ///@todo better solution for solvent???
  //! Ignores if monomer position is missing in the !mcs (maybe due to solvent)
  bool ignore_missing_monos;

  //! The number of frames/conformations/!mcs-command in file.
  uint32_t nFrames;

 public:
   //constructor: sets ignore_missing_monos=true here!!
  ///@todo Why we ignore here monomers as default?
  //! Empty Constructor, but delegates call to the Feature. Default: ignoring of missing monomers.
  ReadJumpsMcs(IngredientsType& destination)
  :ReadToDestination< IngredientsType > (destination)
  ,ignore_missing_monos(true)
  ,first_call(true)
  ,nFrames(0){}

  //! Empty Destructor
  virtual ~ReadJumpsMcs(){}

  //enable ignoring of incomplete monomer initialzation
  /**
   * @brief Set this if missing monomers should ignore by read-in (maybe due to solvent).
   *
   * @details Specify the variable \a ignore_missing_monos if you need this (dangerous) behavior.
   * If not you get a \a runtime_error if monomers are missing.
   * @param val True if missing monomers should ignored - False otherwise.
   */
  void setIgnoreMissingMonomers(bool val){ignore_missing_monos = val;}

  //get information from file
  void execute();

};

/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !mcs when jumps are present (deprecated bfm-version).
 *
 * @throw <std::runtime_error> if Monte-Carlo-Step could not be parsed or if monomers are missing.
 **/
template < class IngredientsType > void ReadJumpsMcs< IngredientsType >::execute()
{
	  //get writeable reference to the molecules container
	  typename IngredientsType::molecules_type& molecules = this->getDestination().modifyMolecules();


	  if ( this->isFirstCall() ) {
		std::cout << "ReadMcs:execute() : updating ";
		std::cout << "connectivity and ";
		std::cout << "positions.";
	  }


	  //some variables needed for reading
	  int monomerCount=0;
	  int x,y,z,bondId;
	  unsigned long mcs;
	  std::string line;
	  std::streampos previous;

	  //read mcs number from file and update data
	  this->getInputStream()>>mcs;

	  molecules.setAge(mcs);

	  //throw exception if there is a reading error or update data otherwise
	  if(this->getInputStream().fail()){

	    std::stringstream errormessage;
	    errormessage<<"ReadMcs<IngredientsType>::execute()\n"
			<<"Could not read mcs number. Previous mcs number was "<<molecules.getAge();
	    throw std::runtime_error(errormessage.str());

	  }

	  molecules.setAge(mcs);

	  nFrames++;
	  if(nFrames%100 == 0) std::cout<<"Setting age: "<<mcs<<std::endl;
	  //move on to next line, where the positions begin.Save position of get pointer into "previous"
	  //this->getInputStream()->ignore(numeric_limits<streamsize>::max(), '\n');

	  // this line will only hold 'jumps' of the '!mcs'-Read if present (backward-compability)
	  // for the new LeMonADe-file there´re no jumps anymore!
	  std::string jumpString;
	  getline(this->getInputStream(),jumpString);

	  std::list < int32_t> JumpsInX;
	  std::list < int32_t> JumpsInY;
	  std::list < int32_t> JumpsInZ;

	  //tokenize the jump-string
	  std::vector<std::string> jumpStringVector;
	  jumpStringVector = this->tokenize2Parameter(jumpString, ' ', ':');

	  int32_t nPeriodConditions=0; //count the number of jumps if available

	  this->getDestination().isPeriodicX() ? nPeriodConditions++: 0;
	  this->getDestination().isPeriodicY() ? nPeriodConditions++: 0;
	  this->getDestination().isPeriodicZ() ? nPeriodConditions++: 0;

	  int32_t nChainStartsInJumpVector=0;
	  int32_t nChainStartsInJumpVectorX=0;
	  int32_t nChainStartsInJumpVectorY=0;
	  int32_t nChainStartsInJumpVectorZ=0;
	  // check if nr of jumps is correct
	  // only if jumps are present - this is due to the fact that jumps can be avoided in the old bfm-file
	  // we use the '<' operator instead of the '!=', maybe the file is corrupted
	  if(jumpStringVector.size() != 0)
	  {
		  //test if nChainBegins meet the periodic conditions
		  if ((jumpStringVector.size() % nPeriodConditions) != 0)
		      {
		      	  std::stringstream errormessage;
		      	  errormessage<<"ParsingError::ReadMcs<IngredientsType>::execute()\n At MCS:" << molecules.getAge()
	  			  			  <<"\n nElements:" << jumpStringVector.size() << " in jumps"
		      	  		      <<"\n Invalid number of chain jumps in old bfm-file\n Use the new standard!";

		      	  throw std::runtime_error(errormessage.str());
		      }

		  for (std::vector<std::string>::iterator iterjumpString=jumpStringVector.begin() ; iterjumpString < jumpStringVector.end();  )
		  {
			  if (this->getDestination().isPeriodicX() == true)
			  {
				  JumpsInX.push_back(atoi((*iterjumpString).c_str()));
				  //std::cout << "jump in x " << atoi((*iterjumpString).c_str()) << std::endl;
				  iterjumpString++;
				  nChainStartsInJumpVectorX++;
			  }

			  if (this->getDestination().isPeriodicY() == true)
			  {
				  JumpsInY.push_back(atoi((*iterjumpString).c_str()));
				  //std::cout << "jump in y " << atoi((*iterjumpString).c_str()) << std::endl;
				  iterjumpString++;
				  nChainStartsInJumpVectorY++;
			  }

			  if (this->getDestination().isPeriodicZ() == true)
			  {
				  JumpsInZ.push_back(atoi((*iterjumpString).c_str()));
				  //std::cout << "jump in z " << (*iterjumpString) << std::endl;
				  iterjumpString++;
				  nChainStartsInJumpVectorZ++;
			  }

			  nChainStartsInJumpVector++;

		  }
	  }


	  //go on with the all positions etc.
	  previous=this->getInputStream().tellg();
	  getline(this->getInputStream(),line); //get the first line with positions


	  //std::cout<<line<< std::endl;

	  //process input lines in this loop
	  while(!line.empty() && !this->getInputStream().fail()){

	    //if the line contains a bfm Read, stop the procedure and set the get pointer back
	    if(this->detectRead(line)){
	      this->getInputStream().seekg(previous);
	      return;
	    }

	    //initialize stringstream with content for ease of processing
	    std::stringstream stream(line);

	    //read first three coordinates
	    stream>>x>>y>>z;

	    //add/apply jumps if possible
	    if ((this->getDestination().isPeriodicX() == true) && (JumpsInX.size() != 0))
	    {
	    	x += this->getDestination().getBoxX()*JumpsInX.front();
	    	JumpsInX.pop_front(); // erase the first element
	    }

	    if ((this->getDestination().isPeriodicY() == true) && (JumpsInY.size() != 0))
	    {
	    	y += this->getDestination().getBoxY()*JumpsInY.front();
	    	JumpsInY.pop_front(); // erase the first element
	    }

	    if ((this->getDestination().isPeriodicZ() == true) && (JumpsInZ.size() != 0))
	    {
	    	z += this->getDestination().getBoxZ()*JumpsInZ.front();
	    	JumpsInZ.pop_front(); // erase the first element
	    }

	    if(!stream.fail())
	    {
	    	molecules[monomerCount].setAllCoordinates(x,y,z);
			++monomerCount;
	    }

	    //throw exception if first coordinates of chain cannot be extracted from file
	    if(stream.fail())
	    {

	      std::stringstream errormessage;
	      errormessage<<"ParsingError::ReadMcs<IngredientsType>::execute()\n"
			  <<"Could not read chain initial coordinates in mcs "<<molecules.getAge();

	      throw std::runtime_error(errormessage.str());

	    }

	    //ignore spaces
	    stream.ignore(1);

	      //read the ASCII coded bond vectors of this chain
	    while(stream.peek()>0 && stream.good())
	    {
	      //read next bond from file
	      bondId=stream.get();
	      x+=this->getDestination().getBondset().getBondVector(bondId).getX();
	      y+=this->getDestination().getBondset().getBondVector(bondId).getY();
	      z+=this->getDestination().getBondset().getBondVector(bondId).getZ();

	      //update molecules
	      molecules[monomerCount].setAllCoordinates(x,y,z);


	      //set up connections at first monte carlo step
	      if ( this->isFirstCall())
	      {
	    	 // std::cout<<"monomer at \t"<< monomerCount << " " <<x<<" "<<y<<" "<<z<<std::endl;

	    	  molecules.connect(monomerCount-1,monomerCount);

	      }

	      //count number of monomers in this mcs
	      ++monomerCount;
	    }

	    //get next line from file
	    getline(this->getInputStream(),line);
	  }

	  //molecules.resize(monomerCount);

	  //if total number of monomers in !mcs differs from previous, throw exception
	  if(monomerCount!=molecules.size()){

	    std::stringstream errormessage;
	    errormessage<<"LogicError::ReadMcs<IngredientsType>::execute()\n"
			<<"Inconsistent number of monomers in bfm file in mcs number "
			<<molecules.getAge()
			<<"\nand chain beginning with coordinates "<<x<<" "<<y<<" "<<z;
	    if ( (this->ignore_missing_monos==true)  || this->isFirstCall() )
	    {
		if(nFrames%100==0)
		{
		std::cout << "Ignoring logic error:\n" << errormessage.str() << std::endl;
		std::cout << "Not all monomers have been initialized.\n" << std::endl;
		}
	    }
	    else
	    {
		throw std::runtime_error(errormessage.str());
	    }

	  }


	  //check if nr of jumps is correct

	  //count nr of chain starts
	  int chainStarts = 1;
	  for (int i=1; i<molecules.size(); ++i){
		  if (!molecules.areConnected(i,i-1)){
			  chainStarts++;
		  }
	  }

	  // only if jumps are present - this is due to the fact that jumps can be avoided in the old bfm-file
	  if ((nChainStartsInJumpVector > 0) && ((nChainStartsInJumpVector != chainStarts) || ((nChainStartsInJumpVectorX+nChainStartsInJumpVectorY+nChainStartsInJumpVectorZ) != nPeriodConditions*chainStarts) ))
	  {
	  	  std::stringstream errormessage;
	  	  errormessage<<"ParsingError::ReadMcs<IngredientsType>::execute()\n at MCS:" << molecules.getAge()
	  			  	  <<"\n single elements:" << (nChainStartsInJumpVectorX+nChainStartsInJumpVectorY+nChainStartsInJumpVectorZ) << " jumps"
	  	  		      <<"\nInvalid number of chain jumps in old bfm-file\n Use the new standard!\n";

	  	  throw std::runtime_error(errormessage.str());
	  }


	  first_call = false;


	}

/***********************************************************************/
/**
 * @class WriteJumpsMcs
 * @brief Handles BFM-File-Write \b !mcs when jumps are present (deprecated bfm-version).
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteJumpsMcs: public AbstractWrite<IngredientsType>
{
public:
  WriteJumpsMcs(const IngredientsType& src):AbstractWrite<IngredientsType>(src){};
  void writeStream(std::ostream& strm);
};

/***********************************************************************/
//! Executes the routine to write \b !mcs when jumps are present (deprecated bfm-version).
template <class IngredientsType>
void WriteJumpsMcs<IngredientsType>::writeStream(std::ostream& strm){

    const typename IngredientsType::molecules_type& molecules(this->getSource().getMolecules());

    std::stringstream contents;

    contents<<"\n!mcs="<< molecules.getAge();

    //count nr of chain starts
    unsigned int chainStarts = 1;
    std::vector<unsigned int> chainStartsVector;

    chainStartsVector.push_back(0); //first monomer is always chain start

    for (int i=1; i<molecules.size(); ++i)
    {
    	if (!molecules.areConnected(i,i-1))
    	{
    		chainStarts++;
    		chainStartsVector.push_back(i);
    	}
    }

	for (unsigned int counterChainStart = 0; counterChainStart < chainStarts; counterChainStart++) {
		if ( ((this->getSource()).isPeriodicX() == true) || ((this->getSource()).isPeriodicY() == true) || ((this->getSource()).isPeriodicZ()	== true))
			contents << ":";

		if ((this->getSource()).isPeriodicX())
			contents << (molecules[chainStartsVector.at(counterChainStart)].getX()-((molecules[chainStartsVector.at(counterChainStart)].getX()%(this->getSource()).getBoxX())+(this->getSource()).getBoxX())%(this->getSource()).getBoxX())/(this->getSource()).getBoxX();

		if (((this->getSource()).isPeriodicX()) && (((this->getSource()).isPeriodicY()) || ((this->getSource()).isPeriodicZ())))
			contents << " ";


		if ((this->getSource()).isPeriodicY())
			contents << (molecules[chainStartsVector.at(counterChainStart)].getY()-((molecules[chainStartsVector.at(counterChainStart)].getY()%(this->getSource()).getBoxY())+(this->getSource()).getBoxY())%(this->getSource()).getBoxY())/(this->getSource()).getBoxY();


		if (((this->getSource()).isPeriodicY()) && ((this->getSource()).isPeriodicZ()))
			contents << " ";


		if ((this->getSource()).isPeriodicZ())
			contents << (molecules[chainStartsVector.at(counterChainStart)].getZ()-((molecules[chainStartsVector.at(counterChainStart)].getZ()%(this->getSource()).getBoxZ())+(this->getSource()).getBoxZ())%(this->getSource()).getBoxZ())/(this->getSource()).getBoxZ();

	}

    //write position of first monomer
	//c++ (ISO 1998)
	//there is a problem with the modulo-operator with negative numbers (april 2014)
	//in ISO 1998 5.6/4 (C++03), the sign is implementation defined :(
	//in ISO (2011) C++11 the sign has the dividend sign, but up to now we using C++03
    contents<<"\n"<<((molecules[0].getX()%(this->getSource()).getBoxX())+(this->getSource()).getBoxX())%(this->getSource()).getBoxX();
    contents<<" " <<((molecules[0].getY()%(this->getSource()).getBoxY())+(this->getSource()).getBoxY())%(this->getSource()).getBoxY();
    contents<<" " <<((molecules[0].getZ()%(this->getSource()).getBoxZ())+(this->getSource()).getBoxZ())%(this->getSource()).getBoxZ();

    //insert space after first monomer
    contents<<" ";

    //write bonds and other monomers
    for(size_t n=1;n< molecules.size();n++){

    	//single unconnected monomer
    	if (molecules.getNumLinks(n) == 0) {
    		contents<<"\n"<<((molecules[n].getX()%(this->getSource()).getBoxX())+(this->getSource()).getBoxX())%(this->getSource()).getBoxX();
    		contents<<" " <<((molecules[n].getY()%(this->getSource()).getBoxY())+(this->getSource()).getBoxY())%(this->getSource()).getBoxY();
    		contents<<" " <<((molecules[n].getZ()%(this->getSource()).getBoxZ())+(this->getSource()).getBoxZ())%(this->getSource()).getBoxZ();
    		contents<<" ";

    		}

      //if the actual monomer is a chainstart, start a new subchain in !mcs Read
    	else if (!molecules.areConnected(n,n-1)){
    		contents<<"\n"<<((molecules[n].getX()%(this->getSource()).getBoxX())+(this->getSource()).getBoxX())%(this->getSource()).getBoxX();
    		contents<<" " <<((molecules[n].getY()%(this->getSource()).getBoxY())+(this->getSource()).getBoxY())%(this->getSource()).getBoxY();
    		contents<<" " <<((molecules[n].getZ()%(this->getSource()).getBoxZ())+(this->getSource()).getBoxZ())%(this->getSource()).getBoxZ();
    		contents<<" ";

      }

      //otherwise write the bond
      else{
    	  contents<<char(this->getSource().getBondset().getBondIdentifier(
		(molecules[n].getX())-(molecules[n-1].getX()),
		(molecules[n].getY())-(molecules[n-1].getY()),
		(molecules[n].getZ())-(molecules[n-1].getZ()) ) );

    	  contents.flush();

      }
    }

    contents<<"\n\n";
    contents.flush();

    strm << contents.str();
    strm.flush();
}



/***********************************************************************/
/**
 * @class FeatureJumps
 * @brief Enables support for !mcs with jump configuration (deprecated bfm-version).
 *
 * @details This FeatureJumps replaces the ReadMCS in MoleculesRead!
 * 		  This FeatureJumps replaces the WriteMCS in MoleculesWrite!
 **/
class FeatureJumps:public Feature
{
public:

	//! This Feature requires a box.
	typedef LOKI_TYPELIST_1(FeatureBox) required_features_front;

/**
 * @brief Export the relevant functionality for reading bfm-files to the responsible reader object
 *
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !mcs with jumps (deprecated bfm-version).
 *
 * @param fileReader File importer for the bfm-file
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
  template <class IngredientsType>
  void exportRead(FileImport <IngredientsType>& fileReader)
  {
    fileReader.replaceRead("!mcs",new ReadJumpsMcs <IngredientsType> (fileReader.getDestination()));
  }

  /**
   * @brief Export the relevant functionality for writing bfm-files to the responsible writer object
   *
   * @details The function is called by the Ingredients class when an object of type Ingredients
   * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
   * taken care automatically when it becomes necessary.\n
   * Registered Write-Out Commands:
   * * !mcs with jumps (deprecated bfm-version).
   *
   * @param fileWriter File writer for the bfm-file.
   * @tparam IngredientsType Features used in the system. See Ingredients.
   */
  template <class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const
  {
    fileWriter.replaceWrite("!mcs",new WriteJumpsMcs <IngredientsType> (fileWriter.getIngredients_()));
  }

};

#endif /* LEMONADE_FEATURE_FEATUREJUMPS_H */
