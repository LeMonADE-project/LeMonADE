

#include <cstring>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/Version.h>

int main(int argc, char* argv[])
{
  try{
	std::string infile;
	std::string outfile;
	uint32_t max_mcs=0;
	uint32_t save_interval=0;

	outfile="outfile.bfm";

	if(!(argc==4 || argc==5 )|| (argc==2 && strcmp(argv[1],"--help")==0 ))
	{
		std::string errormessage;
		errormessage="usage: ./SimpleSimulator input_filename max_mcs save_interval(mcs) [output_filename]\n";
		errormessage+="\nSimple Simulator for the ScBFM with Ex.Vol and BondCheck\n";
		errormessage+="maximum number of connections per monomer is 6\n";
		errormessage+="If output_filename specified, the results are written to the new file\n";
		errormessage+="otherwise the results are appended to the old input file\n";
		errormessage+="Features used: FeatureFixedMonomers, FeatureBondset,FeatureAttributes,FeatureExcludedVolume<FeatureLattice<bool> >\n";
		errormessage+="Updaters used: ReadFullBFMFile, SimpleSimulator\n";
		errormessage+="Analyzers used: WriteBfmFile\n";
		throw std::runtime_error(errormessage);

	}
	else
	{
		infile=argv[1];
		max_mcs=atoi(argv[2]);
		save_interval=atoi(argv[3]);

		if(argc==5) outfile=argv[4];
		else outfile=argv[1];
	}
	  
	std::cout << APPLICATION_NAME " is free software: you can redistribute it and/or modify" << std::endl
              << "it under the terms of the GNU General Public License as published by" << std::endl
              << "the Free Software Foundation, either version 3 of the License, or" << std::endl
              << "(at your option) any later version." << std::endl
              << APPLICATION_NAME " is distributed in the hope that it will be useful," << std::endl
              << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl
              << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl
              << "GNU General Public License for more details." << std::endl
              << "You should have received a copy of the GNU General Public License" << std::endl
              << "along with " << APPLICATION_NAME << ". If not, see <http://www.gnu.org/licenses/>." << std::endl << std::endl;
             
    	std::cout << APPLICATION_NAME << " " << APPLICATION_VERSION_STRING << std::endl
              << "Copyright (C) " << APPLICATION_COPYRIGHT_YEARS << std::endl
              << "LeMonADE Version " <<  LEMONADE_VERSION << std::endl << std::endl;


	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();

	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;

	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));

	taskmanager.initialize();
	taskmanager.run(max_mcs/save_interval);
	taskmanager.cleanup();

	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;

}

