

#include <cstring>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/feature/FeatureNNInteraction.h>

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
		errormessage="usage: ./SimpleSimulator input_file max_mcs save_interval(mcs) [output_file]\n";
		errormessage+="\nSimple Simulator for the ScBFM with Ex.Vol and BondCheck\n";
		errormessage+="maximum number of connections per monomer is 6\n";
		errormessage+="If output_filename specified, the results are written to the new file\n";
		errormessage+="otherwise the results are appended to the old input file\n";
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
	
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureNNInteractionSc<FeatureLatticePowerOfTwo>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	//taskmanager.addAnalyzer(new LatticeOccupationAnalyzer<Ing>(myIngredients));
	
	taskmanager.initialize();
	taskmanager.run(max_mcs/save_interval);
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

