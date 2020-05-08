/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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


/*********************************************************************
 * written by      : Toni Müller
 * email           : mueller-toni@ipfdd.de
 *********************************************************************/

#include <iostream>
#include <exception>

#include <boost/program_options.hpp> // for the command line handling / options
using namespace boost::program_options;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>
#include <LeMonADE/analyzer/AnalyzerMonomerMSD.h>
#include <LeMonADE/analyzer/AnalyzerSystemMSD.h>

int main (int argc, char* argv[])
{
  //read arguments from the command line
  std::string input,output;
  try
  {
    options_description desc{"Analyze an equilibrated melt of chains in the BFM file format.\nAllowed options"};
    desc.add_options()
      ("help,h"      ,                                                          "produce help message")
      ("input,i"     , value<std::string>(&input)->default_value("config.bfm"), "input"               )
      ("output,o"    , value<std::string>(&output)->default_value("output")   , "prefix for the output"              );
    variables_map options_map;
    store(parse_command_line(argc, argv, desc), options_map);
    notify(options_map); 
    
    // help option
    if (options_map.count("help")) 
    {
      std::cout << desc << "\n";
      return 1;
    }
  }
  catch(std::exception& e ){std::cerr<<"Error:\n" << e.what()<< std::endl;}
  catch(...){std::cerr<<"Error: unknown exception\n";}
  typedef LOKI_TYPELIST_1(FeatureMoleculesIO) Features; 
  typedef ConfigureSystem<VectorInt3,Features, 7> Config;
  typedef Ingredients<Config> Ing;
  Ing myIngredients;
  
  //set up the sytem parameters
  TaskManager taskmanager;
  taskmanager.addUpdater ( new UpdaterReadBfmFile      <Ing>(input,myIngredients, UpdaterReadBfmFile<Ing>::READ_STEPWISE) );
  taskmanager.addAnalyzer( new AnalyzerMonomerMSD      <Ing>(myIngredients, 0, output) ) ;
  taskmanager.addAnalyzer( new AnalyzerSystemMSD       <Ing>(myIngredients, 0, output) ) ;
  
  taskmanager.initialize();
  taskmanager.run();
  taskmanager.cleanup();
  
  return 0;
}
