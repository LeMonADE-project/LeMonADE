/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Hauke Rabbel)
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

/* *********************************************************************
* This example demonstrates how you can write your own updaters and
* analyzers and use them together with the task manager. The example
* program is not much different to the previous one. The difference here
* is that we use an updater to create particles, and we write an analyzer
* by ourself.
* The example updater and analyzer are in the files
* ex5updater.h and ex5analyzer.h
*
* After finishing this example tutorial, take a look into the analyzers
* and updaters provided in src/updaters and src/analyzers. You now
* have the tools to understand their structure and how they work.
* *********************************************************************/

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

#include "ex5analyzer.h"
#include "ex5updater.h"



int main(int, char**)
{
    /* ****************************************************************
      * the first part looks exactly like in the last example.
      * ***************************************************************/

    //first set up the random number generator
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    //now set up the system
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO,FeatureBox,FeatureBondset<>) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> MyIngredients;
    MyIngredients mySystem;

    mySystem.setBoxX(64);
    mySystem.setBoxY(64);
    mySystem.setBoxZ(64);

    mySystem.setPeriodicX(true);
    mySystem.setPeriodicY(true);
    mySystem.setPeriodicZ(true);

    mySystem.modifyBondset().addBFMclassicBondset();

    mySystem.synchronize(mySystem);


    /* ****************************************************************
      * as compared to the previous example, we now don't create the
      * particles by hand, but we use an updater for this in the
      * task manager. the advantage of this approach is that you can
      * reuse this updater in different programs, so you don't have to
      * write the code for creating particles every time.
      *
      * we want the particles to be created once before the simulation.
      * so we have to tell the task manager this when adding the
      * updater. this is done by giving the number 0 as execution
      * frequency, as shown below.
      * the largest part of the code below is exactly the same as in the
      * previous example.
      * ***************************************************************/

    //create the task manager
    TaskManager taskmanager;

    //add the ex5updater, which creates the particles
    //again, the addUpdater function takes two arguments. The first one
    //is a pointer to the updater, and the second one is the execution
    //frequency. here we use 0 as frequency, because we want the updater
    //to be executed only once at the beginning.
    //what would happen, if you set a number different from 0?
    uint32_t numberOfParticles=100;
    taskmanager.addUpdater(new Ex5Updater<MyIngredients>(mySystem,numberOfParticles),0);


    //add the simulator
    //here we set the frequency to one, because we want to execute it in every
    //cycle
    taskmanager.addUpdater(new
    UpdaterSimpleSimulator<MyIngredients,MoveLocalSc>(mySystem,1000),1);

    //the following two analyzers are the same as in the last example:
    //on top of this, we add a third one, which is the ex5analyzer
    //written for this example.

    //add the file output, the trajectory file name will be "polymer.bfm"
    taskmanager.addAnalyzer(new
    AnalyzerWriteBfmFile<MyIngredients>("polymer.bfm",mySystem),1);

    //add the radius of gyration calculation. we execute this only
    //every 10th time
    taskmanager.addAnalyzer(new
    AnalyzerRadiusOfGyration<MyIngredients>(mySystem,"Rg2.dat"),10);

    //now add the ex5analyzer. execute it every second cycle
    std::string outputFilename("ex5output.dat");
    taskmanager.addAnalyzer(new Ex5Analyzer<MyIngredients>(mySystem,outputFilename),2);

    //as before, we initialize and run the task manager, and cleanup afterwards
    taskmanager.initialize();
    taskmanager.run(10000);
    taskmanager.cleanup();

    return 0;
}
