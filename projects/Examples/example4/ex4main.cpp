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
* This example demonstrates the basic use of updaters and analyzers.
* A simple simulation program is put together and can be executed.
* In the example, a short polymer chain is placed in the system and
* simulated in a periodic box for 1E7 Monte Carlo steps. The
* configurations are saved to a file so they can be watched later, and
* the radius of gyration is analyzed.
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



int main(int, char**)
{
    /* ****************************************************************
      * as described in the previous examples, we quickly seed the
      * random number generator, and set up a system in a box with the
      * standard BFM bond set. We add 10 monomers and connect them to a
      * linear polymer.
      * ***************************************************************/

    //first set up the random number generator
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    //now set up the system
    //here, we use one additional feature called FeatureMoleculesIO
    //strictly speaking this is not a feature in the sense that it
    //changes the simulation conditions. What it provides is the
    //basic functionalities for writing BFM files. So, in most cases
    //you are going to use this feature in simulations.
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


    //and add a polymer chain of length ten
    uint32_t polymerLength=10;

    mySystem.modifyMolecules().resize(polymerLength);

    //add the polymer chain
    mySystem.modifyMolecules()[0].setX(0);
    mySystem.modifyMolecules()[0].setY(0);
    mySystem.modifyMolecules()[0].setZ(0);

    for(uint32_t i=1;i<polymerLength;i++){
	    mySystem.modifyMolecules()[i].setX(i*2);
	    mySystem.modifyMolecules()[i].setY(i*2);
	    mySystem.modifyMolecules()[i].setZ(i);

	    mySystem.modifyMolecules().connect(i-1,i);
    }


    mySystem.synchronize(mySystem);

    /* ****************************************************************
      * Now we can set up the task manager with the desired tasks.
      * We want to do three things:
      * 1) Simulate the system. For this we use the updater
      * UpdaterSimpleSimulator. The source code of this updater
      * can be found in src/updater/UpdaterSimpleSimulator.h
      * 2) Write the trajectory to a file. For this we use the analyzer
      * AnalyzerWriteBfmFile. The source code of this analyzer can be
      * found in src/analyzer/AnalyzerWriteBfmFile.h
      * 3) Calculate the radius of gyration. For this we use the
      * analyzer AnalyzerRadiusOfGyration, which  can be found in
      * src/analyzer/AnalyzerRadiusOfGyration.h
      *
      * For adding updaters and analyzers, the task manager provides
      * the functions addUpdater, addAnalyzer. These functions take two
      * arguments:
      * 1) a pointer to the updater/analyzer
      * 2) a number indicating how often the particular task should be
      * executed. For example, 1 means every cycle, 2 every second
      * cycle, etc. Setting the second argument to 0 means the task
      * will be executed only once in the first cycle. This is useful
      * if for example a polymer is added in the beginning by an updater.
      * Here, we have created the polymer by hand, so this is not
      * necessary.
      * ***************************************************************/

    //create the task manager
    TaskManager taskmanager;

    //add the simulator
    //the syntax says: the simulator should simulate mySystem
    //every time it is executed, it simulates for 10000 steps
    //The 1 as second argument to addUpdater says that the
    //simulation is to be called in every cycle.
    taskmanager.addUpdater(new
    UpdaterSimpleSimulator<MyIngredients,MoveLocalSc>(mySystem,10000),1);

    //add the file output, the trajectory file name will be "polymer.bfm"
    taskmanager.addAnalyzer(new
    AnalyzerWriteBfmFile<MyIngredients>("polymer.bfm",mySystem),1);

    //add the radius of gyration calculation. we execute this only
    //every 10th time
    //the first argument to the analyzer is the system to be analyzed,
    //the second is a string defining the name of the output file.
    taskmanager.addAnalyzer(new
    AnalyzerRadiusOfGyration<MyIngredients>(mySystem,"Rg2.dat"),10);

    /* ****************************************************************
      * For running the desired tasks, the task manager provides four
      * functions:
      * - initialize() : should always be called at the beginning. It
      * initializes the updaters and analyzers. More on this topic in
      * the next example of the tutorial.
      * - run() : runs the circle of chosen tasks until all updaters
      * are finished
      * - run(nCircles) : runs the circle of chosen tasks nCircles times
      * - cleanup() : should be called at the end. finishes the tasks
      * of updaters and analyzers. More on this topic in the next
      * example of the tutorial.
      * ***************************************************************/

    //this will prepare and run the simulation. look in the directory
    //where you called the program for output files.
    //we run the task manager 1000 times, because in every run, the
    //simulator makes 10000 steps->total of 1E7 steps.
    //when you run this example program, you will see some output
    //indicating the progress of the simulation on the screen.

    taskmanager.initialize();
    taskmanager.run(1000);
    taskmanager.cleanup();

    /* **************************************************************
      * The program should produce among others a file called
      * polymer.bfm. This contains the trajectory. You can visualize
      * the trajectory using the provided LemonadeViewerFLTK from
      * the projects folder.
      * *************************************************************/

    return 0;
}
