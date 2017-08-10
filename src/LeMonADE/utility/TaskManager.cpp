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

#include <LeMonADE/utility/TaskManager.h>

/*****************************************************************************/
/**
 * @file
 * @brief Implementation TaskManager members
 * */
/*****************************************************************************/

//constructor and destructor
TaskManager::TaskManager():running(true),nCircles(0){}

TaskManager::~TaskManager()
{
    vector <UpdaterObject*>::iterator uIterator(updater.begin());
    vector <AnalyzerObject*>::iterator aIterator(analyzer.begin());

    //free memory allocated for UpdaterObjects
    while(uIterator!=updater.end()){
      delete *uIterator;
      *uIterator=0;
      ++uIterator;
    }
    //free memory allocated for AnalyzerObjects
    while(aIterator!=analyzer.end()){
      delete *aIterator;
      *aIterator=0;
      ++aIterator;
    }
}

/*****************************************************************************/
/**
 * @param a Pointer to new AnalyzerObject
 * @param period execution period for the new analyzer. If 0 is used, the analyzer
 *        executes only once at the beginning. Default is 1 - every cycle.
 **/
void TaskManager::addAnalyzer(AbstractAnalyzer* a, int period)
{
  analyzer.push_back(new AnalyzerObject(a,period));
}

/*****************************************************************************/
/**
 * @param a Pointer to new UpdaterObject
 * @param period execution period for the new updater. If 0 is used, the updater
 *        executes only once at the beginning. Default is 1 - every cycle.
 **/
void TaskManager::addUpdater(AbstractUpdater* u, int period)
{
  updater.push_back(new UpdaterObject(u,period));
}

/*****************************************************************************/
/**
 * @details The function periodically executes first all updaters and then all analyzers.
 * The loop stops when the run() functions of ALL updaters return false. In case
 * there are no updaters, the analyzers are executed once.
 **/
void TaskManager::run()
{
  vector <UpdaterObject*>::iterator uIterator;
  vector <AnalyzerObject*>::iterator aIterator;

  //if there are no updaters, execute all analyzers and return
  if(updater.size()==0){
    aIterator=analyzer.begin();
    while(aIterator!=analyzer.end() ){
      if((*aIterator)->shouldExecute(nCircles)) (*aIterator)->run();
      ++aIterator;
    }
    return;
  }
  //else run as long as there are updates coming in from updater objects
  else{
    while(true){
      ++nCircles;
      running=false;

      //call all updaters
      uIterator=updater.begin();
      while(uIterator!=updater.end() ){
	if((*uIterator)->shouldExecute(nCircles)){
		bool currentUpdaterRunning=(*uIterator)->run();
		running=(running ||currentUpdaterRunning);

	}
	++uIterator;
      }
      //if all updaters return false, exit the loop
      if(!running) break;

      //call all analyzers
      aIterator=analyzer.begin();
      while(aIterator!=analyzer.end() ){
	if((*aIterator)->shouldExecute(nCircles)) (*aIterator)->run();
	++aIterator;
      }

    }
  }



}

/*****************************************************************************/
/**
 * @details The function periodically executes first all updaters and then all analyzers.
 * The loop stops after nPeriods circles.
 *
 * @param nPeriods number of execution circles
 **/
void TaskManager::run(int nPeriods)
{

  vector <UpdaterObject*>::iterator uIterator;
  vector <AnalyzerObject*>::iterator aIterator;
  int actualCircles=nCircles;
  while(nCircles<nPeriods+actualCircles){
    ++nCircles;

    //call all updaters
    uIterator=updater.begin();
    ///@todo exception, if all updaters return false and loop keeps running??
    while(uIterator!=updater.end() ){
      if((*uIterator)->shouldExecute(nCircles)) (*uIterator)->run();
      ++uIterator;
    }

    //call all analyzers
    aIterator=analyzer.begin();
    while(aIterator!=analyzer.end() ){
      if((*aIterator)->shouldExecute(nCircles)) (*aIterator)->run();
      ++aIterator;
    }

  }

}

/*****************************************************************************/
/**
 * @todo initialize updaters necessary?
 **/
void TaskManager::initialize()
{
    vector <UpdaterObject*>::iterator uIterator;
    vector <AnalyzerObject*>::iterator aIterator;

    uIterator=updater.begin();
    while(uIterator!=updater.end()){
      (*uIterator)->initialize();
      ++uIterator;
    }

    aIterator=analyzer.begin();
    while(aIterator!=analyzer.end() ){
      (*aIterator)->initialize();
      ++aIterator;
    }
}

/*****************************************************************************/
/**
 * @todo cleanup updaters necessary?
 **/
void TaskManager::cleanup()
{
    vector <UpdaterObject*>::iterator uIterator;
    vector <AnalyzerObject*>::iterator aIterator;

    uIterator=updater.begin();
    while(uIterator!=updater.end()){
      (*uIterator)->cleanup();
      ++uIterator;
    }

    aIterator=analyzer.begin();
    while(aIterator!=analyzer.end() ){
      (*aIterator)->cleanup();
      ++aIterator;
    }
}




