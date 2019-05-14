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

#ifndef LEMONADE_UTILITY_TASKMANAGER_H
#define LEMONADE_UTILITY_TASKMANAGER_H

/*****************************************************************************/
/**
 * @file
 * @brief Definitions of class TaskManager
 **/
/*****************************************************************************/

#include <iostream>
#include <list>
#include <vector>

#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/updater/AbstractUpdater.h>

using std::vector;

/*****************************************************************************/
/**
 * @class TaskManager
 *
 * @brief Manages Updaters and Analyzers
 *
 * @todo rename execution routine
 *
 * @todo use pure virtual routine of AbstractAnalyzer and AbstractUpdater
 *
 * @todo reimplement the add-routines avoiding direct pointer
 *
 * @todo add synchronize between init, execute, cleanup
 **/
class TaskManager
{
public:
  TaskManager();
  ~TaskManager();

  //! Add an analyzer to the list of tasks
  void addAnalyzer(AbstractAnalyzer*, int period=1);

  //! Add an updater to the list of tasks
  void addUpdater(AbstractUpdater*, int period=1);

  int getNCircles(){return nCircles;}

  //! Start execution circle of updaters and analyzers
  void run();

  //! Start execution circle of updaters and analyzers with given number of cycles.
  void run(int nPeriods);

  //! Calls initialize() routine on all updater- and analyzer-objects
  void initialize();

  //! Calles cleanup() routine on all updater- and analyzer-objects
  void cleanup();

private:

  //! Holds an analyzer and executes it with period execution_period
  class AnalyzerObject;

  //! Holds an updater and executes it with period execution_period
  class UpdaterObject;



  vector <AnalyzerObject*> analyzer;
  vector <UpdaterObject*> updater;

  //! Running criterion for running or stopping the execution loop
  bool running;

  //! nCircles counter for execution loops
  int nCircles;
};

/*****************************************************************************/
/**
 * @class TaskManager::AnalyzerObject
 *
 * @brief Holds an analyzer and executes it with period execution_period
 *
 * @details The execution period (in units of TaskManager::nCircles) is specified
 * as constructor argument. If 0 is used, the analyzer executes only once at
 * the beginning.
 **/
class TaskManager::AnalyzerObject{
public:
  AnalyzerObject(AbstractAnalyzer* a,unsigned int period)
    :myAnalyzer(a),execution_period(period),isFirstExecution(true){};

  ~AnalyzerObject(){delete myAnalyzer;};

  /**
   * @brief Executes the analyzer
   *
   * @todo rename to execute
   **/
  void run(){ myAnalyzer->execute();};

  //! Calls the analyzer's initialize routine
  void initialize(){myAnalyzer->initialize();}

  //! Calls the analyzer's cleanup routine
  void cleanup(){myAnalyzer->cleanup();}

  //! Decides whether or not to execute in execution circle given as argument
  bool shouldExecute(int nCircles){
    if(execution_period!=0) return (nCircles%execution_period==0);
    else{
      if(isFirstExecution) {isFirstExecution=false; return true;}
      else return false;
    }
  };
private:
  AbstractAnalyzer* myAnalyzer;
  unsigned int execution_period;
  bool isFirstExecution;
};

/*****************************************************************************/
/**
 * @class TaskManager::UpdaterObject
 *
 * @brief Holds an updater and executes it with period execution_period
 *
 * @details The execution period (in units of TaskManager::nCircles) is specified
 * as constructor argument. If 0 is used, the updater executes only once at
 * the beginning.
 **/
class TaskManager::UpdaterObject{
public:
  UpdaterObject(AbstractUpdater* u, unsigned int period)
    :myUpdater(u),execution_period(period),isFirstExecution(true){};

  ~UpdaterObject(){delete myUpdater;};

  /**
   * @brief Executes the updater
   *
   * @todo rename to execute
   **/
  bool run(){return myUpdater->execute();}

  //! Calls the updater's initialize routine
  void initialize(){myUpdater->initialize();}

  //! Calls the updater's cleanup routine
  void cleanup(){myUpdater->cleanup();}

  //! Decides whether or not to execute in execution circle given as argument
  bool shouldExecute(int nCircles){
    if(execution_period!=0) return (nCircles%execution_period==0);
    else{
      if(isFirstExecution) {isFirstExecution=false; return true;}
      else return false;
    }
  }

private:
  AbstractUpdater* myUpdater;
  unsigned int execution_period;
  bool isFirstExecution;
};

#endif /* LEMONADE_UTILITY_TASKMANAGER_H */
