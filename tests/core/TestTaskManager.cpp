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

/*
 * TestTaskManager.cpp
 *
 *  Created on: Jul 18, 2013
 *      Author: Hauke Rabbel
 */
#include "gtest/gtest.h"

#include <sstream>
#include <string>

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/TaskManager.h>

using namespace std;
/*
 * for the following tests a few dummy analyzers and updaters are needed.
 * these are defined in the following section
 */

class DummyUpdater:public AbstractUpdater
{
public:
	DummyUpdater(string& dest):counter(0),destination(&dest),isInitialized(false),isCleanedUp(false){}
	bool execute(){
		destination->append("1");
		counter++;
		if(counter<5) return true;
		else return false;
	}

	void initialize(){isInitialized=true;}
	void cleanup(){isCleanedUp=true;}

	bool isInitialized;
	bool isCleanedUp;
private:
	int counter;
	string* destination;


};

class DummyAnalyzer: public AbstractAnalyzer
{
public:
	DummyAnalyzer(string& str):source(&str),isInitialized(false),isCleanedUp(false){}
	bool execute(){
		result.append(*source);
		return true;
	}

	void initialize(){isInitialized=true;}
	void cleanup(){isCleanedUp=true;}

	string result;
	bool isInitialized;
	bool isCleanedUp;
private:
	string* source;

};


/*
 * this test checks if new updaters and analyzers can be added correctly.
 * in particular it is checked if the execution periods are set and work.
 */

TEST(TaskManagerTest,AddUpdaters)
{
	string textstring;
	TaskManager taskmanager;
	taskmanager.addUpdater(new DummyUpdater(textstring));
	taskmanager.addUpdater(new DummyUpdater(textstring),2);
	taskmanager.run(4);
	//since the second updater executes only every second time,
	//the textstring should now contain 6 letters.
	EXPECT_EQ(textstring.size(),6);

}

TEST(TaskManagerTest,AddAnalyzers)
{
	string textstring;
	TaskManager taskmanager;
	DummyAnalyzer* analyzer1=new DummyAnalyzer(textstring);
	DummyAnalyzer* analyzer2=new DummyAnalyzer(textstring);
	taskmanager.addUpdater(new DummyUpdater(textstring));
	taskmanager.addAnalyzer(analyzer1);
	taskmanager.addAnalyzer(analyzer2,2);
	taskmanager.run(4);
	//since the first analyzer executes every time, its string "result"
	//should now contain 10 characters (1+2+3+4=10).
	//similarily, since the second analyzer executes only every second time,
	//starting with the second cirlce, its result string should contain
	//6 characters (2+4=6)
	EXPECT_EQ(analyzer1->result.size(),10);
	EXPECT_EQ(analyzer2->result.size(),6);

}

TEST(TaskManagerTest, Run)
{
  string textstring1,textstring2;
  TaskManager taskmanager1,taskmanager2;
  DummyUpdater* updater1=new DummyUpdater(textstring1);
  DummyUpdater* updater2=new DummyUpdater(textstring2);
  DummyAnalyzer* analyzer1=new DummyAnalyzer(textstring1);
  DummyAnalyzer* analyzer2=new DummyAnalyzer(textstring2);

  taskmanager1.addUpdater(updater1);
  taskmanager1.addAnalyzer(analyzer1);
  taskmanager2.addUpdater(updater2);
  taskmanager2.addAnalyzer(analyzer2);

  taskmanager1.run();
  taskmanager2.run(2);

  EXPECT_EQ(taskmanager1.getNCircles(),5);
  EXPECT_EQ(taskmanager2.getNCircles(),2);
  taskmanager2.run(1);
  EXPECT_EQ(taskmanager2.getNCircles(),3);
  taskmanager2.run(2);
  EXPECT_EQ(taskmanager2.getNCircles(),5);
}

TEST(TaskManagerTest, InitAndCleanup)
{
  string textstring;
  TaskManager taskmanager;

  DummyUpdater* updater1=new DummyUpdater(textstring);
  DummyUpdater* updater2=new DummyUpdater(textstring);
  DummyAnalyzer* analyzer1=new DummyAnalyzer(textstring);
  DummyAnalyzer* analyzer2=new DummyAnalyzer(textstring);

  taskmanager.addUpdater(updater1);
  taskmanager.addAnalyzer(analyzer1);
  taskmanager.addUpdater(updater2);
  taskmanager.addAnalyzer(analyzer2);

  EXPECT_FALSE(updater1->isInitialized);
  EXPECT_FALSE(updater1->isCleanedUp);
  EXPECT_FALSE(analyzer1->isInitialized);
  EXPECT_FALSE(analyzer1->isCleanedUp);

  taskmanager.initialize();
  EXPECT_TRUE(updater1->isInitialized);
  EXPECT_FALSE(updater1->isCleanedUp);
  EXPECT_TRUE(analyzer1->isInitialized);
  EXPECT_FALSE(analyzer1->isCleanedUp);

  EXPECT_TRUE(updater2->isInitialized);
  EXPECT_FALSE(updater2->isCleanedUp);
  EXPECT_TRUE(analyzer2->isInitialized);
  EXPECT_FALSE(analyzer2->isCleanedUp);

  taskmanager.cleanup();

  EXPECT_TRUE(updater1->isInitialized);
  EXPECT_TRUE(updater1->isCleanedUp);
  EXPECT_TRUE(analyzer1->isInitialized);
  EXPECT_TRUE(analyzer1->isCleanedUp);

  EXPECT_TRUE(updater2->isInitialized);
  EXPECT_TRUE(updater2->isCleanedUp);
  EXPECT_TRUE(analyzer2->isInitialized);
  EXPECT_TRUE(analyzer2->isCleanedUp);
}



