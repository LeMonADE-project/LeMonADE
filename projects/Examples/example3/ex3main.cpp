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

/* ***************************************************************************
 * As a first example of how the classes Ingredients and Molecules can be used
 * together, we define a group of particles that have three integer coordinates.
 * This is just the same as in example 2. Here, we additionally place the
 * particles into an Ingredients environment defined by features. This example
 * program shows how to access the Molecules object contained within Ingredients,
 * and how to use features to change the environment.
 * **************************************************************************/

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/utility/Vector3D.h>


int main(int, char**)
{
  /* ***************************************************************************
   * The definition of the Ingredients and the features used can conveniently done
   * using some typedefs at the beginning of the main function. These usually make
   * the first part of the main() function of a LeMonADe program.
   * **************************************************************************/
  /* ***************************************************************************
   * In the following four lines we will define an Ingredients object that contains
   * our molecules in a box. The details of these lines are explained in the comment
   * block below.
   * ***************************************************************************/

  typedef LOKI_TYPELIST_1(FeatureBox) Features1;
  typedef ConfigureSystem<VectorInt3,Features1> Config1;
  typedef Ingredients<Config1> SimulationSystem1;

  SimulationSystem1 system1;

  /* ***************************************************************************
   * The first three lines of the code above are only typedefs, that is they only
   * define alias shortcut names to make the code more readable. In fact, all four
   * lines of code above could also be written as
   * Ingredients<ConfigureSystem<VectorInt3,LOKI_TYPELIST_1(FeatureBox)> > system;
   * Here, we use the separate typedefs to make the different parts easier to
   * explain.
   *
   * The list of features to be used must be specified in form of a Loki-typelist.
   * This is what happens in the first of the three lines above. If more than
   * just one feature is used, one uses a typelist of the respective length. For
   * example, for for features named Feature1...4, one would define the list of
   * features as LOKI_TYPELIST_4(Feature1,Feature2,Feature3,Feature4). An example
   * of this follows at the bottom of this example program.
   * For detailed explanations of what a typelist does internally, you can have
   * a look at the documentation of the loki library at
   * http://loki-lib.sourceforge.net,
   * or other sources explaining the loki library. For using LeMonADe, this is not
   * necessary to know in detail.
   *
   * In the second and third lines, the data type of Ingredients is defined.
   * Ingredients uses a Configuration type as template parameter (third line)
   * from which it gets all information about the used features.
   * This type is defined in the second line and specifies internally the Molecules
   * object to be used, and a data type containing all information from the features.
   * Again, the details are not too important for the beginning. It is mostly
   * important to know that the type ConfigureSystem uses two template arguments,
   * the first one being the vector type for Molecules as in example 2, and the
   * second one being the list of features.
   *
   * In the third line, finally, the actual Ingredients object is created.
   *
   * One of the later tutorial programs deals more with the details of features.
   * If the syntax above seems confusing at the beginning, don't worry. It will
   * be explained later. For now we only show how to use features. For this
   * purpose, the only thing worth mentioning is that a feature in this context
   * is in principle just another C++ class. What the ConfigureSystem class does
   * is creating a new class that contains all functions and data of all the
   * features in the list. This way, if FeatureBox has a function called
   * setBoxX(boxsize), also Ingredients will have this function.
   * In the next part of this program this principle is demonstrated, and
   * and particles are added to the system like in the previous example.
   * ***************************************************************************/

  /* ***************************************************************************
   * To create particles, one can get access to the Molecules object inside
   * Ingredients by using the two methods
   * - modifyMolecules()  (for editing the particles)
   * - getMolecules()     (for read only access)
   * These two functions return the Molecules object, which can then be used as
   * explained in example 2. The next lines show how to access the particles.
   * ***************************************************************************/

  system1.modifyMolecules().resize(10);   //resize as before
  system1.modifyMolecules()[0].setX(123); //set x coordinate of 0th particle
  system1.modifyMolecules()[9].setZ(2); //set z coordinate of 9th particle

  std::cout<<"**************************************************\n";
  std::cout<<"EXAMPLE 3: Ingredients\n"
	   <<"The 0th particle is at position "
	   <<system1.getMolecules()[0].getX()<<" "
	   <<system1.getMolecules()[0].getY()<<" "
	   <<system1.getMolecules()[0].getZ()<<std::endl;

  std::cout<<"The 9th particle is at position "
	   <<system1.getMolecules()[9].getX()<<" "
	   <<system1.getMolecules()[9].getY()<<" "
	   <<system1.getMolecules()[9].getZ()<<std::endl;

  /* *************************************************************************
   * Since we have used the FeatureBox, we also want to be able to modify the
   * properties of the box. This can simply be done directly from the
   * Ingredients object. FeatureBox, just like any other feature, is a class,
   * and the ConfigureSystem used above defines Ingredients in such a way that
   * all functions provided by the features, are now contained in Ingredients.
   * In practice, FeatureBox provides functions for setting the size and
   * periodicity of the box. So now these can be set in Ingredients, as shown
   * below. The functions setBoxX() etc. are defined in the source file of
   * FeatureBox, namely in src/feature/FeatureBox.h .
   * *************************************************************************/

  system1.setBoxX(64);         //define the box size 64 in x direction
  system1.setBoxY(32);         //define the box size 32 in y directionw
  system1.setBoxZ(128);        //define the box size 128 in z directionw

  system1.setPeriodicX(true);  //make the box periodic in x,y direction
  system1.setPeriodicY(true);  //and non-periodic in z-direction
  system1.setPeriodicZ(false);

  /* *************************************************************************
   * Ingredients provides a function to make sure that the information from
   * all features leads to a consistent simulation system. This function is
   * called synchronize().
   * In general, it is a good idea, and often even necessary, to call the
   * synchronize function after making changes in ingredients.
   * Consider for example the case above, where the box is not periodic in
   * z direction. If we now set the coordinate of some particle to be larger
   * than 128 (the size of the box in z-direction), or negative, the system
   * becomes inconsistent. Calling synchronize would in that case give an
   * appropriate error message. Below it is shown how to use synchronize.
   * To see the effect, you can change the coordinate of some particle to
   * an invalid value, and then compile and run the example program.
   * ************************************************************************/

  system1.synchronize(system1);

  /* *************************************************************************
   * Using the same principle as above, we can easily define a more complex
   * environment, which is defined by a box, but additionally a set of allowed
   * bond vectors (for the BFM model). The handling of this set of bond vectors
   * is taken care of by FeatureBondset. Setting up such an environment works
   * exactly the same way as before. Note that we now use LOKI_TYPELIST_2,
   * because we use two features to define the system.
   * in addition the the settings for the box, system2 also offers functions for
   * defining the set of allowed bond vectors:
   * *************************************************************************/

  typedef LOKI_TYPELIST_2(FeatureBox,FeatureBondset<>) Features2;
  typedef ConfigureSystem<VectorInt3,Features2> Config2;
  typedef Ingredients<Config2> SimulationSystem2;
  SimulationSystem2 system2;

  system2.modifyMolecules().resize(100);

  system2.setBoxX(64);
  system2.setBoxY(64);
  system2.setBoxZ(128);

  system2.setPeriodicX(true);
  system2.setPeriodicY(true);
  system2.setPeriodicZ(true);

  //add two bond vectors. the first three arguments are the vector
  //the fourth is a unique bond ID
  system2.modifyBondset().addBond(1,2,2,100);
  system2.modifyBondset().addBond(2,1,2,101);

  system2.synchronize(system2);

  /*************************************************************************
   * After this example it has hopefully become clear, how LeMonADe can be
   * customized for more special and/or complicated simulation setups. To
   * add functionality to Ingredients, like an external potential, or other
   * things, the user can write a Feature containing this functionality,
   * and simply add the feature to the typelist. There is no need to change
   * existing features.
   *
   * Example 6 explains how the features are automatically taken into account
   * in the simulation process and how to write your own features.
   *
   ***************************************************************/

  return 0;
}
