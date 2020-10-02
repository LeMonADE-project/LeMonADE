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

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/Vector3D.h>

int main(int, char**)
{
  /* ***************************************************************************
   * As a first example of how the class Molecules can be used we define a
   * group of particles that have three integer coordinates x,y,z. We add some
   * particles and connect them.
   *
   * The class Molecules is in fact a class template with three template
   * parameters.
   *
   * The first parameter has to be given explicitly and defines the underlying
   * type of a monomer (vertex_type). This means, if we want to use particles with
   * three integer coordinates, the monomer type could be VectorInt3. This is the
   * normal case for the lattice based BFM simulations. However, if we wanted to
   * to use coordinates based on double values, or maybe a type containing more
   * information about the monomer than just three coordinates, we would use the
   * name of that type as first template argument.
   *
   * The second template argument specifies the maximum number of bonds a particle
   * can have. The default value here is 7.
   *
   * The third template parameter specifies a type of information that can be
   * stored on a bond. This defaults to an integer value, but is not used in most
   * applications.
   *
   * The use of the three template parameters, as well as how to add and connect
   * particles, is demonstrated in the following.
   *
   * While you go through this code, you may want to look at the code documentation
   * for the class Molecules to familiarize yourself with the documentation.
   * There you find more detailed information on the different functions.
   * **************************************************************************/

  //we define a Molecules object with integer coordinates, a maximum of 7 allowed
  //bonds per particle, and an integer value stored on a bond

  Molecules<VectorInt3,7,int> molecules1;

  //this is in fact equivalent to the following, because of the default values
  //for the allowed bonds and the integer type for the bond
  Molecules<VectorInt3> molecules2;

  //just to demonstrate the idea, the next line defines a molecules object with
  //only 4 allowed bonds, three coordinates with double values, and an integer stored
  //with every bond
  Molecules<VectorDouble3,4,int> molecules3;

  //now we have three different molecules objects. we can start adding particles
  //to them and connect them. The following lines demonstrate how to add,remove,
  //connect, and disconnect particles.

  //to set the number of particles in the system, Molecules offers the function
  //resize(n).
  //similarly, to obtain the number of particles, use the functon size()
  //here we change the number of particles and print the current size to the
  //standard output

  //at the beginning there are no particles inside
  std::cout<<"**************************************************\n";
  std::cout<<"EXAMPLE 2:\n";
  std::cout<<"number of particles in molecules1: "<<molecules1.size()<<std::endl;

  //now change the number of particles to 1000
  molecules1.resize(1000);
  std::cout<<"number of particles in molecules1: "<<molecules1.size()<<std::endl;

  //now delete the particles again
  molecules1.resize(0);
  std::cout<<"number of particles in molecules1: "<<molecules1.size()<<std::endl;
  std::cout<<"**************************************************\n";

  //to add a single particle at a position x,y,z, one can use the following function
  molecules1.addMonomer(10,12,23);

  //or, one can add a monomer by giving the correct type as an argument to the
  //addMonomer function, as shown below for molecules1 and molecules3:

  //molecules1 contains integer vector coordinates
  VectorInt3 a(1,1,2);
  molecules1.addMonomer(a);

  //molecules3 contains double vector coordinates
  VectorDouble3 b(-123.9,2.2,0.19);
  molecules3.addMonomer(b);

  //now we can add some more particles to molecules1:
  molecules1.addMonomer(2,3,4);
  molecules1.addMonomer(12,13,14);
  molecules1.addMonomer(14,14,17);

  //by now, molecules1 contains 5 particles. the number of particles in a Molecules
  //object can be found using the function size().
  //we print the sizes of all the molecules objects to the standard output:
  std::cout<<"**************************************************\n";
  std::cout<<"size of molecules1:"<<molecules1.size()<<std::endl;
  std::cout<<"size of molecules2:"<<molecules2.size()<<std::endl;
  std::cout<<"size of molecules3:"<<molecules3.size()<<std::endl;
  std::cout<<"**************************************************\n";
  //we can not start to connect particles. we connect particles 1 (at position 1,1,2)
  //and 2 (at position 2,3,4), as well as particles 3 and 4
  molecules1.connect(1,2);
  molecules1.connect(3,4);

  //for demonstrating the principle, we also connect particles 4 and 0 and store some
  //value on the bond. In the definition we chose to allow for integer values to be
  //stored on the bonds, so we store a value of 5 here:
  molecules1.connect(4,0,5);

  //this value on the bond can also be set afterwards by using the function
  //setLinkInfo. here we set the bond information on the bond between monomers
  //1 and 2 to the value of 10;
  molecules1.setLinkInfo(1,2,10);

  //if we want to disconnect two monomers, we can use the function disconnect
  molecules1.disconnect(4,0);

  //the next paragraph shows how to access information about particles. again, we
  //print the information to the standard output:

  //first, the monomer itself can be accessed using the operator[]:
  //This operator returns a reference to the type we specified when setting
  //up molecules1, in this case VectorInt3. Here, we get the position
  //information for the 0th monomer:
  std::cout<<"**************************************************\n";
  std::cout<<"monomer 0 - x coordinate:"<<molecules1[0].getX()<<std::endl;
  std::cout<<"monomer 0 - y coordinate:"<<molecules1[0].getY()<<std::endl;
  std::cout<<"monomer 0 - z coordinate:"<<molecules1[0].getZ()<<std::endl;
  std::cout<<"**************************************************\n";

  //similarly, one can alter the positions using VectorInt3's functions
  //setX(),setY(),setZ():
  //now we change the position of the 0th monomer and print the new position
  //to the screen
  molecules1[0].setX(2);
  molecules1[0].setY(2);
  molecules1[0].setZ(3);
  std::cout<<"**************************************************\n";
  std::cout<<"new position of 0th monomer\n";
  std::cout<<"monomer 0 - x coordinate:"<<molecules1[0].getX()<<std::endl;
  std::cout<<"monomer 0 - y coordinate:"<<molecules1[0].getY()<<std::endl;
  std::cout<<"monomer 0 - z coordinate:"<<molecules1[0].getZ()<<std::endl;
  std::cout<<"**************************************************\n";

  //to obtain bonding information, Molecules offers the following functions:

  //check if two particles are connected, you can use  areConnected(i,j)
  //monomers 0 and 1 are not connected, so the function returns false
  //monomers 1 and 2 are connected, so the function returns true for them
  std::cout<<"**************************************************\n";
  std::cout<<"connection of monomer 0 and 1:"
           <<molecules1.areConnected(0,1)<<std::endl;

  std::cout<<"connection of monomer 0 and 1:"
           <<molecules1.areConnected(2,1)<<std::endl;

  //to get the number of bond partners of a certain particles, use the function
  //getNumLinks(index)
  std::cout<<"monomer number 4 has "<<molecules1.getNumLinks(4)<<" bonds\n";

  //to get the indices of the bond partners of a particle, use the function
  //getNeighborIdx(monomer,bondNumber). here we loop over all bond partners
  //of monomer 4 to find the indices
  for(size_t n=0;n<molecules1.getNumLinks(4);n++){
    std::cout<<"bond partner no "<<n
	     <<" of particle 4 is particle "<<molecules1.getNeighborIdx(4,n)
	     <<std::endl;
  }

  //finally, to obtain the numbers stored on the different bonds, one can use
  //the function getLinkInfo.
  //Here we print the information stored on the bond between monomers 1 and 2,
  //which we previously set to 10:
  std::cout<<"value stored on bond between monomers 1 and 2: "
	   <<molecules1.getLinkInfo(1,2)<<std::endl;


  //to get more information on functions provided by Molecules, you can take
  //a look into the source code documentation.

  return 0;
}
