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
 * This example program demonstrates the use of the class template
 * Vector3D, providing a simple 3 dimensional vector class.
 *
 * Vector3D is defined in src/utility/Vector3D.h
 * *********************************************************************/

#include <iostream>

#include <LeMonADE/utility/Vector3D.h>

int main(int, char**)
{
  /* ********************************************************************
   * Vector3D is a class template with the template parameter specifying
   * the underlying basic data type. For example:
   *
   * Vector3D<int32_t> vector32; //creates a 3d vector with 32bit integer numbers
   * Vector3D<int64_t> vector64; //creates a 3d vector with 64bit integer numbers
   * Vector3D<double> vectorDouble; //creates a 3d vector with double numbers
   *
   * there are typedefs for the most used template parameters:
   * VectorInt3 is equivalent to Vector<int32_t>
   * VectorUint3 is equivalent to Vector<uint32_t>
   * VectorLong3 is equivalent to Vector<int64_t>
   * VectorUlong3 is equivalent to Vector<uint64_t>
   * VectorInt3 is equivalent to Vector<float>
   * VectorInt3 is equivalent to Vector<double>
   * and some more (look into the source code documentation for more details)
   * *********************************************************************/

  //create some vectors for later use in the example
  VectorInt3 vector32;
  VectorInt3 vector32_a;
  VectorInt3 vector32_b;

  VectorDouble3 vectorDouble;
  //one can also give the position arguments to the constructor
  VectorInt3 vector32_c(1,0,0);

  /* ********************************************************************
   * the vector coordinates can be accessed and modified with
   * getX(),getY(),getZ()
   * setX(value),...
  * *********************************************************************/
  //change the coordinates of the integer vector
  vector32.setX(10);
  vector32.setY(-12);
  vector32.setZ(23);
  //change the coordinates of the double vector
  vectorDouble.setX(0.123);
  vectorDouble.setY(-12.356);
  vectorDouble.setZ(3.4);

  //print the coordinates to the standard output
  std::cout<<"EXAMPLE 0: Vector3D\n";
  std::cout<<"vector32 coordinates: "
	   <<vector32.getX()<<" "
	   <<vector32.getY()<<" "
	   <<vector32.getZ()<<std::endl;
  std::cout<<"vectorDouble coordinates: "
	   <<vectorDouble.getX()<<" "
	   <<vectorDouble.getY()<<" "
	   <<vectorDouble.getZ()<<std::endl;


  /* ********************************************************************
   * Vector3D provides basic mathematical, assignment and logical operations:
   * operators +,-,* (scalar product), * (vector times scalar), /,
   * +=, *=, /=
   * operator ==
   * there exists also a function to calculate the vector product
   * the following lines shows some examples:
   **********************************************************************/
  vector32_a=vector32;                          //assignment
  int32_t scalarProduct=vector32*vector32_a;    //scalar product
  vector32_b=-10*vector32;                      //product with scalar value

  VectorInt3 vectorProduct=crossProduct(vector32,VectorInt3(1,0,0));

  //print results to screen standard output
  std::cout<<"scalar product of vector32 and vector32_a: "
	   <<scalarProduct<<std::endl;
  std::cout<<"vector product of vector32 and (1,0 0):"
	   <<vectorProduct.getX()<<" "
	   <<vectorProduct.getY()<<" "
	   <<vectorProduct.getZ()
	   <<std::endl;

  /* ********************************************************************
   * One important feature of Vector3D is that it does not allow implicit
   * conversion with loss of data. That means you can assign the values of
   * an integer vector to a double vector, but not the other way around.
   * Trying the latter will result in a compiler error
   **********************************************************************/

  VectorDouble3 vectorDouble_a=vector32;  //this works
  std::cout<<"Assignment of VectorInt3 to VectorDouble3 works:\n"
		<<"\tDouble\tInt\n"
	   <<"x\t"<<vectorDouble_a.getX()<<"\t"<<vector32.getX()<<"\n"
	   <<"y\t"<<vectorDouble_a.getY()<<"\t"<<vector32.getY()<<"\n"
	   <<"z\t"<<vectorDouble_a.getZ()<<"\t"<<vector32.getZ()<<"\n"
	   <<std::endl;

  //uncommenting the following line will result in a compiler error
  //VectorInt3 vector32_d=vectorDouble;

  /* ********************************************************************
   * The class Vector3D provides some more useful functions for getting the
   * length of the vector, normalizing the vector, etc. To find out more,
   * look into the source code documentation.
   ***********************************************************************/

  return 0;
}
