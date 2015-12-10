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

#ifndef LEMONADE_UTILITY_NUMERICFUNCTIONS_H
#define LEMONADE_UTILITY_NUMERICFUNCTIONS_H

#include <fstream>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/DistanceCalculation.h>


/**
 * @file
 * @brief A collection of functions that may be useful for data analysis
 *
 * @deprecated
 *
 * @todo remove?
 **/


/*************************************************************************
 * global function template squaredRadiusOfGyration
 * ***********************************************************************/
/**
 * @brief Calculation of the radius of gyration Rg2
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo delete?
 **/
template < class MoleculesType > double squaredRadiusOfGyration(const MoleculesType& m)
{
    double sum_sqr = 0.0; VectorDouble3 CoM_sum;
    
    for ( uint i = 0; i < m.size(); ++i)
    {
      CoM_sum.setX( CoM_sum.getX() + m[i].getX() );
      CoM_sum.setY( CoM_sum.getY() + m[i].getY() );
      CoM_sum.setZ( CoM_sum.getZ() + m[i].getZ() );
    }
    
    double inv_N = 1.0 / double ( m.size() );
    
    VectorDouble3 CoM (
	double ( CoM_sum.getX() ) * inv_N, 
	double ( CoM_sum.getY() ) * inv_N,
	double ( CoM_sum.getZ() ) * inv_N);
      
    for ( uint i = 0; i < m.size(); ++i)
    {
      double diff;
      
      diff = double(m[i].getX()) - CoM.getX();
      sum_sqr += diff*diff;
      
      diff = double(m[i].getY()) - CoM.getY();
      sum_sqr += diff*diff;
      
      diff = double(m[i].getZ()) - CoM.getZ();
      sum_sqr += diff*diff;

    }
    
    return sum_sqr / double ( m.size() );
 
}

/*************************************************************************
 * global function template squaredRadiusOfGyrationComponents
 * returns (RgX^2,RgY^2,RgZ^2)  ( RgX^2+RgY^2+RgZ^2=Rg^2)
 * ***********************************************************************/
/**
 * @brief Calculation of the compounds of radius of gyration Rg2
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo delete?
 **/
template < class MoleculesType > VectorDouble3 squaredRadiusOfGyrationComponents(const MoleculesType& m)
{
    VectorDouble3 sum_sqr; VectorDouble3 CoM_sum;
    
    for ( uint i = 0; i < m.size(); ++i)
    {
      CoM_sum.setX( CoM_sum.getX() + m[i].getX() );
      CoM_sum.setY( CoM_sum.getY() + m[i].getY() );
      CoM_sum.setZ( CoM_sum.getZ() + m[i].getZ() );
    }
    
    double inv_N = 1.0 / double ( m.size() );
    
    VectorDouble3 CoM (
	double ( CoM_sum.getX() ) * inv_N, 
	double ( CoM_sum.getY() ) * inv_N,
	double ( CoM_sum.getZ() ) * inv_N);
      
    for ( uint i = 0; i < m.size(); ++i)
    {
      double diffX,diffY,diffZ;
      
      diffX = double(m[i].getX()) - CoM.getX();
      
      
      diffY = double(m[i].getY()) - CoM.getY();
      
      
      diffZ = double(m[i].getZ()) - CoM.getZ();
      
      sum_sqr +=VectorDouble3(diffX*diffX,diffY*diffY,diffZ*diffZ);

    }
    
    return sum_sqr / double ( m.size() );
 
}

/**
 *
 * @brief Calculates the center of mass of all monomers in the argument m.
 *
 * @details Calculates the center of mass of all monomers in the argument m.
 * argument can also be some subgroup of monomers. the function uses the absolute
 * coordinates of the monomers. this means that for example for self-organized clusters 
 * like micelles the result might not be the expected one, because single molecules
 * in the cluster can have diffused into a periodic image of the box before joining
 * the micelle. see also function clusterCenterOfMass(const ClusterType& cluster, const BoxType& box)
 *
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo delete?
 **/
template < class MoleculesType > VectorDouble3 centerOfMass(const MoleculesType& m)
{
    VectorDouble3 CoM_sum;
    
    for ( uint i = 0; i < m.size(); ++i)
    {
      CoM_sum.setX( CoM_sum.getX() + m[i].getX() );
      CoM_sum.setY( CoM_sum.getY() + m[i].getY() );
      CoM_sum.setZ( CoM_sum.getZ() + m[i].getZ() );
    }
    
    double inv_N = 1.0 / double ( m.size() );
    
    VectorDouble3 CoM (
	double ( CoM_sum.getX() ) * inv_N, 
	double ( CoM_sum.getY() ) * inv_N,
	double ( CoM_sum.getZ() ) * inv_N);
      
    return CoM;
 
}

/**
 *
 * @brief calculates the folded center of mass position of a cluster of molecules
 *
 * @details calculates the folded center of mass position of a cluster of molecules.
 * the cluster given as argument should be something like a vector of monomer groups,
 * such that the i-th monomer in the j-th molecule of the cluster can be accessed
 * as cluster[j][i]. IMPORTANT: the function may not give the exact result in some
 * cases, because it is not always clear in a periodic box, if the folded position
 * or the non-folded position of a molecules needs to be used. Therefore, the
 * function first makes an estimate of the center of mass (that can in many cases
 * be wrong by some percent), and then uses this estimate to refine the result.
 *(There may be better ways to implement this...)
 *
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo delete?
 **/
template<class ClusterType, class BoxType>
VectorDouble3 clusterCenterOfMass(const ClusterType& cluster, const BoxType& box)
{
    throw std::runtime_error("function clusterCenterOfMass was disabled because it did not work correctly\n");
//	VectorDouble3 result;
//	VectorDouble3 referencePos;
//	double mass=0.0;
	
//	//first make an estimate of the center of mass of the cluster.
//	//the difficulty is to find out how far a coordinate has to be folded, if
//	//a molecule has diffused into a periodic image of the box before joining
//	//the cluster, when at the same time the size of the cluster is larger
//	//than boxsize/2 in the same dimension.
	
//	//randomly use first molecule of cluster as reference point and calculate
//	//the estimate for the center of mass relative to this molecule. then, at
//	//the end, fold back the estimated COM into the box. the estimate may be
//	//not the correct COM, because if the cluster is larger than boxsize/2, it
//	//is not clear which distance to use between molecules of the cluster
//	//(reduced in periodic space, or true distance)
//	result+=centerOfMass(cluster[0])*double(cluster[0].size());
//	referencePos=centerOfMass(cluster[0]);
	
//	for(size_t n=1;n<cluster.size();++n)
//	{
//		//calculate the reduced distance between the reference point
//		VectorDouble3 currentMolecule=centerOfMass(cluster[n]);
//		VectorDouble3 tmp=Lemonade::calcDistanceVector3D(referencePos,currentMolecule,box);
//		tmp=tmp+referencePos;
//		//add the position of each molecule to the result (estimate)
//		//wehghted by its size
//		result+=double(cluster[n].size())*tmp;
//		mass+=double(cluster[n].size());
		
//		//if(currentMolecule.getZ() != tmp.getZ()) std::cout<<"current "<<currentMolecule<<" folded "<<tmp<<std::endl;
//	}
	
//	//get average
//	result/=mass;
//	//fold back into box
//	result.setX(result.getX()-floor(result.getX()/double(box.getBoxX()))*double(box.getBoxX()));
//	result.setY(result.getY()-floor(result.getY()/double(box.getBoxY()))*double(box.getBoxY()));
//	result.setZ(result.getZ()-floor(result.getZ()/double(box.getBoxZ()))*double(box.getBoxZ()));

//	//now use result as reference position and do as before. now the probability
//	//is relatively small that the distance between a molecule in the cluster
//	//and the reference point is not well defined, because the reference
//	//point is already close to the center of mass. it can still happen that
//	//the result is not exactly correct though.
	
//	//set new reference point
//	referencePos=result;
//	//reset result and mass
//	result=VectorDouble3(0.0,0.0,0.0);
//	mass=0.0;
	
//	//calculate center of mass
//	for(size_t n=0;n<cluster.size();++n)
//	{
//		VectorDouble3 currentMolecule=centerOfMass(cluster[n]);
//		VectorDouble3 tmp=Lemonade::calcDistanceVector3D(referencePos,currentMolecule,box);
//		tmp=tmp+referencePos;
		
//		result+=double(cluster[n].size())*tmp;
//		mass+=double(cluster[n].size());
		
//		//if(currentMolecule.getZ() != tmp.getZ()) std::cout<<"current "<<currentMolecule<<" folded "<<tmp<<std::endl;
//	}
	
//	result/=mass;
			
//	result.setX(result.getX()-floor(result.getX()/double(box.getBoxX()))*double(box.getBoxX()));
//	result.setY(result.getY()-floor(result.getY()/double(box.getBoxY()))*double(box.getBoxY()));
//	result.setZ(result.getZ()-floor(result.getZ()/double(box.getBoxZ()))*double(box.getBoxZ()));

//	return result;
		
		
};



#endif /* LEMONADE_UTILITY_NUMERICFUNCTIONS_H */
