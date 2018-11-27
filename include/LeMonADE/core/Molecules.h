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

#ifndef LEMONADE_CORE_MOLECULES_H
#define LEMONADE_CORE_MOLECULES_H

#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>
#include <stdexcept>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/core/ConnectedDecorator.h>
#include <LeMonADE/core/MoleculesRead.h>
#include <LeMonADE/core/MoleculesWrite.h>
#include <LeMonADE/io/FileImport.h>

/**
 * @file
 *
 * @class Molecules
 *
 * @brief \b Basic \b structure \b holding \b all \b monomers \b in \b the \b system \b organized \b in \b a \b graph.
 *
 * This is the basic structure holding all monomers of the system.
 * The structure is organized as a graph consisting of vertices and edges.
 * Typically, the single monomers are the vertices, while the edges provide
 * the connectivity information (However, since this is a class template with the
 * vertex- and edge-types as template parameters, on could in principle put
 * any type of information on vertices and edges.)
 * Example: Molecules<5,Position> is a graph with "Position" objects as monomers,
 * each of which can have a maximum of five bonds.
 *
 * @tparam max_connectivity maximum allowed connectivity of the sigle monomers
 * @tparam Vertex vertex type of the graph, typically associated with single
 * monomer type, for instance "Position".
 * @tparam Edge edge type of the graph, by default of type int. Any other type
 * could be used to store additional information on the edges
 *
 **/


/*****************************************************************************/
//class definition (definitions of members see below)
/*****************************************************************************/

//template < class V, uint MC, class E > class Molecules;

template < class Vertex, uint max_connectivity = 7, class Edge = int > class Molecules
{


  //! Adds a neighbor list to the vertices
  typedef Connected < Vertex, max_connectivity > internal_vertex_type;

  public:

  //! Definition of the vertex type
  typedef Vertex vertex_type;

  //! Definition of the edge type
  typedef Edge edge_type;

  /*****************************************************************************/
  //constructors
  /*****************************************************************************/
  //! Standard constructor that initialize with zero vertices at age equal 0
  Molecules():vertices(0),myAge(0){}

  //! Conversion constructor
  template < class V, uint m, class E> Molecules (const Molecules<V,m,E>& src);

  Molecules<Vertex,max_connectivity, Edge>& operator=  (const Molecules<Vertex,max_connectivity, Edge>& src);

  Molecules<Vertex,max_connectivity, Edge>& operator+=  (const Molecules<Vertex,max_connectivity, Edge>& src);

  /*****************************************************************************/
  //access to general properties of the graph
  /*****************************************************************************/
  /**
   * @brief Returns the maximum connectivity of each vertex (not the actual one).
   *
   * @return The maximum number of edges each vertex (not the actual one).
   */
  static uint 	getMaxConnectivity() 	{return max_connectivity;}

  /**
   * @brief Resizes the graph to a new size creating/destroying vertices.
   *
   * Resizes the graph to \a newGraphSize number of elements. Initially, all newly created vertex
   * positions are set to (0,0,0) as default. If \a newGraphSize is smaller than the actual size you
   * lose vertex (monomer) information. Note, that the new graph size is set to \a newGraphSize and
   * all indexing (e.g. addMonomer() ) refers to the new size.
   *
   * @param newGraphSize the number of vertices the graph should hold
   */
  void     	resize(uint32_t newGraphSize) {vertices.resize(newGraphSize);}

  /**
   * @brief Returns the actual size of the graph, esp. the number of vertices.
   *
   * @return The actual size of the graph, esp. the number of vertices.
   */
  uint32_t   	size() const {return vertices.size();}

  /**
   * @brief Returns the actual time in the graph.
   *
   * @return The actual time in the graph.
   */
  const uint64_t 	getAge() const 			{return myAge;}


  /**
   * @brief Sets the actual time to a given value \a newAge in the graph.
   *
   * @param newAge The time to set in the graph.
   */
  void			setAge(const uint64_t newAge)      {myAge = newAge ;}

  /*****************************************************************************/
  //access monomers
  /*****************************************************************************/
  /**
   * @brief Adds a vertex (monomer) of \a vertex_type to the graph.
   *
   * This function increases the size of the graph by one and
   * maybe increases the capacity of the graph to add a monomer, too.
   * It returns the index of the new vertex (monomer) in the graph.
   * This should be useful for direct access of vertices (molecules) or
   * for the establishing the connectivity.
   *
   * @param [in] monomer A \a vertex_type monomer to add to the graph.
   * @return The index of the new vertex (monomer) in the graph.
   */
  uint32_t addMonomer(vertex_type monomer)
  {
	  vertices.push_back(monomer);
	  return (vertices.size()-1);
  }

  /**
   * @brief Adds a vertex (monomer) with position (\a x, \a y, \a z) to the graph.
   *
   * A vertex (monomer) with \a vertex_type with position (\a x, \a y, \a z) is created
   * and added the graph.
   * This function increases the size of the graph by one and
   * maybe increases the capacity of the graph to add a monomer, too.
   * It returns the index of the new vertex (monomer) in the graph.
   * This should be useful for direct access of vertices (molecules) or
   * for the establishing the connectivity.
   *
   * @param[in] x x-coordinate in the Cartesian space
   * @param[in] y y-coordinate in the Cartesian space
   * @param[in] z z-coordinate in the Cartesian space
   *
   * @tparam CoordinateType Type of the coordinate specified in vertex_type/Vector3D
   *
   * @return The index of the new vertex (monomer) in the graph.
   */
  template < class CoordinateType >
  uint32_t addMonomer(const CoordinateType& x, const CoordinateType& y, const CoordinateType& z)
  {
	  vertex_type monomer;
	  monomer.setX(x);
	  monomer.setY(y);
	  monomer.setZ(z);
	  vertices.push_back(monomer);

	  return (vertices.size()-1);
  }

  //! Connect two vertices (monomers) with indices a and b in the graph.
  void connect(uint32_t a, uint32_t b, const Edge& edgeVal = Edge() );

  //! Disconnect two vertices (monomers) with indices a and b in the graph.
  void disconnect(uint32_t a, uint32_t b);


  /**
   * @brief Returns the number of connections of vertex (monomer) with index idx.
   *
   * @throw <runtime_error> if the vertex (monomer) with index \a idx doesn´t exist in the graph.
   *
   * @param idx The index of vertex (monomer) in the graph.
   * @return The number of connections of vertex (monomer) with index \a idx.
   */
  uint32_t   getNumLinks(uint32_t idx) const 	{return vertices.at(idx).getNumLinks();}


  //! Returns the index in the graph of the j-th connection (bond partner) of the vertex with index \a idx.
  uint32_t   getNeighborIdx(uint32_t idx, uint32_t j) const;


  /**
   * @brief Provides direct access to the vertex (monomer) with index \a idx in the graph.
   *
   * @param idx The index of vertex (monomer) in the graph.
   * @return The vertex (monomer) in the graph with index \a idx.
   *
   * @todo change return type to internal_vertex_type?
   *
   * @todo use vertices[idx] instead of vertices.at(idx) to avoid boundary check.
   *
   * @todo catch exception and use DEBUG-mode?
   */
  const Vertex& operator[](uint32_t idx) const {return vertices.at(idx);}

  /**
   * @brief Provides direct access to the vertex (monomer) with index \a idx in the graph.
   *
   * @param idx The index of vertex (monomer) in the graph.
   * @return The vertex (monomer) in the graph with index \a idx.
   *
   * @todo change return type to internal_vertex_type?
   *
   * @todo use vertices[idx] instead of vertices.at(idx) to avoid boundary check.
   *
   * @todo catch exception and use DEBUG-mode?
   */
  Vertex& operator[](uint32_t idx)	{return vertices.at(idx);}

  //! Returns the information \a Edge stored on the connection (bond) between vertices with indices a and b
  const Edge& getLinkInfo(uint32_t a, uint32_t b) const;

  //! Set the information on the \a Edge on the connection (bond) between vertices with indices a and b
  void setLinkInfo(uint32_t a, uint32_t b, Edge edge);

  //! Returns the total number of bonds/edges in graph
  uint32_t getTotalNumLinks() const;

  //! Returns if vertices (monomers) with index a and b are connected.
  bool areConnected(uint32_t a, uint32_t b) const;


  //! Clear the graph destroying all vertices&edges and setting the age to zero.
  void clear(){resize(0); edges.clear(); myAge=0;}

  /** Delete all the edges (bonds) in the graph. This does not destroy the vertices&edges.
   * @todo check where this might be used?
   */
  void clearBonds();

  //! returns the edges of the graph 
  std::map < std::pair < uint32_t, uint32_t > , Edge > getEdges() const {return edges;}

private:

  /**
   * @brief Pair of connected vertices (monomers), esp. their indices in the graph.
   *
   * @typedef IndexPair
   */
  typedef std::pair < uint32_t, uint32_t > IndexPair;

  //! Stores the vertices (monomers) of the graph,
  std::vector < internal_vertex_type > vertices;


  //! Stores the connection information
  std::map < IndexPair, Edge > edges;

  //! Age of the configuration in Monte-Carlo steps (MCS)
  uint64_t myAge;
};


/*****************************************************************************/
//members of class Molecules
/*****************************************************************************/

/**
 * Copys and Converts one type of Molecules object into another, i.e. copies vertex,
 * connectivity and age information, and does type-conversion
 *
 * @param src A reference to another graph to copy from.
 * @tparam m Maximum allowed connectivity of the source graph
 * @tparam V Vertex type of the source graph
 * @tparam E Edge type of the source graph
 */
template < class Vertex, uint max_connectivity, class Edge>
template < class V, uint m, class E >
Molecules<Vertex,max_connectivity, Edge>::Molecules (const Molecules<V,m,E>& src)
{

	clear();
	vertices.resize(src.size());

	for ( uint i = 0 ; i < vertices.size(); ++i) {

		//this initializes a new vertex, all connectivity information is lost
		this->vertices[i]=src[i];
		//copy the connectivity and edges, too
		uint32_t nNeighbours=src.getNumLinks(i);
		try{
			for(uint j=0; j<nNeighbours;++j){
				connect(i, src.getNeighborIdx(i,j), src.getLinkInfo(i,src.getNeighborIdx(i,j)));
			}
		}
		catch(std::runtime_error& e){
			std::stringstream errormessage;
			errormessage<<"Error when copying connectivity information in "
			<<"Molecules copy constructor.\n Original error:\n"
			<<e.what();

			throw std::runtime_error(errormessage.str());
		}

	}

	this->myAge = src.getAge();

}

/**
 * Assign Molecules object to another, i.e. copies vertex,
 * connectivity and age information, and does type-conversion
 *
 * @param src A reference to another graph to copy from.
 * @tparam m Maximum allowed connectivity of the source graph
 * @tparam V Vertex type of the source graph
 * @tparam E Edge type of the source graph
 */
template < class Vertex, uint max_connectivity, class Edge>
Molecules<Vertex,max_connectivity, Edge>& Molecules<Vertex, max_connectivity, Edge>::operator=  (const Molecules<Vertex,max_connectivity, Edge>& src)
{

	clear();
	vertices.resize(src.size());

	for ( uint i = 0 ; i < vertices.size(); ++i) {
		//copy all vertices.
		//this does not copy the connectivity info, but it is copied
		//in the loop below
		this->vertices[i]=src[i];

		//copy the connectivity and edges, too
		uint32_t nNeighbours=src.getNumLinks(i);

		for(uint j=0; j<nNeighbours;++j){
			// copy connections and edges contents
			connect(i, src.getNeighborIdx(i,j), src.getLinkInfo(i,src.getNeighborIdx(i,j)));
		}

	}

	this->myAge = src.getAge();
	return *this;
}

/**
 * Adding Molecules object to another, i.e. copies vertex,
 * connectivity and age information, and does type-conversion
 *
 * @todo testing!!!
 *
 * @param src A reference to another graph to copy from.
 * @tparam m Maximum allowed connectivity of the source graph
 * @tparam V Vertex type of the source graph
 * @tparam E Edge type of the source graph
 */
template < class Vertex, uint max_connectivity, class Edge>
Molecules<Vertex,max_connectivity, Edge>& Molecules<Vertex, max_connectivity, Edge>::operator+=  (const Molecules<Vertex,max_connectivity, Edge>& src)
{
	uint32_t oldsize =  vertices.size();

#ifdef DEBUG
	std::cout << "old size: " << oldsize<< std::endl;
	std::cout << "src size: " << src.size()<< std::endl;
	std::cout << "new size: " << (oldsize+src.size())<< std::endl;
#endif //DEBUG

	vertices.resize(oldsize+src.size());

#ifdef DEBUG
	std::cout << "new size after resize: " << (vertices.size())<< std::endl;
#endif //DEBUG

	for ( uint i = 0 ; i < src.size(); ++i) {
		//copy all vertices.
		this->vertices[oldsize+i]=static_cast<const typename Molecules<Vertex,max_connectivity, Edge>::vertex_type&>(src[i]);

		//copy the edges, too
		uint32_t nNeighbours=src.getNumLinks(i);//vertices[i].getNumLinks();

		for(uint j=0; j<nNeighbours;++j){

			connect(oldsize+i, oldsize+src.getNeighborIdx(i,j));

		}

	}

	this->myAge = src.getAge();
	return *this;
}
/**
 * Optionally a value \a edgeVal can be assigned to the edge - default is value-initialized (most Zero).
 *
 * @throw <std::range_error> If one index exceed the boundary this displays detailed error messages.
 *
 * @param a The index \a a of vertex (monomer) in the graph.
 * @param b The index \a b of vertex (monomer) in the graph.
 * @param edgeVal Optional value for the edge - default is value-initialized (most Zero).
 */
template < class Vertex, uint max_connectivity, class Edge>
void Molecules <Vertex,max_connectivity,Edge>::connect(uint32_t a, uint32_t b, const Edge& edgeVal)
{
    //add connection information to both vertices
    //display detailed errormessage, if errors occur

    // check if it is me
    if (a==b) {
    	std::cout << "trying to connect " << a << " to " << b << std::endl;
      return;
    }

    try{
      vertices.at(a).connect(b);
      vertices.at(b).connect(a);
    }
    //method vertices.at(int) might throw out_of_range exception.
    catch(std::out_of_range& exception){
    	std::stringstream messagestream;
      messagestream<<"Molecules::connect(int a, int b): a="<<a<<" and b="<<b<<std::endl;

      if(a>=size()){
	messagestream<<"a is out of range"<<std::endl;
      }
      if(b>=size()){
	messagestream<<"b is out of range"<<std::endl;
	//undo changes in a's connectivity
	vertices.at(a).disconnect(b);
      }
      //re-throw exception
      throw std::range_error(messagestream.str());
    }
    //catch other errormessages here. we normally use runtime_error
    catch(std::runtime_error& e){
      std::stringstream messagestream;
      messagestream<<"Molecules::connect(int a, int b): a="<<a<<" and b="<<b<<std::endl;
      messagestream<<"Indices were: a="<<a<<" ,b="<<b<<std::endl;
      messagestream<<e.what();
      //re-throw exception
      throw std::runtime_error(messagestream.str());
    }

    //store connection information in edges map, if connecting went fine
    IndexPair edge_key(std::min(a,b),std::max(a,b));
    edges[edge_key] = edgeVal;
}


/**
 * Also the value of the edge between the vertices (monomers) with index a and b will be erased.
 *
 * @param a The index \a a of vertex (monomer) in the graph.
 * @param b The index \a b of vertex (monomer) in the graph.
 */
template < class Vertex,uint max_connectivity,  class Edge>
void Molecules <Vertex,max_connectivity, Edge>::disconnect(uint32_t a, uint32_t b)
{
  //erase connection from both vertices and the edges map
   vertices.at(a).disconnect(b);
   vertices.at(b).disconnect(a);
   IndexPair idx(std::min(a,b),std::max(a,b));
   edges.erase(idx);
}



  /**
   * @param idx The index of vertex (monomer) in the graph.
   * @param j The j-th connection of the vertex (monomer) in the graph.
   *
   * @throw <runtime_error> if no \a j-th connection exist in the graph.
   *
   * @return Returns the index of the connected vertex (monomer) in the graph to vertex with index \a idx.
   */
template < class Vertex, uint max_connectivity, class Edge>
uint32_t Molecules <Vertex,max_connectivity,Edge>::getNeighborIdx(uint32_t idx, uint32_t j) const
{
//in the debug mode we give more useful output in case of error,
//than the usual exception coming from the return statement can provide
#ifdef DEBUG
	if(j>=getNumLinks(idx)){
		std::stringstream messagestream;
		messagestream<<"getNeighborIdx(i, j): j="<<j<<" is out of range for monomer i="<<idx<<std::endl;
		throw std::runtime_error(messagestream.str());
	}
#endif
	return vertices.at(idx).getNeighborIdx(j);
}


/**
 * @param a The index \a a of vertex (monomer) in the graph.
 * @param b The index \a b of vertex (monomer) in the graph.
 *
 * @throw <std::runtime_error> If no link (connection, bond) between the vertices exist.
 *
 * @return Value of the edge between the vertices.
 */
template < class Vertex, uint max_connectivity, class Edge>
const Edge& Molecules <Vertex,max_connectivity,Edge>::getLinkInfo(uint32_t a, uint32_t b) const
{
  IndexPair idx(std::min(a,b),std::max(a,b));
  //try to find the information in edges map. if bond does not exist, display
  //detailed error and throw exception
  try{
    return edges.at(idx);
  }
  catch(std::out_of_range& exception){
    std::stringstream errormessage;
    errormessage <<"Molecules::getLinkInfo(uint a, uint b): with a="<<a<<" and b="<<b
      <<". Bond does not exist."<<std::endl;
    throw std::runtime_error(errormessage.str());
  }
}


/**
 * @param a The index \a a of vertex (monomer) in the graph.
 * @param b The index \a b of vertex (monomer) in the graph.
 * @param edge Value of the edge between the vertices.
 *
 * @throw <std::runtime_error> If no link (connection, bond) between the vertices exist.
 */
template < class Vertex, uint max_connectivity, class Edge>
void Molecules <Vertex,max_connectivity,Edge>::setLinkInfo(uint32_t a, uint32_t b, Edge edge)
{
  IndexPair idx(std::min(a,b),std::max(a,b));
  //try to set the information in edges map. if bond does not exist, display
  //detailed error and throw exception
  try{
	  edges.at(idx)=edge;
  }
  catch(std::out_of_range& exception){
    std::stringstream errormessage;
    errormessage <<"Molecules::setLinkInfo(uint a, uint b, Edge edge): with a="<<a<<" and b="<<b
      <<". Bond does not exist."<<std::endl;
    throw std::runtime_error(errormessage.str());
  }
}


/**
 * @return The overall number of all edges (connections/bonds) in the graph.
 */
template < class Vertex, uint max_connectivity, class Edge>
uint32_t Molecules <Vertex,max_connectivity,Edge>::getTotalNumLinks() const
{
  return edges.size();
}


/*****************************************************************************/
/**
 * @fn bool Molecules::areConnected(int a, int b) const
 * @brief checks if monomers a and b are connected
 * @details
 * */
/*****************************************************************************/
/**
 * This function looks up the edge (connection/bond) in the edges map between the
 * vertices with index a and b in the graph. It only checks if the edge (connection/bond) is existing,
 * not if the indices are valid.
 *
 * @param a The index \a a of vertex (monomer) in the graph.
 * @param b The index \a b of vertex (monomer) in the graph.
 *
 * @return True if the edge exist and False if not.
 */
template<class Vertex, uint max_connectivity, class Edge>
bool Molecules<Vertex, max_connectivity, Edge>::areConnected(uint32_t a, uint32_t b) const
{
	//check for boundaries
	if ((std::min(a, b) < 0) || (std::max(a, b) > size())) {
		return false;
	}

	IndexPair idxPair(std::min(a, b), std::max(a, b));
	if (edges.find(idxPair) == edges.end())
		return false;
	else
		return true;
}


/**
 * This function loops over all vertices in the graph, checks if they are connected,
 * and disconnect them if applicable.
 */
template<class Vertex, uint max_connectivity, class Edge>
void Molecules<Vertex, max_connectivity, Edge>::clearBonds()
{
	for(uint32_t n=0;n<vertices.size();n++)
	{
		for(uint32_t m=0;m<n;m++)
		{
			if(areConnected(m,n)) disconnect(m,n);
		}
	}
}



#endif
