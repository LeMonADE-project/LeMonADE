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

#ifndef LEMONADE_CORE_CONNECTEDDECORATOR_H
#define LEMONADE_CORE_CONNECTEDDECORATOR_H

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>


/**
 * @file
 * @class Connected
 * @brief Adds a neighbor list to an arbitrary class Vertex.
 *
 * @tparam Vertex
 * @brief The class to be decorated
 *
 * @tparam max_connectivity
 * @brief Maximum number of connected objects. Default is \a 7.
 * @todo implement copy operator (or using the copy operator from the base class. But this
 * means one only copies the positions and not the connections.
 **/
template <class Vertex, uint max_connectivity = 7>
class Connected: public Vertex
{

public:

	/**********************************************************************/
	//constructors
	/**********************************************************************/
	/**
	 * @brief Construct an empty neighbor list
	 */
	Connected() :	Vertex(), counter(0) {}

	/**
	 * @brief Construct an empty neighbor list, but copies the Vertex internally.
	 *
	 * @param vertex Reference to Vertex to copy from.
	 */
	Connected(const Vertex& vertex) :Vertex(vertex), counter(0) {}


	//! Conversion constructor (from different Connected type)
	template<class V, uint MC>
	Connected(const Connected<V, MC>& src);

	/**********************************************************************/
	//connectivity
	/**********************************************************************/
	/**
	 * @brief Get the established number of connections/links/bonds of the Vertex.
	 *
	 * @return The number of connections of Vertex.
	 */
	uint32_t getNumLinks() const {
		return counter;
	}

	/**
	* @brief Returns the maximum connectivity of each vertex (not the actual one).
	*
	* @return The maximum number of edges each vertex (not the actual one).
	*/
	static uint getMaxConnectivity() {
		return max_connectivity;
	}

	//! Returns the index of the i-th connection (bond partner) of the vertex.
	uint32_t getNeighborIdx(uint32_t i) const;

	//! Connect this Vertex (monomer) to another Vertex with index \a b.
	void connect(int32_t b);

	//! Disconnect this Vertex (monomer) from another Vertex with index \a b.
	void disconnect(uint32_t b);

	/**********************************************************************/
	//operators
	/**********************************************************************/
	//! Assignment operator doing a type-conversion from a different Connected type
	template<class V, uint MC>
	void operator =(const Connected<V, MC> &val);

private:

	//! Checks, if there is possibility for more connections, throws exception otherwise.
	void checkMaxConnectivity();

	//! Array contains the indices of connected vertices
	uint32_t links[max_connectivity];

	//! Current number of connections
	uint32_t counter;

};

/****************************************************************************/
//implementation of member functions
/****************************************************************************/


/**
 * @details Conversion and copy of the specified neighbor list \a src to this implementation
 *
 * @throw <std::runtime_error> if copying not possible due to limited max_connectivity
 *
 * @param[in] src Neighbor list to convert.
 */
template<class Vertex,uint max_connectivity>
template<class V, uint MC>
Connected<Vertex,max_connectivity>::Connected(const Connected<V, MC>& src)
	:Vertex(src) {

	//if max_connectivity not high enough, throw exception
	if (this->getMaxConnectivity() < src.getNumLinks()) {
		throw std::runtime_error("class Connected - constructor: Tried to import from connected Vertex with larger connectivity.");
	}
	//else, copy content
	else {
		for (uint i = 0; i < src.getNumLinks(); ++i) {
			links[i] = src.getNeighborIdx(i);
		}
		counter = src.getNumLinks();
	}
}

/** *************************************************************************
 *
 * @throws <std::runtime_error> if copying not possible due to limited max_connectivity
 * */
template<class Vertex,uint max_connectivity>
template<class V, uint MC>
void Connected<Vertex,max_connectivity>::operator=(const Connected<V, MC> &val){
	// copy position
	//static_cast<V>(&this) = static_cast<V>(val);
	this->setAllCoordinates(val.getX(), val.getY(), val.getZ());///@todo call base = operator instead

	//copy connectivity. throw exception if not possible due to limited max_connectivity
	if (this->getMaxConnectivity() < val.getNumLinks()) {
		throw std::runtime_error("class Connected - constructor: Tried to import from connected Vertex with larger connectivity.");
	}
	else {
		for (uint i = 0; i < val.getNumLinks(); ++i) {
			links[i] = val.getNeighborIdx(i);
		}
		counter = val.getNumLinks();
	}
}


/**
 * @throw <std::runtime_error> if no \a j-th connection exist in the graph.
 *
 * @param i The i-th connection of this vertex (monomer).
 * @return Returns the index of the connected vertex (monomer).
 */
template<class Vertex,uint max_connectivity>
uint32_t Connected<Vertex,max_connectivity>::getNeighborIdx(uint32_t i) const {
	if (i < this->getNumLinks())
		return links[i];
	else{
		std::stringstream errormessage;
		errormessage<<"Connected::getNeighborIdx( idx ): idx="<<i<<" is exceeding the number of connected neighbors.";
		throw std::runtime_error(errormessage.str());
	}
}


/**
 * @details It extends the list storing the indices of connected Vertices.
 * It´s also test if maximum connectivity is sufficient.
 * This function does nothing if the link already exists.
 *
 * @throw <std::runtime_error> if index is invalid esp. negative.
 * @param b Index of Vertex (monomer) to connect with.
 *
 */
template<class Vertex,uint max_connectivity>
void Connected<Vertex,max_connectivity>::connect(int32_t b) {

	// negative index is not allowed
	if (b<0)
	{
		throw std::runtime_error("Connected::connect(uint b): b < 0");
	}
	// check whether object b is already connected to me
	for (uint32_t i = 0; i < counter; i++)
	{
		if (((uint32_t)b) == links[i]){
/*#ifdef DEBUG
			std::cerr<< "Connected::connect(unit b): Vertices already connected" << std::endl;
#endif	// End DEBUG
			/*  */
			return;
		}
	}

	checkMaxConnectivity();

	links[counter] = uint32_t (b);
	counter++;
}


/**
 * @details This function checks if specified link to the Vertex with index b exists.
 * It (maybe) changes the position of the i-th link of this Vertex in the internal list.
 * If you using \sa getNeighborIdx , make sure you take the reordering into account.
 *
 * @throw <std::runtime_error> if Vertex with index \a b is not linked to this one.
 *
 * @param b Index of Vertex (monomer) to disconnect from.
 */
template<class Vertex,uint max_connectivity>
void Connected<Vertex,max_connectivity>::disconnect(uint32_t b) {
	uint32_t i = 0;

	//test if link with vertex index b exists
	while (i < counter && links[i] != b) {
		++i;
	}
	if (i == counter)

		throw std::runtime_error("Connected::disconnect(): Tried to disconnect unlinked vertex.");
	else {
		for (; i < counter - 1; ++i) {
			links[i] = links[i + 1];
		}
		--counter;
	}
}


/**
 * @throw <std::runtime_error> if Maximum connectivity is reached.
 */
template<class Vertex,uint max_connectivity>
void Connected<Vertex,max_connectivity>::checkMaxConnectivity() {
	if (counter >= max_connectivity) {
		throw std::runtime_error("Connected :: Maximum connectivity reached.");
	}
}


/****************************************************************************/
//! Specializations for max. connectivity 0 consuming zero memory.
/****************************************************************************/
template<class Vertex>
class Connected<Vertex, 0> : public Vertex
{

public:
	Connected() {
	}
	Connected(const Vertex& vertex) :
			Vertex(vertex) {
	}

	template<class V, uint MC>
	Connected(const Connected<V, MC>& src) :
			Vertex(src) {
		if (src.getNumLinks()) {
			throw std::runtime_error("class Connected - constructor: Tried to import from connected Vertex with larger connectivity.");
		}
	}

	uint getNumLinks() const {
		return 0;
	}
	static uint getMaxConnectivity() {
		return 0;
	}
	uint getNeighborIdx(int i) const {
		throw std::runtime_error("Connected::getNeihborIdx(): Connected < ... , 0 > does not hold any neighbors.");
	}
	void connect(uint b) {
		throw std::runtime_error("Connected::connect(): Connected < ... , 0 > does not hold any neighbors.");
	}
	void disconnect(uint b) {
		throw std::runtime_error("Connected::disconnect(): Connected < ... , 0 > does not hold any neighbors.");
	}

};


/****************************************************************************/
//! Specializations for max. connectivity 1, omitting the counter variable.
/****************************************************************************/
template<class Vertex>
class Connected<Vertex, 1> : public Vertex
{
	int link;
	bool isConnected() const {
		return link >= 0;
	}
public:
	Connected() :
			link(-1) {
	}
	Connected(const Vertex& vertex) :
			Vertex(vertex), link(-1) {
	}

	template<class V, uint MC>
	Connected(const Connected<V, MC>& src) :
			Vertex(src) {
		if (this->getMaxConnectivity() < src.getNumLinks()) {
			throw std::runtime_error("class Connected - constructor: Tried to import from connected Vertex with larger connectivity.");
		} else if (src.getNumLinks())
		{
			link = src.getNeighborIdx(0);
		}
		else {
			link = -1;
		}

	}

	uint getNumLinks() const {
		return this->isConnected();
	}
	static uint getMaxConnectivity() {
		return 1;
	}
	uint getNeighborIdx(int i) const {
		if (i < this->getNumLinks())
			return link;
		else
			throw std::runtime_error("Connected::getNeihborIdx( idx ): 'idx' is exceeding the number of connected neighbors.");
	}
	void connect(uint b) {
		if (!this->isConnected())
			link = b;
		else
			throw std::runtime_error("Connected :: connect() : Maximum connectivity reached.");
	}
	void disconnect(uint b) {
		if (b != link)
			throw std::runtime_error("Connected::disconnect(): Tried to disconnect unlinked vertex.");
		else
			link = -1;
	}

};

////////////////////////////////////////////////////////////////////////////////

#endif
