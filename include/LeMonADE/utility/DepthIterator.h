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

#ifndef LEMONADE_UTILITY_DEPTHITERATOR_H
#define LEMONADE_UTILITY_DEPTHITERATOR_H

#include <deque>
#include <set>
#include <stack>

#include <LeMonADE/utility/DepthIteratorPredicates.h>

/******************************************************************************/
/**
 * @file
 *
 * @class GraphIteratorDepthFirst
 *
 * @brief Iterates over a given graph (e.g. Molecules <...>) with priority to follow branches into depth under arbitrary conditions (vertex predicates).
 *
 * @tparam GraphType Class of the graph to iterate over, for instance a Molecules <...> class.
 * GraphType is required to have an index-based interface for connectivity:
 * getNumLinks(int i) and getNeighborIdx(int i,int j).
 * Furthermore a member size() returning the number of vertices is used and the operator[] to access a vertex.
 * The expected behavior of GraphType can be seen in Molecules < ... > as an example.
 *
 * @tparam Predicate This user-specified class defines, which vertices are allowed to be visited.
 * Predicate is required to have an operator() (const GraphType& , int ) returning bool.
 * This operator should return "True" for a pair of graph object and vertex index,
 * if the vertex is allowed to be visited and "False" otherwise.
 *
 * @todo we should reconsider this approach for usability
 **/
template < class GraphType, class Predicate = alwaysTrue > class GraphIteratorDepthFirst
{

  struct IteratorPosition;

public:

/**
 * @brief Constructor.
 * @param graph The graph object to iterator through
 * @param excluded exclude vertices to be visited by a set of integers. If not NULL, this is handled as the set of already visited vertices (also further filled up during iteration).
 * @param pred This can be used to provide a user-created Predicate object (e.g., if more information is needed for evaaluation).
 *
 **/
  GraphIteratorDepthFirst(const GraphType& graph, std::set < int >* excluded = 0, const Predicate& pred = Predicate() )
  :graph(graph)
  ,myPosition(0,0)
  ,predicate(pred)
  ,visited(excluded)
  ,owning_visited(false)
  {
    if ( visited == 0) {
      visited = new std::set < int >;
      owning_visited = true;
    }

    int i = 0; while ( i < graph.size() && !allowed_to_visit(i) ) { ++i; }

    myPosition.currentVertex = i;
    myPosition.currentNeighbor = 0;

    if ( !this->isEnd() ) visit_current() ;

  };

  ///@brief Destructor.
  ~GraphIteratorDepthFirst()
  {
    if ( owning_visited && visited != 0 )
    {
      delete visited; visited = 0;
    }
  }

  /**
   * @brief Advances the iterator.
   *
   * @todo return value is missing
   * @return
   **/
  GraphIteratorDepthFirst& operator++()
  {

	// The following in simple words:
	//
	// If the end has been reached, there is nothing to iterate further.
	// Else:
	//   If one of the neighbors of the current vertex is allowed to be visited: Visit this neighbor.
	//   Else:
	//     Look for return points:
	//     If: there are unvisited branches left: Go to one of them.
	//     Else: there is nothing to be done any more, We label ourselves to be at the end.

      while(!this->isEnd())
      {
	for (; myPosition.currentNeighbor < graph.getNumLinks(myPosition.currentVertex); ++myPosition.currentNeighbor)
	{
	      // Get the vertex idx of the "myPosition.currentNeighbor"-th neighbor of the current vertex.
	      int neighbor = graph.getNeighborIdx(myPosition.currentVertex, myPosition.currentNeighbor );

	      if ( allowed_to_visit(neighbor) )
	      {
		      // If this is not the last "free" neighbor of the current vertex,
		      // we add a "return point" to iterate over remaining neighbors later.
		      if ( myPosition.currentNeighbor+1 < graph.getNumLinks(myPosition.currentVertex) ) return_points.push(myPosition);

		      // Update the current state, i,e, jumping to the neighbor.
		      myPosition.currentVertex = neighbor;
		      myPosition.currentNeighbor = 0;

		      // Exclude the new vertex from visits in the future.
		      visit_current();

		      return *this;
	      }
	}


	if (return_points.empty())
	{
		//If there is no return point left, we prepare ourselves for the end. (this->isEnd() will be true from now on).
		myPosition = IteratorPosition(graph.size(),0);
		return *this;
	}
	else
	{
		//If there is a return point left, continue iteration from there.
		myPosition = return_points.top();
		return_points.pop();
	}
      }
  }

  //! Gives the index of the current vertex, which is visited in the moment.
  const int& getVertexIdx() const { return myPosition.currentVertex; }

//   ///@brief Derefers the current vertex index returning a reference to the actual vertex object. (Think of (*this)->someVertexMemberFunction() ).
//   const typename GraphType::vertex_type& operator*() const { return graph.at(myPosition.currentVertex); }

  //! Returns "True" if iteration cannot advance further.
  bool isEnd() const { return myPosition.currentVertex >= graph.size(); }

private:

  //! Reference to the graph object, over which is iterated.
  const GraphType&	graph;

  //! Internal copy of a Predicate (conditional exclusion of vertices), used in allowed_to_visit() method.
  Predicate 		predicate;

  //! The current state of the iterator.
  IteratorPosition 		myPosition;

  //! Return points oif the recursion, i.e. points in the history, where a new branch might be started.
  std::stack < IteratorPosition >	return_points;

  //! Holds indiced of all visited vertices (not to be visited again) - and additionally user-given excluded vertices (see constructor).
  std::set < int >*		visited;

  //! Should be true, if the user-given "excluded" pointer was NULL --> Then the visited set should be deleted upon destruction.
  bool 			owning_visited;


  /**
   * @brief Position of the iterator consisting in a vertex index and the neighbor pointing to.
   * Objects of this class are used to hold the actual state of the iterator as well as return points for the recursion.
   */
  struct IteratorPosition
  {
	//! Standard constructor.
	IteratorPosition(int i, int n):currentVertex(i),currentNeighbor(n){}

	//! The index of the current vertex, which is visited in the moment.
	int currentVertex;

	/**
	 * @brief The neighbor, which is aimed to be visited at next iteration.
	 *
	 * @details currentNeighbor should have values between 0 and graph.getNumLinks(currentVertex).
	 * The actual vertex index of the neighbor will be retrieved by using graph.getNeighborIdx( currentVertex, currentNeighbor ).
	 *
	 **/
	int currentNeighbor;
  };


  //! Adds the current vertex (i.e. myPosition.currentVertex) to the set of already visited vertices.
  void visit_current() { visited->insert(myPosition.currentVertex);}

  //! Returns true, if "i" has not been visited and the user has not excluded this vertex manually (see constructor).
  bool is_not_visited( int i ) { return visited->find( i ) == visited->end(); }

  //! Returns true, if Vertex "i" is allowed to be visited: has not been visited; predicate returns "True"; user has not excluded it explicitly (see constructor).
  bool allowed_to_visit( int i) { return is_not_visited(i) && predicate(graph,i); }

};


// template < class Graph, template <class,class> class IndexContainer = std::vector   > class VertexGroup
// {
//
// private:
//
//   IndexContainer< int, std::allocator<int> > indices;
//   const Graph& graph;
//
// public:
//
//   typedef typename Graph::vertex_type vertex_type;
//
//   VertexGroup(const Graph& graph):graph(graph){}
//
//   void push_back(int idx){
//     indices.push_back(idx);
//   }
//
//   size_t size() const { return indices.size(); }
//
//   const vertex_type& operator[] (int idx) const
//   {
//     return graph[indices[idx]];
//   }
// };

/**
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 * @todo comment!!!
 **/
template < class Graph, class GroupContainer, class GroupType, class Predicate >
void fill_connected_groups ( const Graph& g , GroupContainer& groups, const GroupType& group_init, const Predicate& pred )
{

  std::set < int > visited;
//   deque < int > tmp_group;

  bool end = false;
  while (!end)
  {

    GraphIteratorDepthFirst< Graph, Predicate > iter( g, &visited, pred );

    if ( iter.isEnd() )
    {
      end = true;
      break;
    }

    groups.push_back ( group_init );

    do
    {
      groups.back().push_back(iter.getVertexIdx());
      ++iter;
    }
    while ( !iter.isEnd() );
  };
}

/**
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 * @todo comment!!!
 **/
template < class Molecules >
void connectMonosWhichAreClose( Molecules& m){

	for (int i=0; i< m.size(); ++i){
		for (int j=0; j< m.size(); ++j){

			double distance;

			double dx = ((m[i].getX())&255) - ((m[j].getX())&255);
					//cout << "dx: " << dx << endl;

					double dy = ((m[i].getY())&255) - ((m[j].getY())&255);
					//cout << "dy: " << dy << endl;
					double dz = ((m[i].getZ())&255) - ((m[j].getZ())&255);
					//cout << "dz: " << dz << endl;
					distance =  sqrt( dx*dx + dy*dy + dz*dz );



			if ( (distance < 4) && (distance != 0) ){   //if ( (m[i].getDistanceTo(m[j]) < 4) && (m[i].getDistanceTo(m[j]) != 0) ){
				m.connect(i,j);

			}

		}

	}

	return;


}

#endif

