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

/**
 * @todo comment
 **/

#ifndef TL_EXT_H
#define TL_EXT_H

namespace Loki
{
namespace TL
{
	template <class Head, class Tail>
	struct TypeAt<Typelist<Head, Tail>, -1>
	{
		typedef NullType Result;
	};

	template <class H , class T >                         // Specialization 1
	struct Erase<Typelist<H,T>, NullType>
	{
		typedef Typelist<H,T> Result;
	};

	template <class TList> struct NoDuplicatesReverse;

	template <> struct NoDuplicatesReverse<NullType>
	{
		typedef NullType Result;
	};

	template <class Head, class Tail>
	struct NoDuplicatesReverse< Typelist<Head, Tail> >
	{
		private:
		typedef typename NoDuplicatesReverse<Tail>::Result NewTail;
		typedef Typelist<Head, NewTail> TmpTlist;

		enum { idx = ( IndexOf < NewTail, Head >::value == -1) ? -1 : 0 };

		public:
		typedef typename Erase < TmpTlist, typename TypeAt < TmpTlist, idx >::Result >::Result Result;
	};

// 	template < class Tlist > struct AppendIfNotNull;

	template < class T, class H > struct AppendIfNotNull //< Typelist< T,H > >
	{
		typedef typename Append< H, T > ::Result Result;
	};

	template < class T > struct AppendIfNotNull < T, NullType >
	{
		typedef T Result;
	};
};
};

#endif
