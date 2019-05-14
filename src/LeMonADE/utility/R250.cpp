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

#include <LeMonADE/utility/R250.h>

/**
 * @file
 * @brief implementation of R250Engine
 * */

/// constructor
//pointers pos,other147 and other250 are set to positions in the array
//according to random number generation algorithm.
//loads predefined state.
R250::R250():
	   arrayEnd(&array[R250_RANDOM_PREFETCH]),
	   pos(&array[250]),
	   other147(pos-147),
	   other250(pos-250)
{
		loadDefaultState();
}


//applies the random number algorithm to the internal state array
//the algorithm applies bitwise ^ operations to the numbers in the array, such
//that new pseudo random numbers are generated. the poiner dice is then set to
//the beginning of the array, from where the new numbers are drawn.
void R250::refresh()
{
	do
	{
	  *(pos++)		= *(other147++)^(*(other250++));
	}
	while(pos != arrayEnd);
	pos			= array;

	do
	{
	  *(pos++)		= *(other147++)^(*(other250++));
	}
	while(other147!=arrayEnd);

	other147			= array;
	do
	{
	  *(pos++)		= *(other147++)^(*(other250++));
	}
	while(other250!=arrayEnd);

	dice 					= array;
	other250			= array;

}

//initializes the internal state array from /dev/urandom
void R250::loadRandomState()
{
		std::cout << "size of rng = " << sizeof(array) << std::endl;
		std::cout << "loading r250 state from /dev/urandom...";

		std::ifstream in("/dev/urandom");
		in.read((char*)array,sizeof(array));

		std::cout << "ready.\n";

		//shuffle array and initialize pointers
		pos=&array[250];
		other147=(pos-147);
		other250=(pos-250);
		refresh();
}

//this allows to set the internal state array to given values.
void R250::setState( uint32_t const * stateArray )
{
		for(size_t i=0;i<R250_RANDOM_PREFETCH;i++)
		{
			array[i]=stateArray[i];
		}
		//shuffle array and initialize pointers
		pos=&array[250];
		other147=(pos-147);
		other250=(pos-250);
		refresh();
}

/**
 * this resets the state array to a predefined value, so that a defined state
 * can always be obtained if needed. the numbers were generated from /dev/urandom
 * once and then simply copied here.
 * The array can e.g. be generated with:
 *   head -c 1024 /dev/urandom | hexdump -v -e '4 4 "0x%08X, " "\n"'; echo
 */
void R250::loadDefaultState()
{
	static uint32_t const tmp[R250_RANDOM_PREFETCH] =
    {
        128861083 ,1578465172 ,1132109222 ,1536640935 ,488265105 ,249545814 ,
		931473061 ,816964859 ,296369473 ,818986032 ,530371272 ,295126204 ,
		1516346231 ,3026009 ,245630516 ,671061864 ,791099058 ,852977361 ,
		635292087 ,1059222032 ,288528100 ,877554999 ,1941668641 ,407091003 ,
		1633608675 ,87217833 ,1372169418 ,1162729165 ,1579988388 ,25035216 ,
		267956911 ,424734689 ,672799208 ,1461350327 ,905462485 ,467686926 ,
		1072513216 ,932097491 ,1656505938 ,815809826 ,389301423 ,1083577174 ,
		1496526548 ,1568324204 ,1631549701 ,527246014 ,1660056902 ,1892593678 ,
		458159066 ,989601043 ,1370765921 ,1486161662 ,315061300 ,752007777 ,
		172558206 ,1314025472 ,239686405 ,2074002083 ,1794956553 ,929472681 ,
		722167620 ,1639077828 ,1342901080 ,1390526957 ,850546521 ,1840137370 ,
		2010025976 ,1342378906 ,362570076 ,1184522593 ,278059745 ,245252898 ,
		120920959 ,2122315972 ,1801249454 ,1624621127 ,2126148557 ,689581345 ,
		886719575 ,197326348 ,1692832343 ,1691880553 ,1121407416 ,1231431139 ,
		540550150 ,924384068 ,446830463 ,663805735 ,11964822 ,440537339 ,
		1800879640 ,1530404315 ,938488208 ,851126733 ,903853048 ,27314800 ,
		2003780392 ,443160143 ,244208778 ,427777353 ,447520999 ,1453776985 ,
		424274945 ,1928014498 ,1532427678 ,443111457 ,1020877813 ,1420058598 ,
		426111972 ,1836100456 ,1547943653 ,331276469 ,546494085 ,1488526425 ,
		649850648 ,109698414 ,49682325 ,583184179 ,1561190358 ,1088257452 ,
		2111576298 ,2144830684 ,1160934926 ,1692925416 ,438197862 ,491798340 ,
		727308436 ,594748985 ,1342345046 ,925081293 ,1691764903 ,843310425 ,
		1200838282 ,1444721217 ,2128110335 ,60522306 ,95354429 ,2026845056 ,
		1369901402 ,753733242 ,719461246 ,1203997431 ,780174478 ,806089039 ,
		371934522 ,1925937855 ,505485435 ,1827061014 ,37806837 ,148646884 ,
		1292540716 ,2059142989 ,394384341 ,65226295 ,934617354 ,745940390 ,
		614072781 ,601418493 ,613537240 ,1156936041 ,895410893 ,224492918 ,
		648568013 ,1329964435 ,1798116056 ,1257249331 ,974027352 ,642549434 ,
		1093233997 ,1538944681 ,1508303997 ,90457605 ,925894041 ,105072093 ,
		1537361614 ,544371048 ,624170083 ,1986381084 ,1840135544 ,923836043 ,
		340434338 ,248360032 ,442370652 ,1261015906 ,868832477 ,301762886 ,
		1845376447 ,1037773296 ,961432347 ,1331157624 ,953366723 ,285066463 ,
		125584791 ,1280701633 ,2059297229 ,1898695339 ,422570809 ,1421759215 ,
		1042390228 ,1940638910 ,298849264 ,882619357 ,1502641883 ,994170722 ,
		758164584 ,2059360332 ,1110464839 ,1950696339 ,538103478 ,263502002 ,
		1905740234 ,684347334 ,1868238890 ,1788991082 ,645021092 ,456474838 ,
		1572869344 ,197764990 ,598604653 ,737670288 ,1827143439 ,466334365 ,
		265970336 ,1197591797 ,1604053999 ,732466110 ,1843046562 ,2111327269 ,
		1902834539 ,634999463 ,303074817 ,184445115 ,1094921037 ,357510967 ,
		866089755 ,313184427 ,619508128 ,759157274 ,627786029 ,640921970 ,
		5035326 ,351442498 ,310023587 ,713733322 ,1567488348 ,225865858 ,
		1379860041 ,1091951759 ,952392682 ,650949040 ,2000007399 ,82102099 ,
		1654883130 ,394023729 ,2069457331 ,1270363068
    };
	for (size_t i=0;i<R250_RANDOM_PREFETCH;i++)
	{
		array[i]=tmp[i];
	}

	//shuffle
	pos=&array[250];
	other147=(pos-147);
	other250=(pos-250);
	refresh();
}

//prints the internal state array
void R250::printState()
{
	for (size_t i = 0; i<R250_RANDOM_PREFETCH; i++)
	{
		std::cout << "\n array["
			 << std::setw(4)  << std::right << i << "] = "
			 << std::setw(14) << std::right << array[i] << "; ";
	}
}
