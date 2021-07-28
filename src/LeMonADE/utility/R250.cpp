/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015,2021 by
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
	
		printState();

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
	
		std::cout << "loaded r250 with given values..."<< std::endl;
	
		printState();
		//shuffle array and initialize pointers
		pos=&array[250];
		other147=(pos-147);
		other250=(pos-250);
		refresh();
}

/**
 * this resets the state array to a predefined value, so that a defined state
 * can always be obtained if needed. the numbers were generated from /dev/urandom
 * once and then simply copied here. (replaced 2021 in Release 2.2.2)
 * The array can e.g. be generated with:
 *   head -c 1024 /dev/urandom | hexdump -v -e '4 4 "%u, " "\n"'
 */
void R250::loadDefaultState()
{
	static uint32_t const tmp[R250_RANDOM_PREFETCH] =
    {
        2183635719, 2190244759, 1838922374, 2492928253,
        1691620525, 4174268300, 2925036378, 1686150634,
        4143523061, 92006323, 4177728194, 1411981943,
        1245853850, 1521440829, 2963330919, 856290965,
        2647997580, 2981053599, 1099760386, 3284599097,
        1153699365, 604329032, 2976899913, 3519137272,
        4184227608, 2072347806, 731765320, 3731617640,
        2011825539, 3519091383, 777911520, 3593825682,
        2213846057, 2818450139, 1518574725, 3087387789,
        2839883236, 2917794674, 911197067, 211765470,
        428397452, 2446693707, 921325450, 135181643,
        2748792785, 1430813481, 2819005913, 1794303866,
        1684166381, 2236579585, 183216389, 22600318,
        424653856, 3742353884, 3645276761, 2061915755,
        3915077421, 600520186, 556839987, 734918579,
        2120203728, 2899243738, 1481433027, 2344617992,
        1599513085, 1398259692, 1670606037, 648503549,
        1539868690, 128537250, 3997202055, 2749540959,
        327991835, 2963440711, 3023804961, 4077646674,
        1658275422, 2090997227, 3780061296, 2729972273,
        1926252795, 3202049657, 3879207852, 1900755841,
        2480152485, 561994355, 1380333390, 2134798066,
        2128088413, 14348827, 1584549315, 2472893490,
        2427233072, 3471354632, 3092802058, 3537096898,
        3131526979, 600558991, 729675831, 712609463,
        2977016554, 3732375413, 2921287250, 2698425775,
        57604637, 2706609512, 632385541, 2012671337,
        2010953452, 824814095, 1619673816, 767630393,
        195784944, 662004853, 3374593664, 1653069366,
        1323411388, 2229537257, 2587783893, 447703779,
        1687964944, 3616195188, 4227645990, 3873208902,
        2630042627, 718425856, 706712992, 2272163330,
        2123548242, 3052380041, 1286900347, 3030225619,
        3990904814, 460716308, 685339369, 834928970,
        2129756172, 495105944, 4161716583, 4257896697,
        1728169764, 777344186, 3790271563, 569734115,
        4136669083, 3156717764, 1113962257, 1250225371,
        666738466, 3660168312, 2582701982, 2927556822,
        2309677514, 4228759946, 4275226748, 1917759123,
        2301422716, 3325595869, 3600918959, 1507742373,
        2189228943, 1354410878, 1893866291, 3405482202,
        258140360, 3412838799, 789120517, 1774390742,
        1276112849, 2971021240, 661478379, 2366454877,
        1350258451, 3269346974, 1568983620, 366528176,
        1512560683, 428458871, 2837985004, 1145283036,
        3488172914, 2830841878, 578995039, 1783082104,
        763626672, 3609683706, 997748244, 291508732,
        1082754828, 1793724907, 445425176, 306544940,
        1141049923, 2340771713, 2695001510, 1039741825,
        25338870, 1961306592, 1970584190, 1932328866,
        70484709, 2794213096, 1947251298, 2384361454,
        1532598777, 1413594045, 3741770196, 3341090742,
        3425538499, 1077106632, 679120682, 548535452,
        1785460316, 2833553614, 2641172924, 3482103930,
        2922868348, 3452689737, 778677809, 1852367327,
        3185222054, 2475234521, 1105940197, 188880018,
        3701380499, 4075520071, 3197904393, 3346375591,
        2092415115, 3541136001, 1609674945, 2334264928,
        442922817, 3551460990, 96703578, 355623315,
        1838176107, 2630741652, 2614239307, 947656070,
        1720998510, 1260963745, 412104712, 2897568188,
        504494002, 409223266, 13584767, 4006567101,
        1443312917, 1985791450, 1428736195, 3128819800,
        4106263747, 2661162000, 2440997736, 3111816589
    };
	for (size_t i=0;i<R250_RANDOM_PREFETCH;i++)
	{
		array[i]=tmp[i];
	}
	
	std::cout << "size of rng = " << sizeof(array) << std::endl;
    	std::cout << "loaded r250 default state..."<< std::endl;
    	printState();

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
	std::cout << std::endl << std::endl;
}
