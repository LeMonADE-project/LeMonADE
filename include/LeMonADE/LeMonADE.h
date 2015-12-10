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

#ifndef LEMONADE_LEMONADE_H
#define LEMONADE_LEMONADE_H

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/core/MoleculesRead.h>
#include <LeMonADE/utility/DistanceCalculation.h>
#include <LeMonADE/utility/NumericFunctions.h>
#include <LeMonADE/utility/NumericTools.h>
#include <LeMonADE/Version.h>

/**
 * @file
 *
 * @namespace LeMonADE
 *
 * @brief The namespace <b>LeMonADe</b> of the Monte-Carlo-Algorithm-project
 *
 * @details The name is an abbreviation of
 * * <b>L</b>attice-based <b>e</b>xtensible <b>Mon</b>te-Carlo <b>A</b>lgorithm and <b>D</b>evelopment <b>E</b>nvironment"
 * 
 * @todo do we use the namespace?
 **/
using namespace Lemonade;

#endif

