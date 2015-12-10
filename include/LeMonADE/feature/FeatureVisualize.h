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

#ifndef LEMONADE_FEATURE_FEATUREVISUALIZE_H
#define LEMONADE_FEATURE_FEATUREVISUALIZE_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/utility/Vector3D.h>

/*****************************************************************************/
/**
 * @file FeatureVisuaize.h
 * @brief Feature stores information needed for uniform visualization
 * */
/*****************************************************************************/



/*****************************************************************************/
/**
 * @class VisualizeMonomer
 * @brief Holds color, opacity and visibility tag for a single monomer
 * */
/*****************************************************************************/
class VisualizeMonomer {
public:
	VisualizeMonomer() :
		visible(true),color(VectorFloat3(0.0f,0.0f,0.54296875f)), smoothCoordinate(VectorFloat3(0.0f,0.0f,0.0f)), opacity(1.0), radius(1.0) {
	}

	VectorFloat3 getColor() const
	{
		return color;
	}

	void setColor(VectorFloat3 color)
			{
		this->color = color;
	}

	void setColor(float red, float green, float blue)
	{
			this->color.setAllCoordinates(red, green, blue);
	}

	bool isVisible() const
	{
		return visible;
	}

	void setVisible(bool visible)
			{
		this->visible = visible;
	}

	float getOpacity() const
	{
		return opacity;
	}

	void setOpacity(float opacity)
			{
		this->opacity = opacity;
	}

	VectorFloat3 getSmoothCoordinate() const {
		return smoothCoordinate;
	}

	void setSmoothCoordinate(VectorFloat3 smoothCoordinate) {
		this->smoothCoordinate = smoothCoordinate;
	}

	VectorFloat3 getSmoothCoordinateBetween() const {
		return smoothCoordinateBetween;
	}

	void setSmoothCoordinateBetween(VectorFloat3 smoothCoordinateBetween) {
		this->smoothCoordinateBetween = smoothCoordinateBetween;
	}

	float getRadius() const {
		return radius;
	}

	void setRadius(float radius) {
		this->radius = radius;
	}

private:
	bool visible; // monomer is visible or not
	VectorFloat3 color; // color in rgb-code  with r/g/b/ between 0.0...1.0
	float opacity; // opacity between 0.0...1.0

	float radius;

	VectorFloat3 smoothCoordinate;
	VectorFloat3 smoothCoordinateBetween;
};


/*****************************************************************************/
/**
 * @class FeatureVisualize
 * @brief extends monomer type with VisualizeMonomer and holds global visualization parameters
 * */
/*****************************************************************************/
class FeatureVisualize: public Feature {
public:

	enum FoldBackType {ABSOLUTE, CHAINSTARTS_IN_BOX, EVERYTHING_IN_BOX};

	typedef LOKI_TYPELIST_1(VisualizeMonomer)monomer_extensions;

	FeatureVisualize(): cameraPosition(VectorInt3(1000,1000,1000)), visualizeBonds(true), foldBack(EVERYTHING_IN_BOX), isSmoothingCoordinates(false), maxNumSmoothingFrame(2), Translation(VectorInt3(0,0,0)), drawingMonomersAsSpheres(false), subdivisionSpheres(4) {}
	virtual ~FeatureVisualize() {}

public:
	VectorInt3 getCameraPosition() const
	{
		return cameraPosition;
	}

	void setCameraPosition(VectorInt3 cameraPosition)
			{
		this->cameraPosition = cameraPosition;
	}

	FoldBackType getFoldBack() const
	{
		return foldBack;
	}

	void setFoldBack(FoldBackType foldBack)
			{
		this->foldBack = foldBack;
	}

	bool isVisualizeBonds() const
	{
		return visualizeBonds;
	}

	void setVisualizeBonds(bool visualizeBonds)
			{
		this->visualizeBonds = visualizeBonds;
	}

	bool isSmoothing() const
	{
		return isSmoothingCoordinates;
	}

	void setSmoothing(bool isSmoothing)
	{
		isSmoothingCoordinates = isSmoothing;
	}

	uint8_t getMaxNumSmoothingFrame() const
	{
		return maxNumSmoothingFrame;
	}

	void setMaxNumSmoothingFrame(uint8_t maxNumSmoothingFrame)
	{
		this->maxNumSmoothingFrame = maxNumSmoothingFrame;
	}

	void setTranslationInX(int32_t transInX) {
		Translation.setX(transInX);
	}

	void setTranslationInY(int32_t transInY) {
			Translation.setY(transInY);
		}

	void setTranslationInZ(int32_t transInZ) {
				Translation.setZ(transInZ);
			}

	VectorInt3 getTranslation() const {
		return Translation;
	}

	int32_t getTranslationInX() const {
			return Translation.getX();
		}

	int32_t getTranslationInY() const {
				return Translation.getY();
			}

	int32_t getTranslationInZ() const {
					return Translation.getZ();
				}

	bool isDrawingMonomersAsSpheres() const {
		return drawingMonomersAsSpheres;
	}

	void setDrawingMonomersAsSpheres(bool drawingMonomersAsSpheres) {
		this->drawingMonomersAsSpheres = drawingMonomersAsSpheres;
	}

	uint8_t getSubdivisionSpheres() const {
		return subdivisionSpheres;
	}

	void setSubdivisionSpheres(uint8_t subdivisionSpheres) {

		if(subdivisionSpheres > 15)
			subdivisionSpheres = 15;

		if(subdivisionSpheres < 1)
			subdivisionSpheres = 1;

		this->subdivisionSpheres = subdivisionSpheres;
	}

private:
	VectorInt3 cameraPosition;
	bool visualizeBonds; // switch visualization on or off
	FoldBackType foldBack; // fold back monomers type, can either be ABSOLUTE (all monomers are displayed at their absolute postion),
	//CHAINSTARTS_IN_BOX (all chain starts are in the box) or EVERYTHING_IN_BOX(amm monomers are in the box)

	bool isSmoothingCoordinates;
	uint8_t maxNumSmoothingFrame;

	VectorInt3 Translation;

	bool drawingMonomersAsSpheres;

	uint8_t subdivisionSpheres;
};


#endif /* LEMONADE_FEATURE_FEATUREVISUALIZE_H */
