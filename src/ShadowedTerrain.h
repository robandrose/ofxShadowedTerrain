/*
 *  ShadowedTerrain.h
 *  LandscapeShadow2
 *
 *  Created by Matthias Rohrbach on 26.08.10.
 *  Copyright 2010 rob & rose grafik. All rights reserved.
 *
 */

#pragma once

#include "basicScreenObject.h"
#include "ofxDisplayList.h"
#include "ofxFBOTexture.h"

#include "container.h"
#include "Image.h"
#include "TerrainDataProvider.h"
#include "contours.h"


class ShadowedTerrain: public BasicScreenObject{
	
public:
	
	ShadowedTerrain();
	virtual ~ShadowedTerrain();
	
	void setup();
	void firstUpdate();
	void update();
	void _draw();
	
	
	void isConfigMode(bool _isconfigmode);
	bool isConfigMode(){return isconfigmode;};
	
	float* getHeightmap(){return terraindata->getHeightmap();};
	int getHeightMapWidth(){return terraindata->getHeightMapWidth();};
	int getHeightMapHeight(){return terraindata->getHeightMapHeight();};
	float getCellsize(){return terraindata->getCellsize();};
	
	void setShadowAlpha(float _shadalpha);	
	void setLightAngles(float xangle, float yangle);
	void setLightDirection(ofxVec3f _dir);
	
	
	
	Container fbocontainer;
	Container container;
	
	TerrainDataProvider* terraindata;
	
private:
	void drawToList();

	ofxDisplayList mylist, diagnosticList;
	ofxFBOTexture fbo;
	
	float shadowalpha;
	bool isconfigmode;
		
	
};