/*
 *  TerrainDataProvider.h
 *  LandscapeShadow4
 *
 *  Created by Matthias Rohrbach on 30.09.10.
 *  Copyright 2010 rob & rose grafik. All rights reserved.
 *
 
Eigentlich das Model für das Terrain, weiss alles, kann alles.
 
 */

#pragma once

#include "basicScreenObject.h"
#include "ofxDisplayList.h"
#include "ofxFBOTexture.h"

#include "container.h"
#include "Image.h"

#include "contours.h"
#include "ofxOpenCv.h"

#define NUM_GRID_POINTS 1024*1024*2

class TerrainDataProvider: public BasicScreenObject, public CRaster{
	
public:
	
	TerrainDataProvider();
	virtual ~TerrainDataProvider();
	
	void loadHeightMap(string _filename);
	
	void firstUpdate();
	void update();
	void _draw();
	
	
	ofxVec3f getNormalAt(int x, int y);
	ofxVec3f getSmoothNormalAt(int x, int y);
	
	
	
	float* getHeightmap(){return heightmap;};
	
	int getHeightMapWidth(){return gridw;};
	int getHeightMapHeight(){return gridh;};
	float getCellsize(){return cellsize;};
	
	void setLightAngles(float xangle, float yangle);
	void setLightDirection(ofxVec3f _dir);	
	
	void drawForTexture();
	void drawForColor();


	ofImage* getShadowImage();
	
	// CRaster Funktionen:
	double value(double x,double y);
	SPoint upper_bound();
	SPoint lower_bound();
	
	
private:	
	
	void updateLightMap(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]);
	void updateLightMapFast(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]);
	
	ofxVec3f calculateNormal(int x, int y);	
	
	// ASC VARS	
	int ncols, nrows;
	float xllcorner, yllcorner;
	float nodata_value;
	int gridw, gridh;
	float cellsize;	
	
	vector<float>heightmapdata;
	
	ofxVec3f* normals;
	ofxVec3f* smoothnormals;
	
	float* heightmap;
	unsigned char *lightmap;
	unsigned char *lightmapforcalc;
	
	int numvalues;
	
	ofxVec3f lightdir;
	
	ofxDisplayList mylist;
	bool heightmaploaded;
	
	ofColor shadowcolor;
	ofColor	lightcolor;
	ofImage shadowimg;	
	
	float shadowalpha;
	float zstretchfact;
	
	
	
	
};

