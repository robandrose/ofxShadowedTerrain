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
#include "ofMain.h"
#include "contours.h"

#define NUM_GRID_POINTS 1024*1024*2

struct heightMapData{
    
    int ncols, nrows;
    float xllcorner, yllcorner;
    float nodata_value;
    int gridw, gridh;
    float cellsize;
    
    vector<float>mydata;
    
    ofVec3f getVertexAt(int x, int y){
        int index=x+(y*ncols);
        if(index>=0 && index< mydata.size()){
            return ofVec3f(x, y,mydata.at(index));
        }return ofVec3f(x,y,0);
    }
};


class ofxShadowedTerrain: public CRaster{

    
    
public:
	
	ofxShadowedTerrain();
	virtual ~ofxShadowedTerrain();
	
	void loadMapFromTextfile(string _filename);
    void loadMapFromTextfileOld(string _filename);
    void loadMapFromImage(string _filename);
	
    
	ofVec3f getNormalAt(int x, int y);
	ofVec3f getSmoothNormalAt(int x, int y);
	
    
	float* getHeightmap(){return heightmap;};
	
	int getHeightMapWidth(){return gridw;};
	int getHeightMapHeight(){return gridh;};
	float getCellsize(){return cellsize;};
	
	void setLightAngles(float xangle, float yangle);
	void setLightDirection(ofVec3f _dir);
	
	void drawForTexture();
	void drawForColor();
    void drawMesh();


	ofImage* getShadowImage();
	
	// CRaster Funktionen:
	double value(double x,double y);
	SPoint upper_bound();
	SPoint lower_bound();
	
	
    
     ofMesh mesh;
    
private:	
	
	void updateLightMap(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]);
	void updateLightMapFast(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]);
    
    
	ofVec3f calculateNormal(int x, int y);
	ofVec3f getVertexFromImg(ofFloatImage& img, int x, int y);
    void addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c);
    void addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c, ofVec3f d);
    
    
    void loadHeightmapData(string _filename);
    
    
    
    
    
	// ASC VARS	
	int ncols, nrows;
	float xllcorner, yllcorner;
	float nodata_value;
	int gridw, gridh;
	float cellsize;	
	
    
    
    
	vector<float>heightmapdata;
	
	ofVec3f* normals;
	ofVec3f* smoothnormals;
	float* heightmap;
    
	unsigned char *lightmap;
	unsigned char *lightmapforcalc;
	
	int numvalues;
	
	ofVec3f lightdir;
	
	bool heightmaploaded;
	
	ofColor shadowcolor;
	ofColor	lightcolor;
	ofImage shadowimg;	
	
	float shadowalpha;
	float zstretchfact;
    
    
    
    
   
    ofFloatImage heightmapimg;
    heightMapData  heightmapdataobj;
    


    
	
	
};

