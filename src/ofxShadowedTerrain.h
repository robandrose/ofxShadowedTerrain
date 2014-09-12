/*
 *  ofxShadowedTerrain
 *  Terrain container that provides a mesh, shadow images depending on direktional light direction and heightlines
 */


#pragma once
#include "ofMain.h"
#include "contours.h"
#include "ofxOpenCv.h"


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
	
    // Loading:
	void loadMapFromTextfile(string _filename);
    void loadMapFromImage(string _filename);
	
    // Getters:
    ofMesh* getMesh(){return &mesh;};
    ofImage* getShadowImage();
	
    ofVec3f getNormalAt(int x, int y);
	ofVec3f getSmoothNormalAt(int x, int y);
	float* getHeightmap(){return heightmap;};
    int getHeightMapWidth(){return gridw;};
	int getHeightMapHeight(){return gridh;};
	float getCellsize(){return cellsize;};
	
	// Setterns
    void setLightAngles(float xangle, float yangle);
	void setLightDirection(ofVec3f _dir);
	void setHeightlineParameters();
    
    void drawMesh();

	
	// CRaster Funktionen:
	double value(double x,double y);
	SPoint upper_bound();
	SPoint lower_bound();
	
	
    
    
private:	
	
    void loadHeightmapData(string _filename);
    
    
	ofVec3f calculateNormal(int x, int y);
	ofVec3f getVertexFromImg(ofFloatImage& img, int x, int y);
    float getValueFromImagePos(ofFloatImage& img, int x, int y); 
    void addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c);
    void addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c, ofVec3f d);
   
    void prepareForShadows();
    void updateLightMap(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]);
	
    
    ofMesh mesh;
    
    ofFloatImage heightmapimg;
    heightMapData  heightmapdataobj;
    
    // TODO: remove these vars, its all in the struct:
	int ncols, nrows;
	float xllcorner, yllcorner;
	float nodata_value;
	int gridw, gridh;
	float cellsize;
    int numvalues;
    bool heightmaploaded;
	float zstretchfact;
	
    // Normals are also in the mesh, maybe have a smoothed mesh? and get the smoothednormals from there?
	ofVec3f* normals;
	ofVec3f* smoothnormals;
	float* heightmap;
    
    unsigned char *lightmap;
	unsigned char *lightmapforcalc;
	
	ofVec3f lightdir;
	ofColor shadowcolor;
	float shadowalpha;
	ofColor	lightcolor;
	ofImage shadowimg;
	ofxCvGrayscaleImage blurimg;
    

    
	
	
};

