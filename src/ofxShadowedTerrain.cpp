/*
 *  ofxShadowedTerrain.cpp
 *  LandscapeShadow4
 *
 *  Created by Matthias Rohrbach on 30.09.10.
 *  Copyright 2010 rob & rose grafik. All rights reserved.
 *
 */

#include "ofxShadowedTerrain.h"
#include <fstream>

ofxShadowedTerrain::ofxShadowedTerrain(){
	heightmaploaded=false;
	
	shadowcolor.r=0;
	shadowcolor.g=0;
	shadowcolor.b=50;
	shadowcolor.a=255;
	
	lightcolor.r=255;
	lightcolor.g=255;
	lightcolor.b=255;
	lightcolor.a=255;
    
	zstretchfact=2;
    
    skipx=3;
    skipy=3;
}

ofxShadowedTerrain::~ofxShadowedTerrain(){
	
}



//***************************************************************** // LOADING:

void ofxShadowedTerrain::setSkipValues(int _skipx, int _skipy){
    skipx=_skipx;
    skipy=_skipy;
    
}

void ofxShadowedTerrain::setStretchfactor(float _stretchfactor){
    zstretchfact=_stretchfactor;
}

void ofxShadowedTerrain::loadMapFromImage(string _filename){
    heightmapimg.loadImage(_filename);
	
    mesh.setMode(OF_PRIMITIVE_TRIANGLES);

	gridw = heightmapimg.getWidth();
	gridh = heightmapimg.getHeight();
	
    int rawindex=0;
    
    heightmapdataobj.ncols=heightmapimg.getWidth();
    heightmapdataobj.nrows=heightmapimg.getHeight();
    
    
    for(int y = 0; y < gridh; y += 1) {
		for(int x = 0; x < gridw ; x += 1) {
            rawindex=x+y*gridw;
            heightmapdataobj.mydata.push_back(getValueFromImagePos(heightmapimg, x, y));
        }
    }
    
    
    triangles.clear();
    
    for(int y = 0; y < gridh - skipy; y += skipy) {
		for(int x = 0; x < gridw - skipx; x += skipx) {
			ofVec3f nw = getVertexFromImg(heightmapimg, x, y);
			ofVec3f ne = getVertexFromImg(heightmapimg, x + skipx, y);
			ofVec3f sw = getVertexFromImg(heightmapimg, x, y + skipy);
			ofVec3f se = getVertexFromImg(heightmapimg, x + skipx, y + skipy);
			
            addFace(mesh, nw, ne, se, sw);
		}
	}
    mesh.setFromTriangles(triangles);

    prepareForShadows();
}


ofVec3f ofxShadowedTerrain::getVertexFromImg(ofFloatImage& img, int x, int y) {
    return ofVec3f(x, y, getValueFromImagePos(img, x, y));
}

float ofxShadowedTerrain::getValueFromImagePos(ofFloatImage& img, int x, int y){
    return zstretchfact * img.getColor(x,y).getBrightness();
}

void ofxShadowedTerrain::addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c) {
	ofVec3f normal = ((b - a).cross(c - a)).normalize();
    
    ofMeshFace face;
    face.setVertex(0, a);
    face.setNormal(0, normal);
    face.setTexCoord(0, a);
    
    
    face.setVertex(1, b);
    face.setNormal(1, normal);
    face.setTexCoord(1, b);
    
    face.setVertex(2, c);
    face.setNormal(2, normal);
    face.setTexCoord(2, c);
    
    triangles.push_back(face);
}

void ofxShadowedTerrain::addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c, ofVec3f d) {
	addFace(mesh, a, b, c);
	addFace(mesh, a, c, d);
}


// LOAD FROM TEXTFILE:

void ofxShadowedTerrain::loadMapFromTextfile(string _filename){
    
    loadHeightmapData(_filename);
    mesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
	gridw = heightmapdataobj.ncols;
	gridh = heightmapdataobj.nrows;
    
    triangles.clear();
    
    for(int y = 0; y < gridh - skipy; y += skipy) {
		for(int x = 0; x < gridw - skipx; x += skipx) {
            ofVec3f nw= heightmapdataobj.getVertexAt(x, y);
            ofVec3f ne=heightmapdataobj.getVertexAt(x+skipx, y);
			ofVec3f sw=heightmapdataobj.getVertexAt(x, y+skipy);
			ofVec3f se=heightmapdataobj.getVertexAt(x+skipy, y+skipy);
            addFace(mesh, nw, ne, se, sw);
		}
	}
    
    mesh.setFromTriangles(triangles);
    prepareForShadows();
}


void ofxShadowedTerrain::loadHeightmapData(string _filename){
    ifstream data;
	string line;
	float zraw;
    heightmapdataobj.mydata.clear();
	
    data.open(ofToDataPath(_filename,true).c_str());
	while(!data.eof()){
		getline(data, line);
		vector<string>linevec=ofSplitString(line, " ");
		if(linevec.size()>0){
			if(linevec.at(0)=="NCOLS"){
				heightmapdataobj.ncols=atoi(linevec[1].c_str());
			}else if(linevec.at(0)=="NROWS"){
				heightmapdataobj.nrows=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="XLLCORNER"){
				heightmapdataobj.xllcorner=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="YLLCORNER"){
				heightmapdataobj.yllcorner=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="CELLSIZE"){
				heightmapdataobj.cellsize=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="NODATA_VALUE"){
				heightmapdataobj.nodata_value=atoi(linevec[1].c_str());
			}
			else{
                
				for(int i=0;i<linevec.size();i++){
					zraw=atof(linevec[i].c_str());
					if(zraw==heightmapdataobj.nodata_value){
						zraw=0;
					}
					zraw/=heightmapdataobj.cellsize;
					zraw*=zstretchfact;
                    heightmapdataobj.mydata.push_back(zraw);
				}
			}
		}
	}
}

//***************************************************************** // END LOADING:


void ofxShadowedTerrain::drawMesh(){
    mesh.draw();
}


//***************************************************************** // Shadowcasting image generator:


void ofxShadowedTerrain::setLightAngles(float _xangle, float _yangle){
	ofVec3f dir;
	dir.set(0,1,0);
	dir.rotate(_xangle, ofVec3f(1,0,0));
	dir.rotate(_yangle, ofVec3f(0,1,0));
	dir.normalize();
	setLightDirection(dir);
}


void ofxShadowedTerrain::setLightDirection(ofVec3f _dir){
    // rounding because we dont have sooo many precalculated pictures available
    _dir*=100.0;
	_dir.x=round(_dir.x);
	_dir.y=round(_dir.y);
	_dir.z=round(_dir.z);
	_dir/=100.0;
	
	if(_dir==ofVec3f(0,1,0))return;
	if(_dir==lightdir)return;
	
	lightdir=_dir;
	int accuracy=1;
	
    string filename="media/shadowmaps/sm"+ofToString(_dir.x, accuracy)+"_"+ofToString(_dir.y, accuracy)+"_"+ofToString(_dir.z, accuracy)+".png";
	
	// Check if File already exists, if so return this file
	ifstream ifile(ofToDataPath(filename, true).c_str());
	
    if(ifile){
		shadowimg.loadImage(filename);
		return;
	}
    
	float lightdirfloat[3];
	lightdirfloat[0]=lightdir.x;
	lightdirfloat[1]=lightdir.y;
	lightdirfloat[2]=lightdir.z;
	
	for(int i=0;i<numvalues;i++){
		lightmapforcalc[i]=255;
		int ind=i*4;
		
		// Lichtfarbe
		lightmap[ind]=lightcolor.r;
		lightmap[ind+1]=lightcolor.g;
		lightmap[ind+2]=lightcolor.b;
		lightmap[ind+3]=lightcolor.a;
	}
	
	updateLightMap(heightmap, lightmapforcalc, gridw, lightdirfloat);
	
	for(int i=0;i<numvalues;i++){
		float bright=lightmapforcalc[i];
		int ind=i*4;
        lightmap[ind]=shadowcolor.r;
		lightmap[ind+1]=shadowcolor.g;
		lightmap[ind+2]=shadowcolor.b;
		lightmap[ind+3]=bright*2;
		
	}
	blurimg.setFromPixels(lightmapforcalc, gridw, gridh);
    blurimg.blurGaussian(3);
    shadowimg.setFromPixels(blurimg.getPixels(), gridw, gridh, OF_IMAGE_GRAYSCALE);
    shadowimg.saveImage(filename);
}



ofImage* ofxShadowedTerrain::getShadowImage(){
	return &shadowimg;
}


void ofxShadowedTerrain::prepareForShadows(){
	numvalues=gridw*gridh;
	heightmap=new float[numvalues];
	lightmapforcalc=new unsigned char[numvalues];
	lightmap=new unsigned char[numvalues*4];
	for(int i=0;i<numvalues;i++){
		heightmap[i]=heightmapdataobj.mydata[i];
	}
    heightmaploaded=true;
    shadowimg.allocate(gridw,gridh, OF_IMAGE_GRAYSCALE);
}

// NOT REALLY NEEDED ANYMORE!?
ofVec3f ofxShadowedTerrain::getNormalAt(int x, int y){
	if(x>=0 && x<gridw && y>=0 && y<gridh){
		int index=x+(y*gridw);
		return normals[index];
	}
	return ofVec3f(0,0,0);
}

// NOT REALLY NEEDED ANYMORE!?
ofVec3f ofxShadowedTerrain::getSmoothNormalAt(int x, int y){
	if(x>=0 && x<gridw && y>=0 && y<gridh){
		int index=x+(y*gridw);
		return smoothnormals[index];
	}
	return ofVec3f(0,0,0);
}

// NOT REALLY NEEDED ANYMORE!?
float ofxShadowedTerrain::getHeightAt(int x, int y){
    if(x>=0 && x<gridw && y>=0 && y<gridh){
        int index=x+(y*gridw);
        return heightmapdataobj.mydata[index];;
    }
    return 0;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  LIGHT MAPPING

void ofxShadowedTerrain::updateLightMap(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]){
	
    const int anzvalues=size*size;
	
	// create flag buffer to indicate where we've been
	float *flagMap = new float[size*size];
	for(int i = 0; i < size*size; i++) 
		flagMap[i] = 0;
	
	
	// calculate absolute values for light direction
	float lightDirXMagnitude = lightDir[0];
	float lightDirZMagnitude = lightDir[2];
	if(lightDirXMagnitude < 0) lightDirXMagnitude *= -1;
	if(lightDirZMagnitude < 0) lightDirZMagnitude *= -1;
	
	float distanceStep = sqrtf(lightDir[0] * lightDir[0] + lightDir[1] * lightDir[1]);
	
	// decide which loop will come first, the y loop or x loop
	// based on direction of light, makes calculations faster
	
	// outer loop
	//
#pragma omp parallel for firstprivate(distanceStep, lightDirXMagnitude, \
lightDirZMagnitude, heightmap, lightmap, flagMap, size, lightDir) schedule(static,4)
	for(int y=0; y<size; y++)
	{
		int *X, *Y;
		int iX, iY;
		int dirX, dirY;
		
		// this might seem like a waste, why calculate it here? you can calculate it before... 
		// that's because threading is really picky about sharing variables. the less you share, 
		// the faster it goes.
		if(lightDirXMagnitude > lightDirZMagnitude)
		{
			Y = &iX;
			X = &iY;
			
			if(lightDir[0] < 0)
				dirY = -1;
			else
				dirY = 1;
			
			if(lightDir[2] < 0)
				dirX = -1;
			else
				dirX = 1;
		}
		else
		{
			Y = &iY;
			X = &iX;
			
			if(lightDir[0] < 0)
				dirX = -1;
			else
				dirX = 1;
			
			if(lightDir[2] < 0)
				dirY = -1;
			else
				dirY = 1;
		}
		// if you decide to just do it single-threaded, 
		// just copy the previous block back just above the for loop 
		
		if(dirY < 0)
			iY = size - y - 1;
		else
			iY = y;
		
		// inner loop
		for(int x=0; x<size; x++)
		{
			if(dirX < 0)
				iX = size - x - 1;
			else
				iX = x;
			
			float px, py, height, distance, origX, origY;
			int index;
			
			// travel along the terrain until we:
			// (1) intersect another point
			// (2) find another point with previous collision data
			// (3) or reach the edge of the map
			px = *X;
			py = *Y;
			origX = px;
			origY = py;
			index = (*Y) * size + (*X);
			distance = 0.0f;
			
			// travel along ray
			while(1)
			{
				px -= lightDir[0];
				py -= lightDir[2];
				
				// check if we've reached the boundary
				if(px < 0 || px >= size-1 || py < 0 || py >= size-1)
				{  
					flagMap[index] = -1;
					break;
				}
				
				// calculate interpolated values
				int x0, x1, y0, y1;
				float du, dv;
				float interpolatedHeight, interpolatedFlagMap;
				//float heights[4];
				//float pixels[4];
				float invdu, invdv;
				float w0, w1, w2, w3;
				
				x0 = int(px);
				//x1 = ceilf(px);
				y0 = int(py);
				//y1 = ceilf(py);
				
				x1 = x0 + 1;
				y1 = y0 + 1;
				
				du = px - x0;
				dv = py - y0;
				
				invdu = 1.0 - du;
				invdv = 1.0 - dv;
				w0 = invdu * invdv;
				w1 = invdu * dv;
				w2 = du * invdv;
				w3 = du * dv;
                
                
				
				// compute interpolated flagmap value from point directly below ray
				int ind1,ind2, ind3, ind4;
				ind1=y0*size+x0;
				ind2=y1*size+x0;
				ind3=y0*size+x1;
				ind4=y1*size+x1;
				
				if(ind1<0)ind1=0;
				if(ind1>=anzvalues)ind1=anzvalues-1;
				if(ind2<0)ind2=0;
				if(ind2>=anzvalues)ind2=anzvalues-1;
				if(ind3<0)ind3=0;
				if(ind3>=anzvalues)ind3=anzvalues-1;
				if(ind4<0)ind4=0;
				if(ind4>=anzvalues)ind4=anzvalues-1;
				
				interpolatedHeight = w0*heightmap[ind1] + w1*heightmap[ind2] + w2*heightmap[ind3] + w3*heightmap[ind4];
				interpolatedFlagMap = w0*flagMap[ind1] + w1*flagMap[ind2] + w2*flagMap[ind3] + w3*flagMap[ind4];
				
				
				// get distance from original point to current point
				//distance = sqrtf( (px-origX)*(px-origX) + (py-origY)*(py-origY) );
				distance += distanceStep;
				
				// get height at current point while traveling along light ray
				height = heightmap[index] + lightDir[1]*distance;
				
				// check intersection with either terrain or flagMap
				// if interpolatedHeight is less than interpolatedFlagMap that means 
				// we need to use the flagMap value instead
				// else use the height value
				float val;
				if(interpolatedHeight < interpolatedFlagMap) val = interpolatedFlagMap;
				else val = interpolatedHeight;
				
				if(height < val)
				{
					flagMap[index] = val - height;
					lightmap[index] = distance;
					break;
				}
				
				// check if pixel we've moved to is unshadowed
				// since the flagMap value we're using is interpolated, we will be in 
				// between shadowed and unshadowed areas
				// to compensate for this, simply define some epsilon value and use 
				// this as an offset from -1 to decide
				// if current point under the ray is unshadowed
				static float epsilon = 0.5f;
				if(interpolatedFlagMap < -1.0f+epsilon && interpolatedFlagMap > -1.0f-epsilon)
				{
					flagMap[index] = -1.0f;
					break;
				}   
			}
		}
	}
	
	delete [] flagMap;
}