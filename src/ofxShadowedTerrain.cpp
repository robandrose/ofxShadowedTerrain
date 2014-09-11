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
	
	zstretchfact=1;
    
}

ofxShadowedTerrain::~ofxShadowedTerrain(){
	
}



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
	int accuracy=2;
	
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
	
	updateLightMapFast(heightmap, lightmapforcalc, gridw, lightdirfloat);
	
	for(int i=0;i<numvalues;i++){
		float bright=lightmapforcalc[i];
		int ind=i*4;
		if(bright==0){
			// Schattenfarbe
			lightmap[ind]=shadowcolor.r;
			lightmap[ind+1]=shadowcolor.g;
			lightmap[ind+2]=shadowcolor.b;
			lightmap[ind+3]=shadowcolor.a;
		}
	}
	
	shadowimg.setFromPixels(lightmapforcalc, gridw, gridh, OF_IMAGE_GRAYSCALE);
	shadowimg.saveImage(filename);
}


double ofxShadowedTerrain::value(double x, double y)
{
	if(x>=0 && x<gridw && y>=0 && y<gridh && heightmaploaded){
		int index=(y*gridw)+x;
		return heightmap[index];
	}
	return 0;
}


SPoint ofxShadowedTerrain::lower_bound()
{
	return SPoint(0,0);
}

SPoint ofxShadowedTerrain::upper_bound()
{
	return SPoint(1023,1023);
}




ofImage* ofxShadowedTerrain::getShadowImage(){
	return &shadowimg;
}


//***************************************************************** // LOADING:

void ofxShadowedTerrain::loadMapFromImage(string _filename){
    heightmapimg.loadImage(_filename);
	
    mesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
    int skip = 1;
    
    
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
    

    for(int y = 0; y < gridh - skip; y += skip) {
		for(int x = 0; x < gridw - skip; x += skip) {
			ofVec3f nw = getVertexFromImg(heightmapimg, x, y);
			ofVec3f ne = getVertexFromImg(heightmapimg, x + skip, y);
			ofVec3f sw = getVertexFromImg(heightmapimg, x, y + skip);
			ofVec3f se = getVertexFromImg(heightmapimg, x + skip, y + skip);
			addFace(mesh, nw, ne, se, sw);
		}
	}
  
    mesh.smoothNormals(40);
    prepareForShadows();
    
}



ofVec3f ofxShadowedTerrain::getVertexFromImg(ofFloatImage& img, int x, int y) {
    return ofVec3f(x, y, getValueFromImagePos(img, x, y));
}

float ofxShadowedTerrain::getValueFromImagePos(ofFloatImage& img, int x, int y){
    return 50 * img.getColor(x,y).getBrightness();
}

void ofxShadowedTerrain::addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c) {
	ofVec3f normal = ((b - a).cross(c - a)).normalize();
    mesh.addNormal(normal);
	mesh.addVertex(a);
	mesh.addTexCoord(a);
	mesh.addNormal(normal);
	mesh.addVertex(b);
    mesh.addTexCoord(b);
    mesh.addNormal(normal);
	mesh.addVertex(c);
    mesh.addTexCoord(c);
}

void ofxShadowedTerrain::addFace(ofMesh& mesh, ofVec3f a, ofVec3f b, ofVec3f c, ofVec3f d) {
	addFace(mesh, a, b, c);
	addFace(mesh, a, c, d);
}


// LOAD FROM TEXTFILE:

void ofxShadowedTerrain::loadMapFromTextfile(string _filename){
  
    loadHeightmapData(_filename);
    mesh.setMode(OF_PRIMITIVE_TRIANGLES);
	int skip = 1;
	gridw = heightmapdataobj.ncols;
	gridh = heightmapdataobj.nrows;
    
    for(int y = 0; y < gridh - skip; y += skip) {
		for(int x = 0; x < gridw - skip; x += skip) {
            ofVec3f nw= heightmapdataobj.getVertexAt(x, y);
            ofVec3f ne=heightmapdataobj.getVertexAt(x+skip, y);
			ofVec3f sw=heightmapdataobj.getVertexAt(x, y+skip);
			ofVec3f se=heightmapdataobj.getVertexAt(x+skip, y+skip);
            addFace(mesh, nw, ne, se, sw);
		}
	}
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






void ofxShadowedTerrain::loadMapFromTextfileOld(string _filename){
	
	ifstream data;
	string line;
	float zraw;
	
	data.open(ofToDataPath(_filename,true).c_str());
	
	while(!data.eof()){
		getline(data, line);
		vector<string>linevec=ofSplitString(line, " ");
		if(linevec.size()>0){
			if(linevec.at(0)=="NCOLS"){
				ncols=atoi(linevec[1].c_str());
			}else if(linevec.at(0)=="NROWS"){
				nrows=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="XLLCORNER"){
				xllcorner=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="YLLCORNER"){
				yllcorner=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="CELLSIZE"){
				cellsize=atoi(linevec[1].c_str());
			}
			else if(linevec.at(0)=="NODATA_VALUE"){
				nodata_value=atoi(linevec[1].c_str());
			}
			else{
				for(int i=0;i<linevec.size();i++){
					zraw=atof(linevec[i].c_str());
					if(zraw==nodata_value){
						zraw=0;
					}
					zraw/=cellsize;
					zraw*=zstretchfact;
					heightmapdata.push_back(zraw);
				}
			}
		}
	}
	

	
    
	// calculate Vertices
	
	int offx=0;
	int offy=0;
	
	gridw=ncols;
	gridh=nrows;
	
	numvalues=gridw*gridh;	
	heightmap=new float[numvalues];
	lightmapforcalc=new unsigned char[numvalues];
	lightmap=new unsigned char[numvalues*4];
	
	
	int xind=0;
	int yind=0;
	int origind;
	float zval=0;
	
	
	for(int i=0;i<numvalues;i++){
		xind=(int)i%gridw;
		yind=(int)floor(i/gridw);
		origind=offx+xind+(offy+yind)*ncols;
		heightmap[i]=heightmapdata[origind];
	}
	
	
	// Calculate Normals:		
	
	normals=new ofVec3f[numvalues];
	smoothnormals=new ofVec3f[numvalues];
	
	
	ofVec3f norm;
	norm.set(0,0,0);
	
	for(int i=0;i<numvalues;i++){
		normals[i]=norm;
		smoothnormals[i]=norm;
	}
		
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
	float x4,y4,z4;
	float c1, c2;
	
	int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
	ofVec3f vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8;
	
	for(int h=1;h<gridh-1; h++){
		for(int i=1;i<gridw-1; i++){
	
			ind0=i+h*gridw;
			ind1=(i)+(h-1)*gridw;
			ind2=(i+1)+(h)*gridw;
			ind3=(i)+(h+1)*gridw;
			ind4=(i-1)+(h)*gridw;
			
            norm.x=heightmap[ind4]-heightmap[ind2];
			norm.y=heightmap[ind1]-heightmap[ind3];
			norm.z=2;
			norm.normalize();
            normals[ind0].set(norm);
		}
	}
	
	float smoothfakt=.6;
    
	for(int k=0;k<10;k++){
        for(int h=1;h<gridh-1; h++){
            for(int i=1;i<gridw-1; i++){
                ind0=i+h*gridw;
                
                ind1=(i)+(h-1)*gridw;
                ind2=(i+1)+(h)*gridw;
                ind3=(i)+(h+1)*gridw;
                ind4=(i-1)+(h)*gridw;
                
                norm=normals[ind0];
                norm+=normals[ind1]*smoothfakt;
                norm+=normals[ind2]*smoothfakt;
                norm+=normals[ind3]*smoothfakt;
                norm+=normals[ind4]*smoothfakt;
                
                norm.normalize();
                smoothnormals[ind0]=norm;
                
            }
        }
	}
	
	heightmaploaded=true;
    shadowimg.allocate(gridw,gridh, OF_IMAGE_GRAYSCALE);
}

ofVec3f ofxShadowedTerrain::getNormalAt(int x, int y){
	
	if(x>=0 && x<gridw && y>=0 && y<gridh){
		int index=x+(y*gridw);
		return normals[index];
	}
	return ofVec3f(0,0,0);
}

ofVec3f ofxShadowedTerrain::getSmoothNormalAt(int x, int y){
	if(x>=0 && x<gridw && y>=0 && y<gridh){
		int index=x+(y*gridw);
		return smoothnormals[index];
	}
	return ofVec3f(0,0,0);
}


void ofxShadowedTerrain::drawMesh(){
    mesh.draw();
}

void ofxShadowedTerrain::drawForColor(){
	int ind1, ind2;
	float x1,y1,z1, x2,y2,z2;
	
	ofColor c1;
	ofColor c2;
	
	
	for(int h=0;h<gridh-1; h++){
		
		glBegin(GL_TRIANGLE_STRIP);
		
		for(int i=0;i<gridw; i++){
			
			ind1=i+(h)*gridw;
			ind2=i+(h+1)*gridw;
			
			x1=i;
			y1=h;
			z1=heightmap[ind1];
			
			x2=i;
			y2=h+1;
			z2=heightmap[ind2];
			
			c1.setHsb(z1*10, 1, 1, 1);
			c2.setHsb(z2*10,1,1,1);
			
			
			glColor3f(c1.r,c1.g,c1.b);
			//glTexCoord2d(i/(double)gridw,h/(double)gridh);
			//glNormal3f(normals[ind1].x, normals[ind1].y, normals[ind1].z);
			glVertex3f(x1,y1,z1);
			
			glColor3f(c2.r,c2.g,c2.b);
			//glTexCoord2d((i)/(double)gridw,(h+1)/(double)gridh);
			//glNormal3f(normals[ind2].x, normals[ind2].y, normals[ind2].z);
			glVertex3f(x2, y2, z2);
			
		}
		glEnd();
	}
	
}


void ofxShadowedTerrain::drawForTexture(){
	int ind1, ind2;
	float x1,y1,z1, x2,y2,z2;
	float c1,c2;
	
	int step=3;
	
	for(int h=0;h<gridh-step; h+=step){
		
		glBegin(GL_TRIANGLE_STRIP);
		
		int i;
		
		for(i=0;i<gridw; i+=step){
			
			ind1=i+(h)*gridw;
			ind2=i+(h+step)*gridw;
			
			x1=i;
			y1=h;
			z1=heightmap[ind1];
			
			x2=i;
			y2=h+step;
			z2=heightmap[ind2];
			
			
            glTexCoord2d(i,h);
			glNormal3f(normals[ind1].x, normals[ind1].y, normals[ind1].z);
			glVertex3f(x1,y1,z1);
			
			glTexCoord2d(i,h+step);
			glNormal3f(normals[ind2].x, normals[ind2].y, normals[ind2].z);
			glVertex3f(x2, y2, z2);
			
			if(i==gridw-1)break;
			if(i+step>gridw-1){
				i=gridw-1-step;
			}
			
		}
		glEnd();
		if(h==gridh-step-1)break;
		if(h+step>gridh-step-1){
			h=gridh-1-step-step;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  LIGHT MAPPING


void ofxShadowedTerrain::updateLightMap(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]){
	
	// create flag buffer to indicate where we've been
	float *flagMap = new float[size*size];
	float *shadowColor=new float[3];
	
	shadowColor[0]=0;
	shadowColor[1]=0;
	shadowColor[2]=0;
	
	for(int i = 0; i < size*size; i++) 
		flagMap[i] = 0;
	
	int *X, *Y;
	int iX, iY;
	int dirX, dirY;
	
	// calculate absolute values for light direction
	float lightDirXMagnitude = lightDir[0];
	float lightDirZMagnitude = lightDir[2];
	if(lightDirXMagnitude < 0) lightDirXMagnitude *= -1;
	if(lightDirZMagnitude < 0) lightDirZMagnitude *= -1;
	
	// decide which loop will come first, the y loop or x loop
	// based on direction of light, makes calculations faster
	if(lightDirXMagnitude > lightDirZMagnitude)
	{
		Y = &iX;
		X = &iY;
		
		if(lightDir[0] < 0)
		{
			iY = size-1;
			dirY = -1;
		}
		else
		{
			iY = 0;
			dirY = 1;
		}
		
		if(lightDir[2] < 0)
		{
			iX = size-1;
			dirX = -1;
		}
		else
		{
			iX = 0;
			dirX = 1;
		}
	}
	else
	{
		Y = &iY;
		X = &iX;
		
		if(lightDir[0] < 0)
		{
			iX = size-1;
			dirX = -1;
		}
		else
		{
			iX = 0;
			dirX = 1;
		}
		
		if(lightDir[2] < 0)
		{
			iY = size-1;
			dirY = -1;
		}
		else
		{
			iY = 0;
			dirY = 1;
		}
	}
	
	// outer loop
	while(1)
	{
		// inner loop
		while(1)
		{
			// travel along the terrain until we:
			// (1) intersect another point
			// (2) find another point with previous collision data
			// (3) or reach the edge of the map
			float px = *X;
			float py = *Y;
			int index = (*Y) * size + (*X);
			
			// travel along ray
			while(1)
			{
				px -= lightDir[0];
				py -= lightDir[2];
				
				// check if we've reached the boundary
				if(px < 0 || px >= size || py < 0 || py >= size)
				{  
					flagMap[index] = -1;
					break;
				}
				
				// calculate interpolated values
				static int x0, x1, y0, y1;
				static float du, dv;
				static float interpolatedHeight, interpolatedFlagMap;
				static float heights[4];
				static float pixels[4];
				static float invdu, invdv;
				static float w0, w1, w2, w3;
				
				x0 = floor(px);
				x1 = ceil(px);
				y0 = floor(py);
				y1 = ceil(py);
				
				du = px - x0;
				dv = py - y0;
				
				invdu = 1.0 - du;
				invdv = 1.0 - dv;
				w0 = invdu * invdv;
				w1 = invdu * dv;
				w2 = du * invdv;
				w3 = du * dv;
				
				// compute interpolated height value from the heightmap direction below ray
				heights[0] = heightmap[y0*size+x0];
				heights[1] = heightmap[y1*size+x0];
				heights[2] = heightmap[y0*size+x1];
				heights[3] = heightmap[y1*size+x1];
				interpolatedHeight = w0*heights[0] + w1*heights[1] + w2*heights[2] + w3*heights[3];
				
				// compute interpolated flagmap value from point directly below ray
				pixels[0] = flagMap[y0*size+x0];
				pixels[1] = flagMap[y1*size+x0];
				pixels[2] = flagMap[y0*size+x1];
				pixels[3] = flagMap[y1*size+x1];
				interpolatedFlagMap = w0*pixels[0] + w1*pixels[1] + w2*pixels[2] + w3*pixels[3];
				
				// get distance from original point to current point
				float distance = sqrt( (px-*X)*(px-*X) + (py-*Y)*(py-*Y) );
				
				// get height at current point while traveling along light ray
				float height = heightmap[index] + lightDir[1]*distance;
				
				// check intersection with either terrain or flagMap
				// if interpolatedHeight is less than interpolatedFlagMap that means we need to use the flagMap value instead
				// else use the height value
				static float val;
				val = interpolatedHeight;
				if(interpolatedHeight < interpolatedFlagMap) val = interpolatedFlagMap;
				if(height < val)
				{
					flagMap[index] = val - height;
					
					lightmap[index*3+0] = shadowColor[0];
					lightmap[index*3+1] = shadowColor[1];
					lightmap[index*3+2] = shadowColor[2];
					
					break;
				}
				
				// check if pixel we've moved to is unshadowed
				// since the flagMap value we're using is interpolated, we will be in between shadowed and unshadowed areas
				// to compensate for this, simply define some epsilon value and use this as an offset from -1 to decide
				// if current point under the ray is unshadowed
				static float epsilon = 0.5f;
				if(interpolatedFlagMap < -1.0f+epsilon && interpolatedFlagMap > -1.0f-epsilon)
				{
					flagMap[index] = -1.0f;
					break;
				}   
			}
			
			// update inner loop variable
			if(dirY < 0)
			{
				iY--;
				if(iY < 0) break;
			}
			else
			{
				iY++;
				if(iY >= size) break;
			}
		}
		
		// reset inner loop starting point
		if(dirY < 0) iY = size - 1;
		else iY = 0;
		
		// update outer loop variable
		if(dirX < 0)
		{
			iX--;
			if(iX < 0) break;
		}
		else
		{
			iX++;
			if(iX >= size) break;
		}
	}
	delete [] flagMap;
}





void ofxShadowedTerrain::updateLightMapFast(float *heightmap, unsigned char *lightmap, int size, float lightDir[3]){
	
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
				if(px < 0 || px >= size || py < 0 || py >= size)
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
				
				// compute interpolated height value from the heightmap direction below ray
				interpolatedHeight = w0*heightmap[y0*size+x0] + w1*heightmap[y1*size+x0] + 
				w2*heightmap[y0*size+x1] + w3*heightmap[y1*size+x1];
				
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
					
					lightmap[index] = 0;
					
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



/*



int intersect_map(const vector3& iv,const ray& r,Image* hm,float fHeightScale){
    int w,hits;
    float d,h,D;
    vector3 v,dir;
    
    v = iv + r.direction;
    w = hm->w;
    
    hits = 0;
    
    while (!(( v.x >= w-1 ) || ( v.x <= 0 ) || ( v.z >= w-1 ) || ( v.z <= 0 ))){
        // length of lightdir's projection
        D = Magnitude(vector3(v.x,0,v.z)-vector3(r.origin.x,0,r.origin.z));
        d = Magnitude(iv-v);            // light direction
        h = iv.y + (d*r.origin.y) / D;  // X(P) point
        
        // check if height in point P is bigger than point X's height
        if (hm->data[ifloor(v.z)* w + ifloor(v.x)] * fHeightScale > h){
            hits++;   // if so, mark as hit, and skip this work point.
            break;
        };
        
        dir = r.direction;
        dir.y = 0;
        v += Normalize(dir);   // fetch new working point
    };
    return hits;
};

Image* genLightmap(char* normal,Image* hm,vector3 fSunDir,int w,float fAmbient){
    int i,j,hits;
    float f,dot;
    vector3 n,fVertex;
    Image* lmap;
    ray r;
    
    float fHeightScale = 10.0f / 255.0f;
    lmap = new Image(w,w,1);
    if (!lmap){printf("(!) Error: cannot alloc lightmap!\n");return 0;};
    
    for (j=0; jdata[j*w+i] * fHeightScale;
         fVertex.z = j;
         
         f = fAmbient ;
         
         r.origin = fVertex + fSunDir * 2000.0f;
         r.direction = fSunDir;
         
         // checks current working point for intersection
         if (!intersect_map(fVertex,r,hm,fHeightScale)){
             // compute the lighting equation
             n.x = (float)(normal[3*(j*w+i)+0]);
             n.y = (float)(normal[3*(j*w+i)+1]);
             n.z = (float)(normal[3*(j*w+i)+2]);
             f += 0.5f*(1.0f+DotProduct(Normalize(n),Normalize(fSunDir)));
             if (f>1.0f) f = 1.0f;
         };
         
         dot = f * 255.0f;
         lmap->data[j*w+i] = (unsigned char)dot;
         };
         };
         return lmap;
         };

*/