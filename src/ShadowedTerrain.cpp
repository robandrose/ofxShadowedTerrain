/*
 *  ShadowedTerrain.cpp
 *  LandscapeShadow2
 *
 *  Created by Matthias Rohrbach on 26.08.10.
 *  Copyright 2010 rob & rose grafik. All rights reserved.
 *
 */

#include "ShadowedTerrain.h"
#include "Globals.h"

#include <fstream>
//#include <omp.h>



ShadowedTerrain::ShadowedTerrain(){
	
	terraindata=new TerrainDataProvider();
	
	
	addChild(&container);

}


ShadowedTerrain::~ShadowedTerrain(){
	
}

void ShadowedTerrain::setup(){
	terraindata->loadHeightMap("media/static/terrain_data.txt");
	
}

void ShadowedTerrain::isConfigMode(bool _isconfigmode){
	isconfigmode=_isconfigmode;
	
}

void ShadowedTerrain::firstUpdate(){
	
	ofDisableArbTex();
	
	fbo.allocate(1024, 1024, GL_RGB, 3, true, false, false);
	fbo.clear(1, 1, 1, 0);
	
	ofEnableArbTex();
	
}

void ShadowedTerrain::setShadowAlpha(float _shadowalpha){
	shadowalpha=_shadowalpha;
}

void ShadowedTerrain::setLightAngles(float _xangle, float _yangle){
	if(_xangle<-99)_xangle=-99;
	if(_xangle>99)_xangle=99;
	
	terraindata->setLightAngles(_xangle, _yangle);
	
}

void ShadowedTerrain::update(){

	float scalefactx=(float)fbo.getWidth()/terraindata->getHeightMapWidth();
	float scalefacty=(float)fbo.getHeight()/terraindata->getHeightMapHeight();
	
	fbo.swapIn();
	fbo.setupScreenForMe();
	fbo.clear(1, 1, 1, 1);
	
	ofScale(scalefactx,-scalefacty,1);
	ofTranslate(0, -terraindata->getHeightMapHeight(), 0);
	ofSetColor(255, 255, 255, 255);
	fbocontainer.draw();
	
	ofSetColor(255, 255, 255,shadowalpha);	
	terraindata->getShadowImage()->draw(0,0,1024,1024);

	fbo.setupScreenForThem();
	fbo.swapOut();
	
}

void ShadowedTerrain::_draw(){
	
	if(ofGetFrameNum()==1){
		drawToList();
	}
	
	if(!globals.isdiagnostic){
		fbo.bind();
		mylist.draw();
		fbo.unbind();
	}else{
		disableLighting();
		diagnosticList.draw();
	}
	
	
}

void ShadowedTerrain::drawToList(){
	int ind1, ind2;
	float x1,y1,z1, x2,y2,z2;
	float c1,c2;
	
	mylist.begin();	
	

	terraindata->drawForTexture();
	mylist.end();
	
	diagnosticList.begin();	
	
	terraindata->drawForColor();
	diagnosticList.end();
}

