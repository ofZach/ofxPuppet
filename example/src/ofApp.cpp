#include "ofApp.h"

ofMesh makeGrid(ofRectangle square, int nHoriz, int nVert) {
	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_TRIANGLES);
	for (int j = 0; j < nVert; j++){
    for (int i = 0; i < nHoriz; i++){
			float x = ofMap(i, 0, nHoriz-1, square.x, square.x + square.width);
			float y = ofMap(j, 0, nVert-1, square.y, square.y + square.height);
			mesh.addVertex(ofPoint(x,y));
		}
	}
	for ( unsigned int y = 0; y < (nVert-1); y++ ) {
		for ( unsigned int x = 0; x < nHoriz-1; x++ ) {
			unsigned int nRow1 = y * nHoriz;
			unsigned int nRow2 = (y+1) * nHoriz;
			mesh.addIndex(nRow1 + x);
			mesh.addIndex(nRow2 + x + 1);
			mesh.addIndex(nRow1 + x + 1);
			mesh.addIndex(nRow1 + x);
			mesh.addIndex(nRow2 + x);
			mesh.addIndex(nRow2 + x + 1);
		}
	}
	return mesh;
}

void ofApp::setup(){
	ofSetVerticalSync(true);
	ofMesh mesh = makeGrid(ofRectangle(100,100,600,600), 10, 10);
	puppet.setup(mesh);
	
	puppet.setControlPoint(0); // pin the top left
	puppet.setControlPoint(9); // pin the top right
}

void ofApp::update(){
	puppet.update();
}

void ofApp::draw(){
	ofBackground(0);
	puppet.drawWireframe();
	puppet.drawControlPoints();
}












