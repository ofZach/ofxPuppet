#include "ofApp.h"

void ofApp::setup(){
	ofSetVerticalSync(true);
	puppet.createSquareMesh(ofRectangle(50,50,600,600), 10, 10);
}

void ofApp::update(){
	puppet.update();
}

void ofApp::draw(){
	puppet.draw();
}

void  ofApp::mouseDragged(int x, int y, int button){
	puppet.mouseDragged(x, y, button);
}

void  ofApp::mousePressed(int x, int y, int button){
	puppet.mousePressed(x, y, button);
}

void  ofApp::mouseReleased(int x, int y, int button){
	puppet.mouseReleased(x, y, button);
}












