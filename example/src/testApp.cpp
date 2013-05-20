#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
	ofSetVerticalSync(true);
	
	puppet.createSquareMesh(ofRectangle(50,50,800,800), 10, 10);
}

//--------------------------------------------------------------
void testApp::update(){
	
	puppet.update();
}

//--------------------------------------------------------------
void testApp::draw(){
	
	
	//puppet.deformedMesh.drawWireframe();
	
	puppet.draw();
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
	
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){
	
}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){
	
}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
	
	puppet.mouseDragged(x, y, button);
	
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
	
	puppet.mousePressed(x, y, button);
	
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
	
	puppet.mouseReleased(x, y, button);;
	
}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){
	
}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){
	
}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 
	
}
















