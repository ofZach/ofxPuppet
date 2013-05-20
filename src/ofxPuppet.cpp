#include "ofxPuppet.h"

#define WML_INSTANTIATE_BEFORE

ofxPuppet::ofxPuppet()
:needsUpdating(false)
,nSelected(0) {
}

void ofxPuppet::setup(ofMesh & mesh){
	mesh.setMode(OF_PRIMITIVE_TRIANGLES);	
	originalMesh = mesh, deformedMesh = mesh;
	deformer.InitializeFromMesh( &originalMesh );
	needsUpdating = true;
}

void ofxPuppet::update(){
	if(needsUpdating) {
		int nConstraints = controlPoints.size();
		set<unsigned int>::iterator cur(controlPoints.begin()), end(controlPoints.end());
		while ( cur != end ) {
			unsigned int nVertex = *cur++;
			ofVec3f vVertex = deformedMesh.getVertices()[nVertex];
			deformer.SetDeformedHandle( nVertex, ofVec2f( vVertex.x, vVertex.y ) );
		}
		deformer.ForceValidation();
		needsUpdating = false;
	}
	deformer.UpdateDeformedMesh( &deformedMesh, true );
}

void ofxPuppet::drawFaces(){	
	deformedMesh.drawFaces();
}

void ofxPuppet::drawWireframe() {
	deformedMesh.drawWireframe();
}

void ofxPuppet::drawControlPoints() {
	ofPushStyle();
	ofNoFill();
	ofSetColor(ofColor::red);
	for(set<unsigned int>::iterator itr = controlPoints.begin(); itr != controlPoints.end(); itr++) {
		ofCircle(deformedMesh.getVertex(*itr), 5); 
	}
	ofPopStyle();
}

void ofxPuppet::setControlPoint(int i) {
	setControlPoint(i, deformedMesh.getVertex(i));
}

void ofxPuppet::setControlPoint(int i, const ofVec2f& position) {
	if (controlPoints.find(i) == controlPoints.end()) {
		controlPoints.insert(i);
	}
	deformedMesh.getVertices()[i].set(position.x, position.y);
	needsUpdating = true; 
}

void ofxPuppet::removeControlPoint(int i) {
	controlPoints.erase(i);
	deformer.RemoveHandle(i);
	needsUpdating = true; 
}

ofMesh& ofxPuppet::getDeformedMesh() {
	return deformedMesh;
}

