#include "ofxPuppetInteractive.h"

ofxPuppetInteractive::ofxPuppetInteractive()
:selected(false)
,selectionRadius(10) {
	setEvents(true);
}

// copy constructor
ofxPuppetInteractive::ofxPuppetInteractive(const ofxPuppetInteractive& other)
: selected(other.selected), selectionRadius(other.selectionRadius),
selectedVertex(other.selectedVertex), selectedVertices(other.selectedVertices)
{
    setEvents(true);
}

ofxPuppetInteractive::~ofxPuppetInteractive()
{
    setEvents(false);
}

float getNearestVertex(const ofMesh& mesh, const ofVec2f& target, int& vertexIndex) {
	float bestDistance = 0;
	for(int i = 0; i < mesh.getNumVertices(); i++) {
		float distance = target.distance(mesh.getVertex(i));
		if(distance < bestDistance || i == 0) {
			bestDistance = distance;
			vertexIndex = i;
		}
	}
	return bestDistance;
}

void ofxPuppetInteractive::setSelectionRadius(float selectionRadius) {
	this->selectionRadius = selectionRadius;
}

void ofxPuppetInteractive::setEvents(bool enableEvents) {
	if(enableEvents) {
		ofAddListener(ofEvents().mousePressed, this, &ofxPuppetInteractive::mousePressed);
		ofAddListener(ofEvents().mouseDragged, this, &ofxPuppetInteractive::mouseDragged);
		ofAddListener(ofEvents().mouseReleased, this, &ofxPuppetInteractive::mouseReleased);
	} else {
		ofRemoveListener(ofEvents().mousePressed, this, &ofxPuppetInteractive::mousePressed);
		ofRemoveListener(ofEvents().mouseDragged, this, &ofxPuppetInteractive::mouseDragged);
		ofRemoveListener(ofEvents().mouseReleased, this, &ofxPuppetInteractive::mouseReleased);
	}
}

void ofxPuppetInteractive::mousePressed(ofMouseEventArgs& e){
	float distance = getNearestVertex(deformedMesh, ofVec2f(e.x, e.y), selectedVertex);
	if(distance < selectionRadius) {
		if(e.button == 0) {
			selected = true;
			setControlPoint(selectedVertex);
		} else if(e.button == 2) {
			if(controlPoints.find(selectedVertex) != controlPoints.end()) {
				removeControlPoint(selectedVertex);
			}
		}
	} else {
		selected = false;
	}
}

void ofxPuppetInteractive::mouseDragged(ofMouseEventArgs& e){
	if(selected) {
		setControlPoint(selectedVertex, ofVec2f(e.x, e.y));
	}
}

void ofxPuppetInteractive::mouseReleased(ofMouseEventArgs& e){
	selected = false;
}
