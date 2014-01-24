#pragma once

#include "ofxPuppet.h"

class ofxPuppetInteractive : public ofxPuppet {
protected:	
	int selectedVertex;
	bool selected;
	float selectionRadius;
	set<unsigned int> selectedVertices;
	
public:
	ofxPuppetInteractive();
    ofxPuppetInteractive(const ofxPuppetInteractive& other);
	~ofxPuppetInteractive();
	
	void draw();
	
	void setSelectionRadius(float selectionRadius);
	void setEvents(bool enableEvents);
	void mousePressed(ofMouseEventArgs& e);
	void mouseDragged(ofMouseEventArgs& e);
	void mouseReleased(ofMouseEventArgs& e);
};