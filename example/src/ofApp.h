#pragma once

#include "ofMain.h"
#include "ofxPuppetInteractive.h"

class ofApp : public ofBaseApp{
public:
	void setup();
	void update();
	void draw();
	
	ofxPuppetInteractive puppet;
};
