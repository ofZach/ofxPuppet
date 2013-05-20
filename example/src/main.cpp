#include "ofAppGlutWindow.h"
#include "ofApp.h"

int main() {
	ofAppGlutWindow window;
	ofSetupOpenGL(&window, 1024, 768, OF_WINDOW);
	ofRunApp(new ofApp());
}
