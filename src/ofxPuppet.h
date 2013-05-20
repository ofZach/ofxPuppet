#pragma once

#include "ofMain.h" 
#include <Set>
#define WML_INSTANTIATE_BEFORE
#include "RigidMeshDeformer2D.h"

class ofxPuppet {
public: 
	
	void setMesh(ofMesh & mesh);
	
	void update();
	void draw();
	void drawWireframe();
	 
	void setVertex(int i, const ofVec2f& position);
	
	ofMesh origMesh, deformedMesh;
	rmsmesh::RigidMeshDeformer2D deformer;
	bool bConstraintsValid;
	std::set<unsigned int> vSelected;
	int nSelected;
	
	void InitializeDeformedMesh();
	void UpdateDeformedMesh();
	void InvalidateConstraints();
	void ValidateConstraints();
	unsigned int FindHitVertex( float nX, float nY );
	
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
};


