#pragma once

#include "ofMain.h"
#define WML_INSTANTIATE_BEFORE
#include "RigidMeshDeformer2D.h"

class ofxPuppet {
protected:
	ofMesh origMesh, deformedMesh;
	rmsmesh::RigidMeshDeformer2D deformer;
	bool bConstraintsValid;
	set<unsigned int> vSelected;
	int nSelected;
	
public: 
	void setup(ofMesh & mesh);
	
	void update();
	void draw();
	void drawWireframe();
	 
	void setVertex(int i, const ofVec2f& position);
	void removeVertex(int i);
	ofMesh& getDeformedMesh();
	
	void InitializeDeformedMesh();
	void UpdateDeformedMesh();
	void InvalidateConstraints();
	void ValidateConstraints();
	unsigned int FindHitVertex( float nX, float nY );
};


