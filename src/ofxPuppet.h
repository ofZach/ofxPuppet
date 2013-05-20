#pragma once

#include "ofMain.h"
#define WML_INSTANTIATE_BEFORE
#include "RigidMeshDeformer2D.h"

class ofxPuppet {
protected:
	ofMesh originalMesh, deformedMesh;
	rmsmesh::RigidMeshDeformer2D deformer;
	bool bConstraintsValid;
	set<unsigned int> controlPoints;
	int nSelected;
	
	void InitializeDeformedMesh();
	void UpdateDeformedMesh();
	void InvalidateConstraints();
	void ValidateConstraints();
	
public: 
	void setup(ofMesh & mesh);
	void update();
	void draw();
	void drawWireframe();
	
	void setControlPoint(int i);
	void setControlPoint(int i, const ofVec2f& position);
	void removeControlPoint(int i);
	
	ofMesh& getOriginalMesh();
	ofMesh& getDeformedMesh();
};


