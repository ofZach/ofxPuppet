#pragma once

#include "ofMain.h" 
#include <Set>
#define WML_INSTANTIATE_BEFORE
#include "RigidMeshDeformer2D.h"



 class ofxPuppet {
 
     public: 
     
         void setOrigMesh(ofMesh & mesh);
         void createSquareMesh( ofRectangle square, int nHoriz, int nVert);
        
         void update();
         void draw();
     
         ofMesh origMesh;
         ofMesh deformedMesh;
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
 
 
