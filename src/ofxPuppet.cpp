#include "ofxPuppet.h"

#define WML_INSTANTIATE_BEFORE
#include "RigidMeshDeformer2D.h"

void ofxPuppet::setOrigMesh(ofMesh & mesh){
	origMesh.clear();
	origMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	deformedMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
	origMesh = mesh;
	deformedMesh = mesh;
	
	InitializeDeformedMesh();
}


void ofxPuppet::createSquareMesh( ofRectangle square, int nHoriz, int nVert){
	origMesh.clear();
	origMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	deformedMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
	
	for (int j = 0; j < nVert; j++){
    for (int i = 0; i < nHoriz; i++){
			
			float x = ofMap(i, 0, nHoriz-1, square.x, square.x + square.width);
			float y = ofMap(j, 0, nVert-1, square.y, square.y + square.height);
			origMesh.addVertex(ofPoint(x,y));
		}
	}
	
	for ( unsigned int y = 0; y < (nVert-1); y++ ) {
		unsigned int nRow1 = y * nHoriz;
		unsigned int nRow2 = (y+1) * nHoriz;
		for ( unsigned int x = 0; x < nHoriz-1; x++ ) {
			
			origMesh.addIndex(nRow1 + x);
			origMesh.addIndex(nRow2 + x + 1);
			origMesh.addIndex(nRow1 + x + 1);
			
			origMesh.addIndex(nRow1 + x);
			origMesh.addIndex(nRow2 + x);
			origMesh.addIndex(nRow2 + x + 1);
			
		}
	}
	
	deformedMesh = origMesh;
	
	InitializeDeformedMesh();    
	
}

void ofxPuppet::InitializeDeformedMesh()
{
	deformedMesh.clear();
	
	unsigned int nVerts = origMesh.getVertices().size();
	
	for ( unsigned int i = 0; i < nVerts; ++i ) {
		ofVec2f vVertex(0,0);
		vVertex = origMesh.getVertices()[i];
		deformedMesh.addVertex(vVertex);
		//origMesh.GetVertex(i, vVertex);
		//deformedMesh.AppendVertex(vVertex);
	}
	
	for ( unsigned int i = 0; i < origMesh.getIndices().size(); ++i ) {
		int index = origMesh.getIndices()[i];
		deformedMesh.addIndex(index);
		//origMesh.GetVertex(i, vVertex);
		//deformedMesh.AppendVertex(vVertex);
	}
	
	
	
	//    // TODO this is problematic
	//	unsigned int nTris = origMesh.GetNumTriangles();
	//	for ( unsigned int i = 0; i < nTris; ++i ) {
	//		unsigned int nTriangle[3];
	//		origMesh.GetTriangle(i,nTriangle);
	//		deformedMesh.AppendTriangleData(nTriangle);
	//	}
	
	deformer.InitializeFromMesh( &origMesh );
	InvalidateConstraints();
	
	//	vector < ofVec2f > vert = origMesh.getVertices();
	//	//vector <unsigned int > tris = origMesh.GetTriangles();
	//	
	//	vert = deformedMesh.getVertices();
	//	//tris = deformedMesh.GetTriangles();
	//	
	//	//cout << " --------- deform ------------ " << endl;
	//    //	for (int i = 0; i < vert.size(); i++){
	//    //		cout << vert[i] << ",";	
	//    //	}
	//    //	cout << endl;
	//    //	
	//    //	for (int i = 0; i < tris.size(); i++){
	//    //		cout << tris[i] << ",";	
	//    //	}
	//    //	cout << endl;
	
	
}



void ofxPuppet::UpdateDeformedMesh() 
{
	ValidateConstraints();
	deformer.UpdateDeformedMesh( &deformedMesh, true );
	vector < ofVec3f > vert = deformedMesh.getVertices();
	
	
}


// deformer stuff
void ofxPuppet::InvalidateConstraints() 
{ 
	bConstraintsValid = false; 
}

void ofxPuppet::ValidateConstraints()
{
	if ( bConstraintsValid )
		return;
	
	size_t nConstraints = vSelected.size();
	std::set<unsigned int>::iterator cur(vSelected.begin()), end(vSelected.end());
	while ( cur != end ) {
		unsigned int nVertex = *cur++;
		ofVec3f vVertex;
		
		///cout << " validate " << nVertex << endl;
		
		vVertex = deformedMesh.getVertices()[nVertex]; //( nVertex, vVertex);
		deformer.SetDeformedHandle( nVertex, ofVec2f( vVertex.x, vVertex.y ) );
	}
	
	deformer.ForceValidation();
	
	bConstraintsValid = true;
}


// selection stuff

unsigned int ofxPuppet::FindHitVertex( float nX, float nY ){
	
	
	vector < ofVec3f > verts = deformedMesh.getVertices();
	
	unsigned int nVerts = deformedMesh.getVertices().size();
	for ( unsigned int i = 0; i < nVerts; ++i ) {
		
		ofVec3f vVertex;
		//deformedMesh.GetVertex(i, vVertex);
		vVertex = verts[i];
		
		ofVec2f temp = ofVec2f(vVertex.x, vVertex.y);
		
		ofVec2f vView = temp;; //WorldToView( temp );
		//cout << vView.X() << endl;
		
		float fX = vView.x;
		float fY = vView.y;
		
		//
		
		double fDist = sqrt((double)((nX - fX)*(nX - fX) + (nY-fY)*(nY-fY) ));
		if ( fDist < 5 ) 
			return i;
	}
	
	return std::numeric_limits<unsigned int>::max();
}


void ofxPuppet::update(){
	UpdateDeformedMesh();
}

void ofxPuppet::draw(){
	ofSetColor(255, 200);
	deformedMesh.drawWireframe();
	
	std::set<unsigned int>::iterator cur(vSelected.begin()), end(vSelected.end());
	while ( cur != end ) {
		unsigned int nSelected = *cur++;
		ofVec3f vSelected;
		vSelected = deformedMesh.getVertices()[nSelected];  //( nSelected, vSelected );
		ofVec2f vView = ofVec2f(vSelected.x,vSelected.y);
		ofRect(vView.x - 5, vView.y - 5, 10, 10);
	}
}

void ofxPuppet::setVertex(int i, const ofVec2f& position) {
	if (vSelected.find(i) == vSelected.end()) {
		vSelected.insert(i);
	}
	deformedMesh.getVertices()[i].set(position.x, position.y);
	InvalidateConstraints();
}

//--------------------------------------------------------------
void  ofxPuppet::mouseDragged(int x, int y, int button){
	if ( nSelected != std::numeric_limits<unsigned int>::max() ) {
		setVertex(nSelected, ofVec2f(x, y));
	}
}

//--------------------------------------------------------------
void  ofxPuppet::mousePressed(int x, int y, int button){
	if ( button == 0 ) {
		nSelected = FindHitVertex( (float)x, (float)(y) );
		//cout << nSelected << endl;
	}
	
	if ( button == 2) {
		unsigned int nHit = FindHitVertex( (float)x, (float)(y) );
		if ( nHit != std::numeric_limits<unsigned int>::max() ) {
			if ( vSelected.find(nHit) == vSelected.end() )
				vSelected.insert(nHit);
			else {
				vSelected.erase(nHit);
				deformer.RemoveHandle(nHit);
				
				// restore position
				ofVec3f vVertex;
				vVertex = origMesh.getVertices()[nHit];
				deformedMesh.getVertices()[nHit] = vVertex; //(nHit, vVertex);
			}
			InvalidateConstraints();
			//glutPostRedisplay();
		}
	}
}

//--------------------------------------------------------------
void  ofxPuppet::mouseReleased(int x, int y, int button){
	
	if ( button == 0 ) {
		nSelected = std::numeric_limits<unsigned int>::max();
	} 
	
}
