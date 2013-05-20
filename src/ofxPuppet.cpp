#include "ofxPuppet.h"

#define WML_INSTANTIATE_BEFORE

void ofxPuppet::setup(ofMesh & mesh){
	origMesh.clear();
	origMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	deformedMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
	origMesh = mesh;
	deformedMesh = mesh;
	
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
	}
	
	for ( unsigned int i = 0; i < origMesh.getIndices().size(); ++i ) {
		int index = origMesh.getIndices()[i];
		deformedMesh.addIndex(index);
	}
	
	deformer.InitializeFromMesh( &origMesh );
	InvalidateConstraints();
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
	deformedMesh.draw();
	
	std::set<unsigned int>::iterator cur(vSelected.begin()), end(vSelected.end());
	while ( cur != end ) {
		unsigned int nSelected = *cur++;
		ofVec3f vSelected;
		vSelected = deformedMesh.getVertices()[nSelected];  //( nSelected, vSelected );
		ofVec2f vView = ofVec2f(vSelected.x,vSelected.y);
		ofRect(vView.x - 5, vView.y - 5, 10, 10);
	}
}

void ofxPuppet::drawWireframe() {
	deformedMesh.drawWireframe();
}

void ofxPuppet::setVertex(int i, const ofVec2f& position) {
	if (vSelected.find(i) == vSelected.end()) {
		vSelected.insert(i);
	}
	deformedMesh.getVertices()[i].set(position.x, position.y);
	InvalidateConstraints();
}

void ofxPuppet::removeVertex(int i) {
	vSelected.erase(i);
	deformer.RemoveHandle(i);
}

ofMesh& ofxPuppet::getDeformedMesh() {
	return deformedMesh;
}

