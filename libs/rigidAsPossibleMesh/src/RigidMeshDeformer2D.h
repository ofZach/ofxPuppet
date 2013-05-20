#ifndef __RMSIMPLICIT_RIGID_MESH_DEFORMER_2D_H__
#define __RMSIMPLICIT_RIGID_MESH_DEFORMER_2D_H__



#define WML_ITEM 
#include <map>
#include <set>
#include <limits>
#include <WmlBox2.h>
#include <WmlGMatrix.h>
#include "WmlLinearSystemExt.h"
#include "ofMain.h"
#include "opencv2/opencv.hpp"

//Vertex ?

namespace rmsmesh {
	
	
	class RigidMeshDeformer2D
	{
	public:
		RigidMeshDeformer2D( );
		~RigidMeshDeformer2D() {};
		
		
		void ForceValidation() { ValidateSetup(); }
		void RemoveHandle( unsigned int nHandle );
		
		/*
		 * interface stuff
		 */
		//unsigned int GetNumHandles();
		
		//const Wml::Vector2f & GetInitialHandle(unsigned int nHandle); 
		//const Wml::Vector2f & GetDeformedHandle( unsigned int nHandle );
		
		//! nHandle is vertex ID
		void SetDeformedHandle( unsigned int nHandle, const ofVec2f & vHandle );
		
		//void TransformPoint( Wml::Vector2f & vTransform );
		void UnTransformPoint( ofVec2f & vTransform );
		
		/*
		 * mesh handling
		 */
		void InitializeFromMesh( ofMesh * pMesh );
		void UpdateDeformedMesh( ofMesh * pMesh, bool bRigid );
		
		
		/*
		 * debug
		 */
	
    protected:
		
		struct Vertex {
			ofVec2f vPosition;
		};
		
	public:
		struct Triangle {
			unsigned int nVerts[3];
			
			// definition of each vertex in triangle-local coordinate system
			ofVec2f vTriCoords[3];
			
			// un-scaled triangle
			ofVec2f vScaled[3];
			
			// pre-computed matrices for triangle scaling step
			Wml::GMatrixd mF, mC;
		};
		
	protected:
        
        
        
		std::vector<Vertex> m_vInitialVerts;
		std::vector<Vertex> m_vDeformedVerts;
		
		std::vector<Triangle> m_vTriangles;
		
		
		struct Constraint {
			unsigned int nVertex;
			ofVec2f vConstrainedPos;
			
			Constraint() { nVertex = 0; vConstrainedPos = ofVec2f(0,0); }
			Constraint( unsigned int nVert, const ofVec2f & vPos ) { nVertex = nVert; vConstrainedPos = vPos; }
			
			bool operator<(const Constraint & c2) const
			{ return nVertex < c2.nVertex; }
		};
		
		std::set<Constraint> m_vConstraints;
		void UpdateConstraint( Constraint & cons );
		
		
		bool m_bSetupValid;
		void InvalidateSetup() { m_bSetupValid = false; }
		void ValidateSetup();
		
		
		Wml::GMatrixd m_mFirstMatrix;
		std::vector<unsigned int> m_vVertexMap;
		Wml::GMatrixd m_mHXPrime, m_mHYPrime;
		Wml::GMatrixd m_mDX, m_mDY;
		
		cv::SVD m_mSVDDecompX, m_mSVDDecompY;
		Wml::LinearSystemExtd::LUData m_mLUDecompX, m_mLUDecompY;
		
		void PrecomputeOrientationMatrix();
		void PrecomputeScalingMatrices( unsigned int nTriangle );
		void PrecomputeFittingMatrices();
		
		void ValidateDeformedMesh( bool bRigid );
		void UpdateScaledTriangle( unsigned int nTriangle );
		void ApplyFittingStep();
		
		ofVec2f GetInitialVert( unsigned int nVert ) 
		{ return ofVec2f( m_vInitialVerts[ nVert ].vPosition.x, m_vInitialVerts[ nVert ].vPosition.y ); }
	};
	
	
} // namespace rmsimplicit


#endif // __RMSIMPLICIT_RIGID_MESH_DEFORMER_2D_H__
