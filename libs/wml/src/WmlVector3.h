// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3VECTOR3_H
#define WM3VECTOR3_H

//#include "Wm3FoundationLIB.h"
#include "WmlMath.h"
#include <assert.h>

#define WM3_ITEM 

namespace Wml
{
	
	template <class Real>
	class Vector3
	{
	public:
		// construction
		Vector3 ();  // uninitialized
		Vector3 (Real fX, Real fY, Real fZ);
		Vector3 (const Real* afTuple);
		Vector3 (const Vector3& rkV);
		
		// coordinate access
		operator const Real* () const;
		operator Real* ();
		Real operator[] (int i) const;
		Real& operator[] (int i);
		Real X () const;
		Real& X ();
		Real Y () const;
		Real& Y ();
		Real Z () const;
		Real& Z ();
		
		// assignment
		Vector3& operator= (const Vector3& rkV);
		
		// comparison
		bool operator== (const Vector3& rkV) const;
		bool operator!= (const Vector3& rkV) const;
		bool operator<  (const Vector3& rkV) const;
		bool operator<= (const Vector3& rkV) const;
		bool operator>  (const Vector3& rkV) const;
		bool operator>= (const Vector3& rkV) const;
		
		// arithmetic operations
		Vector3 operator+ (const Vector3& rkV) const;
		Vector3 operator- (const Vector3& rkV) const;
		Vector3 operator* (Real fScalar) const;
		Vector3 operator/ (Real fScalar) const;
		Vector3 operator- () const;
		
		// arithmetic updates
		Vector3& operator+= (const Vector3& rkV);
		Vector3& operator-= (const Vector3& rkV);
		Vector3& operator*= (Real fScalar);
		Vector3& operator/= (Real fScalar);
		
		// vector operations
		Real Length () const;
		Real SquaredLength () const;
		Real Dot (const Vector3& rkV) const;
		Real Normalize ();
		
		// The cross products are computed using the right-handed rule.  Be aware
		// that some graphics APIs use a left-handed rule.  If you have to compute
		// a cross product with these functions and send the result to the API
		// that expects left-handed, you will need to change sign on the vector
		// (replace each component value c by -c).
		Vector3 Cross (const Vector3& rkV) const;
		Vector3 UnitCross (const Vector3& rkV) const;
		
		// Compute the barycentric coordinates of the point with respect to the
		// tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
		// b0 + b1 + b2 + b3 = 1.
		void GetBarycentrics (const Vector3<Real>& rkV0,
							  const Vector3<Real>& rkV1, const Vector3<Real>& rkV2,
							  const Vector3<Real>& rkV3, Real afBary[4]) const;
		
		// Gram-Schmidt orthonormalization.  Take linearly independent vectors
		// U, V, and W and compute an orthonormal set (unit length, mutually
		// perpendicular).
		static void Orthonormalize (Vector3& rkU, Vector3& rkV, Vector3& rkW);
		static void Orthonormalize (Vector3* akV);
		
		// Input W must be initialized to a nonzero vector, output is {U,V,W},
		// an orthonormal basis.  A hint is provided about whether or not W
		// is already unit length.
		static void GenerateOrthonormalBasis (Vector3& rkU, Vector3& rkV,
											  Vector3& rkW, bool bUnitLengthW);
		
		// Compute the extreme values.
		static void ComputeExtremes (int iVQuantity, const Vector3* akPoint,
									 Vector3& rkMin, Vector3& rkMax);
		
		// special vectors
		WM3_ITEM static const Vector3 ZERO;    // (0,0,0)
		WM3_ITEM static const Vector3 UNIT_X;  // (1,0,0)
		WM3_ITEM static const Vector3 UNIT_Y;  // (0,1,0)
		WM3_ITEM static const Vector3 UNIT_Z;  // (0,0,1)
		WM3_ITEM static const Vector3 ONE;     // (1,1,1)
		
	private:
		// support for comparisons
		int CompareArrays (const Vector3& rkV) const;
		
		Real m_afTuple[3];
	};
	
	// arithmetic operations
	template <class Real>
	Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV);
	
	// debugging output
	//template <class Real>
	//std::ostream& operator<< (std::ostream& rkOStr, const Vector3<Real>& rkV);
	
#include "Wm3Vector3.inl"
	
	typedef Vector3<float> Vector3f;
	typedef Vector3<double> Vector3d;
	
}

#endif
