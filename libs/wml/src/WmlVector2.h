// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3VECTOR2_H
#define WM3VECTOR2_H

#include <assert.h>
//#include "Wm3FoundationLIB.h"
#include "WmlMath.h"

#define WM3_ITEM 

namespace Wml
{
	
	template <class Real>
	class Vector2
	{
	public:
		// construction
		Vector2 ();  // uninitialized
		Vector2 (Real fX, Real fY);
		Vector2 (const Real* afTuple);
		Vector2 (const Vector2& rkV);
		
		// coordinate access
		operator const Real* () const;
		operator Real* ();
		Real operator[] (int i) const;
		Real& operator[] (int i);
		Real X () const;
		Real& X ();
		Real Y () const;
		Real& Y ();
		
		// assignment
		Vector2& operator= (const Vector2& rkV);
		
		// comparison
		bool operator== (const Vector2& rkV) const;
		bool operator!= (const Vector2& rkV) const;
		bool operator<  (const Vector2& rkV) const;
		bool operator<= (const Vector2& rkV) const;
		bool operator>  (const Vector2& rkV) const;
		bool operator>= (const Vector2& rkV) const;
		
		// arithmetic operations
		Vector2 operator+ (const Vector2& rkV) const;
		Vector2 operator- (const Vector2& rkV) const;
		Vector2 operator* (Real fScalar) const;
		Vector2 operator/ (Real fScalar) const;
		Vector2 operator- () const;
		
		// arithmetic updates
		Vector2& operator+= (const Vector2& rkV);
		Vector2& operator-= (const Vector2& rkV);
		Vector2& operator*= (Real fScalar);
		Vector2& operator/= (Real fScalar);
		
		// vector operations
		Real Length () const;
		Real SquaredLength () const;
		Real Dot (const Vector2& rkV) const;
		Real Normalize ();
		
		// returns (y,-x)
		Vector2 Perp () const;
		
		// returns (y,-x)/sqrt(x*x+y*y)
		Vector2 UnitPerp () const;
		
		// returns DotPerp((x,y),(V.x,V.y)) = x*V.y - y*V.x
		Real DotPerp (const Vector2& rkV) const;
		
		// Compute the barycentric coordinates of the point with respect to the
		// triangle <V0,V1,V2>, P = b0*V0 + b1*V1 + b2*V2, where b0 + b1 + b2 = 1.
		void GetBarycentrics (const Vector2<Real>& rkV0,
							  const Vector2<Real>& rkV1, const Vector2<Real>& rkV2,
							  Real afBary[3]) const;
		
		// Gram-Schmidt orthonormalization.  Take linearly independent vectors U
		// and V and compute an orthonormal set (unit length, mutually
		// perpendicular).
		static void Orthonormalize (Vector2& rkU, Vector2& rkV);
		
		// Input V must be initialized to a nonzero vector, output is {U,V}, an
		// orthonormal basis.  A hint is provided about whether or not V is
		// already unit length.
		static void GenerateOrthonormalBasis (Vector2& rkU, Vector2& rkV,
											  bool bUnitLengthV);
		
		// Compute the extreme values.
		static void ComputeExtremes (int iVQuantity, const Vector2* akPoint,
									 Vector2& rkMin, Vector2& rkMax);
		
		// special vectors
		WM3_ITEM static const Vector2 ZERO;    // (0,0)
		WM3_ITEM static const Vector2 UNIT_X;  // (1,0)
		WM3_ITEM static const Vector2 UNIT_Y;  // (0,1)
		WM3_ITEM static const Vector2 ONE;     // (1,1)
		
	private:
		// support for comparisons
		int CompareArrays (const Vector2& rkV) const;
		
		Real m_afTuple[2];
	};
	
	// arithmetic operations
	template <class Real>
	Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV);
	
	// debugging output
//	template <class Real>
//	std::ostream& operator<< (std::ostream& rkOStr, const Vector2<Real>& rkV);
//	
	
#include "Wm3Vector2.inl"
	
	typedef Vector2<float> Vector2f;
	typedef Vector2<double> Vector2d;
	
}

#endif
