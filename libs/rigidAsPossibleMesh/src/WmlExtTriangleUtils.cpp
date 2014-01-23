#include <assert.h>
#include <math.h>

#include "WmlExtTriangleUtils.h"

#include <limits>

using namespace Wml;

template <class Real>
void Wml::BarycentricCoords( const Vector3<Real> & vTriVtx1, 
							 const Vector3<Real> & vTriVtx2,
							 const Vector3<Real> & vTriVtx3,
							 const Vector3<Real> & vVertex,
							 Real & fBary1, Real & fBary2, Real & fBary3 )
{

	Wml::Vector3<Real> kV02 = vTriVtx1 - vTriVtx3;
    Wml::Vector3<Real> kV12 = vTriVtx2 - vTriVtx3;
    Wml::Vector3<Real> kPV2 = vVertex - vTriVtx3;

    Real fM00 = kV02.Dot(kV02);
    Real fM01 = kV02.Dot(kV12);
    Real fM11 = kV12.Dot(kV12);
    Real fR0 = kV02.Dot(kPV2);
    Real fR1 = kV12.Dot(kPV2);
    Real fDet = fM00*fM11 - fM01*fM01;
//    ASSERT( Wml::Math<Real>::FAbs(fDet) > (Real)0.0 );
    Real fInvDet = ((Real)1.0)/fDet;

    fBary1 = (fM11*fR0 - fM01*fR1)*fInvDet;
    fBary2 = (fM00*fR1 - fM01*fR0)*fInvDet;
    fBary3 = (Real)1.0 - fBary1 - fBary2;
}

template <class Real>
Real Wml::Area( const Vector3<Real> & vTriVtx1, 
				 const Vector3<Real> & vTriVtx2,
				 const Vector3<Real> & vTriVtx3 )
{
	Wml::Vector3<Real> edge1( vTriVtx2 - vTriVtx1 );
	Wml::Vector3<Real> edge2( vTriVtx3 - vTriVtx1 );
	Wml::Vector3<Real> vCross( edge1.Cross(edge2) );

	return (Real)0.5 * vCross.Length();	
}


template <class Real>
void Wml::BarycentricCoords( const Vector2<Real> & vTriVtx1, 
							 const Vector2<Real> & vTriVtx2,
							 const Vector2<Real> & vTriVtx3,
							 const Vector2<Real> & vVertex,
							 Real & fBary1, Real & fBary2, Real & fBary3 )
{

	Wml::Vector2<Real> kV02 = vTriVtx1 - vTriVtx3;
    Wml::Vector2<Real> kV12 = vTriVtx2 - vTriVtx3;
    Wml::Vector2<Real> kPV2 = vVertex - vTriVtx3;

    Real fM00 = kV02.Dot(kV02);
    Real fM01 = kV02.Dot(kV12);
    Real fM11 = kV12.Dot(kV12);
    Real fR0 = kV02.Dot(kPV2);
    Real fR1 = kV12.Dot(kPV2);
    Real fDet = fM00*fM11 - fM01*fM01;
//    ASSERT( Wml::Math<Real>::FAbs(fDet) > (Real)0.0 );
    Real fInvDet = ((Real)1.0)/fDet;

    fBary1 = (fM11*fR0 - fM01*fR1)*fInvDet;
    fBary2 = (fM00*fR1 - fM01*fR0)*fInvDet;
    fBary3 = (Real)1.0 - fBary1 - fBary2;
}


template <class Real>
Real Wml::Area( const Vector2<Real> & vTriVtx1, 
				 const Vector2<Real> & vTriVtx2,
				 const Vector2<Real> & vTriVtx3 )
{
	Wml::Vector2<Real> edge1( vTriVtx2 - vTriVtx1 );
	Wml::Vector2<Real> edge2( vTriVtx3 - vTriVtx1 );
	//Real fCross( edge1.cross(edge2) );

	
	Real fCross = edge1.X() * edge2.Y() - edge1.Y() * edge2.X();
	
	//m_afTuple[0]*rkV.m_afTuple[1] - m_afTuple[1]*rkV.m_afTuple[0];
	
	
	
	
	return (Real)0.5 * (Real)fabs(fCross);	
}


template <class Real>
Real Wml::Angle( const Vector2<Real> & vTriVtxAngle,
				 const Vector2<Real> & vTriV1,
				 const Vector2<Real> & vTriV2 )
{
	Wml::Vector2<Real> vEdgePrev( vTriV1 - vTriVtxAngle );
	Wml::Vector2<Real> vEdgeNext( vTriV2 - vTriVtxAngle );
	Real fDot = vEdgePrev.Dot(vEdgeNext);
	Real fLen1 = vEdgePrev.Length();
	Real fLen2 = vEdgeNext.Length();
	return (Real)acos( fDot / (fLen1 * fLen2) );
}


template <class Real>
void Wml::Scale( Vector2<Real> & vTriV0,
				 Vector2<Real> & vTriV1,
				 Vector2<Real> & vTriV2,
				 Real fScale )
{
	// find center of mass
	Wml::Vector2<Real> vCentroid( vTriV0 + vTriV1 + vTriV2 );
	vCentroid *= (Real)1.0 / (Real)3.0;

	// convert to vectors, scale and restore
	vTriV0 -= vCentroid;	vTriV0 *= fScale;	vTriV0 += vCentroid;
	vTriV1 -= vCentroid;	vTriV1 *= fScale;	vTriV1 += vCentroid;
	vTriV2 -= vCentroid;	vTriV2 *= fScale;	vTriV2 += vCentroid;
}



template <class Real>
void Wml::StretchMetric1( const Vector3<Real> & q1, 
						 const Vector3<Real> & q2,
						 const Vector3<Real> & q3,
						 const Vector2<Real> & p1,
						 const Vector2<Real> & p2,
						 const Vector2<Real> & p3,
						 Real & MaxSV, Real & MinSV, Real & L2Norm, Real & LInfNorm )
{
	Real s1 = p1.X();
	Real t1 = p1.Y();
	Real s2 = p2.X();
	Real t2 = p2.Y();
	Real s3 = p3.X();
	Real t3 = p3.Y();

	Real A = (Real)0.5 * ( (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if ( A > 0 ) {

		Wml::Vector3<Real> Ss = 
			(q1 * (t2-t3) + q2 * (t3-t1) + q3 * (t1-t2)) / (2*A);
		Wml::Vector3<Real> St = 
			(q1 * (s3-s2) + q2 * (s1-s3) + q3 * (s2-s1)) / (2*A);

		Real a = Ss.Dot(Ss);
		Real b = Ss.Dot(St);
		Real c = St.Dot(St);

		Real discrim = (Real)sqrt( (a-c)*(a-c) + 4*b*b );

		MaxSV = (Real)sqrt( (Real)0.5 * ( (a+c) + discrim ) );
		MinSV = (Real)sqrt( (Real)0.5 * ( (a+c) - discrim ) );

		L2Norm = (Real)sqrt( (Real)0.5 * (a+c)  );
		LInfNorm = MaxSV;
	} else {
		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<Real>::max();
	}

}

template <class Real> 
bool Wml::PtInTri2D( const Vector2<Real> & vTriVtx1,
					 const Vector2<Real> & vTriVtx2,
					 const Vector2<Real> & vTriVtx3,
					 Real fTestX, Real fTestY )
{
	int nSignSum = 0;

	Wml::Vector2<Real> ve = 
		Wml::Vector2<Real>( vTriVtx2[0] - vTriVtx1[0], vTriVtx2[1] - vTriVtx1[1] );
	Wml::Vector2<Real> vp = 
		Wml::Vector2<Real>( fTestX - vTriVtx1[0], fTestY - vTriVtx1[1] );
	Real fSign = ve[0]*vp[1] - ve[1]*vp[0];
	nSignSum += (fSign > 0) ? 1 : -1;

	ve = Wml::Vector2<Real>( vTriVtx3[0] - vTriVtx2[0], vTriVtx3[1] - vTriVtx2[1] );
	vp = Wml::Vector2<Real>( fTestX - vTriVtx2[0], fTestY - vTriVtx2[1] );
	fSign = ve[0]*vp[1] - ve[1]*vp[0];
	nSignSum += (fSign > 0) ? 1 : -1;

	ve = Wml::Vector2<Real>( vTriVtx1[0] - vTriVtx3[0], vTriVtx1[1] - vTriVtx3[1] );
	vp = Wml::Vector2<Real>( fTestX - vTriVtx3[0], fTestY - vTriVtx3[1] );
	fSign = ve[0]*vp[1] - ve[1]*vp[0];
	nSignSum += (fSign > 0) ? 1 : -1;

	return ( abs(nSignSum) == 3 );
}




namespace Wml
{
template  void BarycentricCoords( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3, const Vector3<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template  void BarycentricCoords( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3, const Vector3<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );
template  void BarycentricCoords( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3, const Vector2<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template  void BarycentricCoords( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3, const Vector2<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );


template  float Area( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3 );
template  double Area( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3 );
template  float  Area( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3 );
template  double Area( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3 );

template  float Angle( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3 );
template  double Angle( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3 );
//template  float Angle( const Vector<3,float> & TriVtx1, const Vector<3,float> & TriVtx2,
//											 const Vector<3,float> & TriVtx3 );
//template  double Angle( const Vector<3,double> & TriVtx1, const Vector<3,double> & TriVtx2,
//											 const Vector<3,double> & TriVtx3 );

template  void Scale(  Vector2<float> & TriVtx1,  Vector2<float> & TriVtx2, Vector2<float> & TriVtx3, float fScale );
template  void Scale(  Vector2<double> & TriVtx1,  Vector2<double> & TriVtx2, Vector2<double> & TriVtx3, double fScale );
//template  void Scale(  Vector<3,float> & TriVtx1,  Vector<3,float> & TriVtx2, Vector<3,float> & TriVtx3, float fScale );
//template  void Scale(  Vector<3,double> & TriVtx1,  Vector<3,double> & TriVtx2, Vector<3,double> & TriVtx3, double fScale );

template  void StretchMetric1( const Vector3<float> & q1, const Vector3<float> & q2, const Vector3<float> & q3,
										  const Vector2<float> & p1, const Vector2<float> & p2, const Vector2<float> & p3,
										  float & MaxSV, float & MinSV, float & L2Norm, float & LInfNorm );
template  void StretchMetric1( const Vector3<double> & q1, const Vector3<double> & q2, const Vector3<double> & q3,
										  const Vector2<double> & p1, const Vector2<double> & p2, const Vector2<double> & p3,
										  double & MaxSV, double & MinSV, double & L2Norm, double & LInfNorm );


template  bool PtInTri2D( const Vector2<float> & vTriVtx1, const Vector2<float> & vTriVtx2, const Vector2<float> & vTriVtx3,
									 float fTestX, float fTestY );
template  bool PtInTri2D( const Vector2<double> & vTriVtx1, const Vector2<double> & vTriVtx2, const Vector2<double> & vTriVtx3,
									 double fTestX, double fTestY );
//template  bool PtInTri2D( const Vector<3,float> & vTriVtx1, const Vector<3,float> & vTriVtx2, const Vector<3,float> & vTriVtx3,
//									 float fTestX, float fTestY );
//template  bool PtInTri2D( const Vector<3,double> & vTriVtx1, const Vector<3,double> & vTriVtx2, const Vector<3,double> & vTriVtx3,
//									 double fTestX, double fTestY );

}
