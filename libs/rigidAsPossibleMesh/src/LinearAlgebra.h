#ifndef __RMS_LINEAR_ALGEBRA_HH__
#define __RMS_LINEAR_ALGEBRA_HH__

#include <stdio.h>
#include <string>
#include <iostream>
#include <strstream>
#include <cmath>
using namespace std;


// numerical error tolerances
#define RMSUTIL_TOLERANCE 0.000001

namespace rmsmesh {

/* vector and point types. They are not the same, so they do
   not share a data structure. That is all. */
struct Vector {
	float v[4];

	inline void zero() 
		{ v[0] = v[1] = v[2] = v[3] = 0.0f; }

	inline void init(float x, float y, float z) 
		{ v[0] = x; v[1] = y; v[2] = z; v[3] = 0.0f;}

	inline void init(float const * from)
		{ v[0] = from[0]; v[1] = from[1]; v[2] = from[2], v[3] = 0.0f; }

	inline void negate() 
		{  v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2]; }


	inline void operator-=(float const * sub)
		{ v[0] -= sub[0]; v[1] -= sub[1]; v[2] -= sub[2];}

	inline void operator+=(float const * add)
		{ v[0] += add[0]; v[1] += add[1]; v[2] += add[2]; }

	inline void operator*=(float c)
		{ v[0] *= c; v[1] *= c; v[2] *= c; }


	inline float magnitude() const
		{ return (float)sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

	inline float magnitude2() const
		{ return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

	inline float dot(float const * with) const
		{ return (v[0]*with[0] + v[1]*with[1] + v[2]*with[2]); }

	inline float dot4(float const * with) const
		{ return (v[0]*with[0] + v[1]*with[1] + v[2]*with[2] + v[3]*with[3]); }

	inline Vector cross(float const * with) const
		{ Vector t; t.init( v[1]*with[2] - v[2]*with[1],
								-(v[0]*with[2] - v[2]*with[0]),
								v[0]*with[1] - v[1]*with[0]); return t; }

	inline void normalize()
		{ float f = this->magnitude(); v[0] /= f; v[1] /= f; v[2] /= f; }

	inline string toString() const
		{ char buf[256];
		  sprintf(buf, "[ %5.5f %5.5f %5.5f %5.5f]", v[0], v[1], v[2], v[3]);
		  return string(buf); }

	inline int equal(float const * with) const
		{ return (v[0]==with[0] && v[1]==with[1] && v[2]==with[2] ); }

	inline int equivalent(float const * with) const
		{ return (fabs(v[0]-with[0]) < RMSUTIL_TOLERANCE &&
					fabs(v[1]-with[1]) < RMSUTIL_TOLERANCE &&
					fabs(v[2]-with[2]) < RMSUTIL_TOLERANCE ); }

	inline int iszero() const
		{ return (fabs(v[0]) < RMSUTIL_TOLERANCE && 
					fabs(v[1]) < RMSUTIL_TOLERANCE &&
					fabs(v[2]) < RMSUTIL_TOLERANCE); }

	inline void projectVectorOnto( Vector const & v2)
		{ float dp = this->dot(v2.v); dp /= magnitude2(); (*this) *= dp; }

	inline float projectVectorOntoDistance( Vector const & v2)
		{ return (this->dot(v2.v) / magnitude2()); }

	inline float calculateDistanceAlongVector( Vector const & v2)
		{ 	Vector tmp; 
			tmp.init(this->v); 
			tmp.normalize(); 
			tmp.projectVectorOnto(v2); 
			return (tmp.magnitude() * ((this->projectVectorOntoDistance(v2) > 0.0f) ? 1.0f : -1.0f)); }

    // access operators
//	inline float & operator[] (unsigned int n) { return v[n]; }
//	inline const float & operator[] (unsigned int n) const { return v[n]; }

	// casting operators
	inline operator float *() { return v; }
	inline operator const float *() const { return v; }

};



struct Point {
	float p[4];

	inline void zero() 
		{ p[0] = p[1] = p[2] = 0.0f; p[3] = 1.0f; }

	inline void init(float x, float y, float z) 
		{ p[0] = x; p[1] = y; p[2] = z; p[3] = 1.0f; }

	inline void init(float const * from) 
		{ p[0] = from[0]; p[1] = from[1]; p[2] = from[2]; p[3] = 1.0f; }

	inline void negate()
		{  p[0] = -p[0]; p[1] = -p[1]; p[2] = -p[2]; }

	inline void operator-=(float const * sub)
		{ p[0] -= sub[0]; p[1] -= sub[1]; p[2] -= sub[2]; }

	inline void operator+=(float const * add)
		{ p[0] += add[0]; p[1] += add[1]; p[2] += add[2]; }

	inline void operator*=(float c)
		{ p[0] *= c; p[1] *= c; p[2] *= c; }

	inline float magnitude() const
		{ return (float)sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); }

	inline float magnitude2() const
		{ return (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); }

	inline float dot(float const * with) const
		{ return (p[0]*with[0] + p[1]*with[1] + p[2]*with[2]); }

	inline float dot4(float const * with) const
		{ return (p[0]*with[0] + p[1]*with[1] + p[2]*with[2] + p[3]*with[3]); }

	inline Vector cross(float const * with) const
		{ Vector t; t.init( p[1]*with[2] - p[2]*with[1],
							  -(p[0]*with[2] - p[2]*with[0]),
							  p[0]*with[1] - p[1]*with[0]); return t; }

	inline void normalize()
		{ float f = this->magnitude(); p[0] /= f; p[1] /= f; p[2] /= f; }

	inline string toString() const
		{ char buf[256];
		  sprintf(buf, "[ %5.5f %5.5f %5.5f %5.5f]", p[0], p[1], p[2], p[3]);
		  return string(buf); }

	inline int equal(float const * with) const
		{ return( p[0]==with[0] && p[1]==with[1] && p[2]==with[2] ); }

	inline int equivalent(float const * with) const
		{ return( fabs(p[0]-with[0]) < RMSUTIL_TOLERANCE &&
				fabs(p[1]-with[1]) < RMSUTIL_TOLERANCE &&
				fabs(p[2]-with[2]) < RMSUTIL_TOLERANCE ); }

	inline Vector directionTo(float const * to) const
		{ Vector v; 
	      v.init(to[0]-p[0], to[1]-p[1], to[2]-p[2]); 
		  return v; }

	inline float distanceTo(float const * to) const
		{ return (float)sqrt( (to[0]-p[0])*(to[0]-p[0]) 
							 + (to[1]-p[1])*(to[1]-p[1]) 
							 + (to[2]-p[2])*(to[2]-p[2])); } 

    // access operators
//	inline float & operator[] (unsigned int n) { return p[n]; }
//	inline const float & operator[] (unsigned int n) const { return p[n]; }

	// casting operators
	inline operator float *() { return p; }
	inline operator const float *() const { return p; }
};


class RMSLine
{
public:
	RMSLine(Point const & p1, Point const & p2) 
		{ l_p1 = p1; l_p2 = p2; }

	inline float length()
		{ return l_p1.distanceTo(l_p2.p); }

	inline Point const &
	endpoint1() const { return l_p1; }

	inline Point const &
	endpoint2() const { return l_p2; }

protected:
	Point l_p1, l_p2;
};




/* Note: this is _not_ a general matrix class. It is optimized to
   perform 3D linear transformations (scale/translate/rotate) as
   quickly as possible. Some operations
   do not do what you would expect (ie transpose+transpose will
   not give you back the original matrix (!))
*/

class Matrix
{
public:
	Matrix();
	~Matrix();

	Matrix(Matrix const & m2);
	Matrix & operator=(const Matrix & m2);

	void identity();
	void transpose();  /* loses translation (last column) */

	// angles are in radians (!)
	void toScale(float sx, float sy, float sz);
	void toTranslate(float tx, float ty, float tz);
	void toRotateX(float theta);
	void toRotateY(float theta);
	void toRotateZ(float theta);
	void toRotate(float theta, float x, float y, float z);

	void multiply(const Matrix & m2);

	Point & multiply(const Point & v, Point & dest) const;
	Vector & multiply(const Vector & v, Vector & dest) const;
	void multiply(float * vec);

	int invert();     /* may fail */

	string toString() const;

	inline float 
	elem(int r, int c) const { return m_elem[r][c]; }

	inline float const *
	row(int r) const { return m_elem[r]; }

protected:
	float m_elem[3][4];  /* last row is useless, do not store */
	int m_type;          /* types defined in LinearAlgebra.cpp */
};



class Transformation
{
public:
  Transformation();
  ~Transformation();

  bool addMatrix(const Matrix & m);
  bool addMatrixPremultiply(const Matrix & m);

  inline Matrix const * 
  transform() 
	{ return &ml_current_transform; }

  inline Matrix const *
  inverse() 
	{ return &ml_current_inverse; }

  inline void
  transform(float * v) 
	{ ml_current_transform.multiply(v); }

  inline void
  inverse(float * v) 
	{ ml_current_inverse.multiply(v); }

  inline void
  identity() 
	{ ml_current_transform.identity(); ml_current_inverse.identity(); }

protected:
  Matrix ml_current_transform;
  Matrix ml_current_inverse;
};


}   // end namespace rms

#endif






