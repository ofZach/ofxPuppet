#include "LinearAlgebra.h"

using namespace rmsmesh;

// Matrix Types 
// You cannot change these types to an enum. They are used in bit operations.
#define MT_NONE        0
#define MT_ROTATION    1
#define MT_TRANSLATION 1<<1
#define MT_SCALE       1<<2
#define MT_UNKNOWN     1<<10


/*
 * things to do
 *  - identity type
 *  - optimize routines based on matrix type (??)
 *
 *
 */


Matrix::Matrix()
  : m_type(MT_NONE)
{
  identity();
}
 

Matrix::Matrix(Matrix const & m2)
{
	memcpy(m_elem, m2.m_elem, sizeof(float)*12);
	m_type = m2.m_type;
}

Matrix & 
Matrix::operator=(Matrix const & m2)
{
	memcpy(m_elem, m2.m_elem, sizeof(float)*12);
	m_type = m2.m_type;
	return *this;
}


Matrix::~Matrix()
{
}


void 
Matrix::identity()
{
	memset(m_elem, 0, sizeof(float)*12);
	m_elem[0][0] = 1.0;
	m_elem[1][1] = 1.0;
	m_elem[2][2] = 1.0;
}


void 
Matrix::transpose()
{
	float t;
	t = m_elem[0][1];
	m_elem[0][1] = m_elem[1][0];
	m_elem[1][0] = t;

	t = m_elem[0][2];
	m_elem[0][2] = m_elem[2][0];
	m_elem[2][0] = t;

	t = m_elem[1][2];
	m_elem[1][2] = m_elem[2][1];
	m_elem[2][1] = t;

	m_elem[0][3] = m_elem[1][3] = m_elem[2][3] = 0.0;
}



void 
Matrix::toScale(float sx, float sy, float sz)
{
	m_elem[0][0] = sx;
	m_elem[1][1] = sy;
	m_elem[2][2] = sz;
	m_type |= MT_SCALE;
}

void 
Matrix::toTranslate(float tx, float ty, float tz)
{
	m_elem[0][3] = tx;
	m_elem[1][3] = ty;
	m_elem[2][3] = tz;
	m_type |= MT_TRANSLATION;
}

void 
Matrix::toRotateX(float theta)
{
	m_elem[1][1] = (float)cos(theta);
	m_elem[1][2] = -(float)sin(theta);
	m_elem[2][1] = -m_elem[1][2];
	m_elem[2][2] = m_elem[1][1];
	m_type |= MT_ROTATION;
}

void 
Matrix::
toRotate(float theta, float x, float y, float z)
{
	float c = (float)cos(theta);
	float s = (float)sin(theta);
	float t = 1-c;
	  
	m_elem[0][0] = t*x*x + c;
	m_elem[0][1] = t*x*y + s*z;
	m_elem[0][2] = t*x*z - s*y;
	m_elem[0][3] = 0.0;

	m_elem[1][0] = t*x*y - s*z;
	m_elem[1][1] = t*y*y + c;
	m_elem[1][2] = t*y*z + s*x;
	m_elem[1][3] = 0.0;

	m_elem[2][0] = t*x*z + s*y;
	m_elem[2][1] = t*y*z - s*x;
	m_elem[2][2] = t*z*z + c;
	m_elem[2][3] = 0.0;
}



void 
Matrix::toRotateY(float theta)
{
	m_elem[0][0] = (float)cos(theta);
	m_elem[0][2] = (float)sin(theta);
	m_elem[2][0] = -m_elem[0][2];
	m_elem[2][2] = m_elem[0][0];
	m_type |= MT_ROTATION;
}


void 
Matrix::toRotateZ(float theta)
{
	m_elem[0][0] = (float)cos(theta);
	m_elem[0][1] = -(float)sin(theta);
	m_elem[1][0] = -m_elem[0][1];
	m_elem[1][1] = m_elem[0][0];
	m_type |= MT_ROTATION;
}


int 
Matrix::invert()
{
	switch(m_type){
	case MT_NONE:
		break;

	case MT_SCALE:
		m_elem[0][0] = 1.0f / m_elem[0][0];
		m_elem[1][1] = 1.0f / m_elem[1][1];
		m_elem[2][2] = 1.0f / m_elem[2][2];
		break;

	case MT_TRANSLATION:
		m_elem[0][3] = -m_elem[0][3];
		m_elem[1][3] = -m_elem[1][3];
		m_elem[2][3] = -m_elem[2][3];
		break;

	case MT_ROTATION:
		transpose();
		break;

	case MT_UNKNOWN:
	default:
		return -1;
	}
	return 0;
}



void 
Matrix::multiply(const Matrix & m2)
{
	float t[4];

	for(int i = 0; i < 3; ++i){
		memcpy(t, m_elem[i], sizeof(float)*4);
		m_elem[i][0] = t[0]*m2.m_elem[0][0] + 
						t[1]*m2.m_elem[1][0] + 
						t[2]*m2.m_elem[2][0];

		m_elem[i][1] = t[0]*m2.m_elem[0][1] + 
						t[1]*m2.m_elem[1][1] + 
						t[2]*m2.m_elem[2][1];

		m_elem[i][2] = t[0]*m2.m_elem[0][2] + 
						t[1]*m2.m_elem[1][2] + 
						t[2]*m2.m_elem[2][2];

		m_elem[i][3] = t[0]*m2.m_elem[0][3] + 
						t[1]*m2.m_elem[1][3] + 
						t[2]*m2.m_elem[2][3] + t[3];
  }
}


Vector & 
Matrix::multiply(const Vector  & v, Vector & dest) const
{
	dest.init( v.dot4(m_elem[0]), v.dot4(m_elem[1]), v.dot4(m_elem[2]) );
	return dest;
}

rmsmesh::Point & 
Matrix::multiply(const rmsmesh::Point & v, rmsmesh::Point & dest) const
{
	dest.init( v.dot4(m_elem[0]), v.dot4(m_elem[1]), v.dot4(m_elem[2]) );
	return dest;
}

void 
Matrix::multiply(float * vec)
{
  float t[4];
  memcpy(t, vec, sizeof(float)*4);
  vec[0] = m_elem[0][0]*t[0] + m_elem[0][1]*t[1] + m_elem[0][2]*t[2] + m_elem[0][3]*t[3];
  vec[1] = m_elem[1][0]*t[0] + m_elem[1][1]*t[1] + m_elem[1][2]*t[2] + m_elem[1][3]*t[3];
  vec[2] = m_elem[2][0]*t[0] + m_elem[2][1]*t[1] + m_elem[2][2]*t[2] + m_elem[2][3]*t[3];
}


string
Matrix::toString() const
{
	char buf[256];
    sprintf(buf, "[ %5.5f %5.5f %5.5f %5.5f]\n[ %5.5f %5.5f %5.5f %5.5f]\n[ %5.5f %5.5f %5.5f %5.5f]\n", 
			 m_elem[0][0], m_elem[0][1], m_elem[0][2], m_elem[0][3],
			 m_elem[1][0], m_elem[1][1], m_elem[1][2], m_elem[1][3],
			 m_elem[2][0], m_elem[2][1], m_elem[2][2], m_elem[2][3] );
	return string(buf);
}


Transformation::Transformation()
  : ml_current_transform(),
    ml_current_inverse()
{
}


Transformation::~Transformation()
{
}


bool 
Transformation::addMatrix(const Matrix & m)
{
  Matrix minverse(m);
  if(minverse.invert())
    return false;

  ml_current_transform.multiply(m);
  minverse.multiply(ml_current_inverse);
  ml_current_inverse = minverse;
  return true;
}
  
bool 
Transformation::addMatrixPremultiply(const Matrix & m)
{
  Matrix minverse(m);
  Matrix mnotinverse(m);
  if(minverse.invert())
    return false;

  mnotinverse.multiply(ml_current_transform);
  ml_current_transform = mnotinverse;
  ml_current_inverse.multiply(minverse);
  return true;
}
