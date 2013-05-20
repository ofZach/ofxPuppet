#ifndef __RMS_WMLEXT_WMLLINEARSYSTEMEXT_H__
#define __RMS_WMLEXT_WMLLINEARSYSTEMEXT_H__

#define WML_ITEM 
#include <WmlGMatrix.h>
#include <WmlGVector.h>
#include <WmlVector2.h>
#include <vector>

namespace Wml
{

template <class Real>
class LinearSystemExt
{
public:

	struct LUData {
		Wml::GMatrix<Real> mLU;
		std::vector<int> vPivots;
	};

	static bool LUDecompose( Wml::GMatrix<Real> & mMatrix, LUData & vDecomposition );
	static bool LUBackSub( LUData & vDecomposition, const Wml::GVector<Real> & vRHS, Wml::GVector<Real> & vSolution );

	//! this function balances mMatrix, changing it in the process
	static bool Balance( Wml::GMatrix<Real> & mMatrix );

	/** QR Algorithm for Real Hessenberg Matrices (mMatrix must be in upper Hessenberg form!). pEigenValues
	    must be same dimension as matrix (which must be square!) **/ 
	static bool QREigenValues( Wml::GMatrix<Real> & mMatrix, Wml::Vector2<Real> * pEigenValues );

protected:


};

typedef LinearSystemExt<double> LinearSystemExtd;

} // namespace Wml


#endif // __RMS_WMLEXT_WMLLINEARSYSTEMEXT_H__