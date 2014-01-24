#include "WmlLinearSystemExt.h"

#include "rmsdebug.h"
#include <limits>
#include <WmlMath.h>
#include <math.h>

using namespace Wml;


template <class Real>
bool LinearSystemExt<Real>::LUDecompose( Wml::GMatrix<Real> & mMatrix, LUData & vDecomposition )
{
	int nRows = mMatrix.GetRows();
	if ( nRows != mMatrix.GetColumns() )
		return false;

	vDecomposition.vPivots.resize(nRows);
	vDecomposition.mLU = Wml::GMatrix<Real>( nRows, nRows );
	std::vector<int> & vPivots = vDecomposition.vPivots;
	Wml::GMatrix<Real> & mLUMatrix = vDecomposition.mLU;

	mLUMatrix = mMatrix;

	Real dRowSwaps = 1;

	// scaling of each row
	std::vector<Real> vScale( nRows, (Real)0 );

	Real dTemp;
	for ( int i = 0; i < nRows; ++i ) {			// find scaling for each row
		Real dLargest = (Real)0;
		for ( int j = 0; j < nRows; ++j ) {
			if ( (dTemp = (Real)fabs( mLUMatrix[i][j] )) > dLargest )
				dLargest = dTemp;
		}
		if ( dLargest == 0 )
			return false;
		vScale[i] = (Real)1.0 / dLargest;
	}

	int niMax = 0;
	for ( int j = 0; j < nRows; ++j ) {		// loop over columns (Crout's method)

		// not entirely sure 
		for ( int i = 0; i < j; ++i ) {
			Real dSum = mLUMatrix[i][j];
			for ( int k = 0; k < i; ++k )
				dSum -= mLUMatrix[i][k] * mLUMatrix[k][j];
			mLUMatrix[i][j] = dSum;
		}

		// find largest pivot element
		Real dLargestPivot = (Real)0;
		for ( int i = j; i < nRows; ++i ) {
			Real dSum = mLUMatrix[i][j];
			for ( int k = 0; k < j; ++k  )
				dSum -= mLUMatrix[i][k] * mLUMatrix[k][j];
			mLUMatrix[i][j] = dSum;
			if ( (dTemp = vScale[i] * (Real)fabs(dSum)) > dLargestPivot ) {
				dLargestPivot = dTemp;
				niMax = i;
			}
		}

		// swap rows if pivot is in another column
		if ( j != niMax ) {
			for ( int k = 0; k < nRows; ++k ) {
				Real dSwap = mLUMatrix[niMax][k];
				mLUMatrix[niMax][k] = mLUMatrix[j][k];
				mLUMatrix[j][k] = dSwap;
			}
			dRowSwaps = -dRowSwaps;
			vScale[niMax] = vScale[j];
		}

		vPivots[j] = niMax;
		if ( mLUMatrix[j][j] == 0 )
			mLUMatrix[j][j] = Wml::Math<Real>::EPSILON;

		if ( j != nRows-1 ) {
			Real dScale = (Real)1.0 / mLUMatrix[j][j];
			for ( int i = j+1; i < nRows; ++i )
				mLUMatrix[i][j] *= dScale;
		}
	}
		
	return true;
}

template <class Real>
bool LinearSystemExt<Real>::LUBackSub( LUData & vDecomposition, 
									   const Wml::GVector<Real> & vRHS,
									   Wml::GVector<Real> & vSolution )
{
	std::vector<int> & vPivots = vDecomposition.vPivots;
	Wml::GMatrix<Real> & mLUMatrix = vDecomposition.mLU;
	int nRows = mLUMatrix.GetRows();
	if ( mLUMatrix.GetRows() != nRows || vPivots.size() != nRows )
		return false;

	// probably there are more efficient ways to do this...
	vSolution = vRHS;

	int nNonVanish = std::numeric_limits<int>::max();
	for ( int i = 0; i < nRows; ++i ) {
		int nPivot = vPivots[i];
		Real dSum = vSolution[nPivot];
		vSolution[nPivot] = vSolution[i];
		if ( nNonVanish != std::numeric_limits<int>::max() ) {
			for ( int j = nNonVanish; j <= i-1; ++j )
				dSum -= mLUMatrix[i][j] * vSolution[j];
		} else if (dSum)
			nNonVanish = i;
		vSolution[i] = dSum;
	}
	for ( int i = nRows-1; i >= 0; --i ) {
		Real dSum = vSolution[i];
		for ( int j = i+1; j < nRows; ++j )
			dSum -= mLUMatrix[i][j] * vSolution[j];
		vSolution[i] = dSum / mLUMatrix[i][i];
	}

	return true;
}



#define RADIX 2.0		// this should be right for IEEE floating point
template <class Real>
bool LinearSystemExt<Real>::Balance( Wml::GMatrix<Real> & mMatrix )
{
	int nRows = mMatrix.GetRows();
	if ( nRows != mMatrix.GetColumns() )
		return false;

	Real dSqrRdx = RADIX*RADIX;
	bool bDone = false;
	while (!bDone) {
		bDone = true;

		for ( int i = 0; i < nRows; ++i ) {

			Real dR = (Real)0.0;
			Real dC = (Real)0.0;

			for ( int j = 0; j < nRows; ++j ) {
				if ( j != i ) {
					dC += (Real)fabs( mMatrix[j][i] );
					dR += (Real)fabs( mMatrix[i][j] );
				}
			}

			if ( dC && dR ) {
				Real dG = dR / (Real)RADIX;
				Real dF = (Real)1.0;
				Real dS = dC + dR;
				while ( dC < dG ) {
					dF *= (Real)RADIX;
					dC *= dSqrRdx;
				}
				dG = dR * (Real)RADIX;
				while ( dC > dG ) {
					dF /= (Real)RADIX;
					dC /= dSqrRdx;
				}
				if ( (dC+dR) / dF < (Real)0.95 * dS ) {
					bDone = false;
					dG = (Real)1.0 / dF;
					for ( int j = 0; j < nRows; ++j )
						mMatrix[i][j] *= dG;
					for ( int j = 0; j < nRows; ++j )
						mMatrix[j][i] *= dF;
				}
			}
		}
	}
	return true;
}


#define SIGN(a,b) ( (b) >= (Real)0 ? (Real)fabs(a) : -(Real)fabs(a) )

template <class Real>
bool LinearSystemExt<Real>::QREigenValues( Wml::GMatrix<Real> & mMatrix, Wml::Vector2<Real> * pEigenValues )
{
	static const int nMaxIters = 30;
	int nRows = mMatrix.GetRows();
	if ( nRows != mMatrix.GetColumns() )
		return false;

	Real dANorm = (Real)0;
	for ( int i = 0; i < nRows; ++i ) {
		for ( int j = std::max(i-1,0); j < nRows; ++j )
			dANorm += (Real)fabs( mMatrix[i][j] );
	}

	int nNN = nRows-1;	// ?? this is # of eigenvalues. but also used for array indexing? trouble...
	int nL = 0;
	Real dT = (Real)0;
	Real dX, dY, dZ, dW, dV, dU, dR, dQ, dP;
	while ( nNN >= 0 ) {
		int nIters = 0;
		do {
			for ( nL = nNN; nL >= 1; nL-- )  {
				Real dS = (Real)fabs( mMatrix[nL-1][nL-1] ) + (Real)fabs( mMatrix[nL][nL] );
				if ( dS == (Real)0 )
					dS = dANorm;
				if ( (Real)fabs(mMatrix[nL][nL-1]) + dS == dS )
					break;
			}

			dX = mMatrix[nNN][nNN];
			if ( nL == nNN ) {					// one root found
				pEigenValues[nNN] = Wml::Vector2<Real>( dX+dT, (Real)0 );
				--nNN;
			} else {
				dY = mMatrix[nNN-1][nNN-1];
				dW = mMatrix[nNN][nNN-1] * mMatrix[nNN-1][nNN];
				if ( nL == (nNN-1)) {			// two roots found
					Real dP = (Real)0.5 * (dY - dX);
					Real dQ = dP * dP + dW;
					Real dZ = (Real)sqrt(fabs(dQ));
					dX += dT;

					if ( dQ > (Real)0 ) {		// real pair
						dZ = dP + SIGN(dZ,dP);
						pEigenValues[nNN-1] = pEigenValues[nNN] = Wml::Vector2<Real>(dX + dZ, (Real)0);
						if ( dZ )
							pEigenValues[nNN].X() = (Real)(dX - dW/dZ);

					} else {					// complex pair
						pEigenValues[nNN] = Wml::Vector2<Real>( dX + dP, dZ );
						pEigenValues[nNN-1] = Wml::Vector2<Real>( dX  + dP, -dZ );
					}
					nNN -= 2;

				} else {				// no roots found - continue iteration
					if  ( nIters == nMaxIters )
						return false;
					if  ( nIters == 10 || nIters == 20 ) {
						dT += dX;
						for ( int i = 0; i < nNN; ++i )
							mMatrix[i][i] -= dX;
						Real dS = (Real)fabs( mMatrix[nNN][nNN-1]) + (Real)fabs( mMatrix[nNN-1][nNN-2] );
						dY = dX = (Real)0.7 * dS;
						dW = (Real)-0.4375 * dS * dS;
					}
					++nIters;
					int m;
					for ( m = (nNN-2); m >= nL; m-- ) {
						dZ = mMatrix[m][m];
						dR = dX - dZ;
						Real dS = dY - dZ;
						dP = (dR*dS - dW) / mMatrix[m+1][m] + mMatrix[m][m+1];
						dQ = mMatrix[m+1][m+1] - dZ - dR - dS;
						dR = mMatrix[m+2][m+1];
						dS = (Real)fabs(dP) + (Real)fabs(dQ) + (Real)fabs(dR);
						dP /= dS;
						dQ /= dS;
						dR /= dS;
						if ( m == nL )
							break;
						dU = (Real)fabs( mMatrix[m][m-1] ) * (Real)( fabs(dQ) + fabs(dR) );
						dV = (Real)fabs(dP) * (Real)( fabs(mMatrix[m-1][m-1]) + fabs(dZ) + fabs(mMatrix[m+1][m+1]) );
						if ( (dU+dV) == dV )
							break;
					}
					for ( int i = m+2; i <= nNN; ++i ) {
						mMatrix[i][i-2] = (Real)0;
						if ( i != m+2 )
							mMatrix[i][i-3] = (Real)0;
					}
					for ( int k = m; k <= nNN-1; ++k ) {
						if ( k != m ) {
							dP = mMatrix[k][k-1];
							dQ = mMatrix[k+1][k-1];
							dR = (Real)0;
							if ( k != (nNN-1) )
								dR = mMatrix[k+2][k-1];
							dX = (Real)( fabs(dP) + fabs(dQ) + fabs(dR) );
							if ( dX != (Real)0 ) {
								dP /= dX;
								dQ /= dX;
								dR /= dX;
							}
						}

						Real dS = SIGN( (Real)sqrt(dP*dP + dQ*dQ + dR*dR), dP );
						if ( dS != (Real)0 ) {
							if ( k == m ) {
								if ( nL != m )
									mMatrix[k][k-1] = -mMatrix[k][k-1];
							} else 
								mMatrix[k][k-1] = -dS * dX;
							dP += dS;
							dX = dP / dS;
							dY = dQ / dS;
							dZ = dR / dS;
							dQ /= dP;
							dR /= dP;
							for ( int j = k; j <= nNN; ++j ) {
								dP = mMatrix[k][j]  + dQ * mMatrix[k+1][j];
								if ( k != (nNN-1) ) {
									dP += dR * mMatrix[k+2][j];
									mMatrix[k+2][j] -= dP*dZ;
								}
								mMatrix[k+1][j] -= dP*dY;
								mMatrix[k][j] -= dP*dX;
							}
							int nMin = ( nNN < k+3 ) ? nNN : k+3;
							for ( int i = nL; i <= nMin; ++i ) {
								dP = dX * mMatrix[i][k] + dY * mMatrix[i][k+1];
								if ( k != (nNN-1)) {
									dP += dZ * mMatrix[i][k+2];
									mMatrix[i][k+2] -= dP*dR;
								}
								mMatrix[i][k+1] -= dP*dQ;
								mMatrix[i][k] -= dP;
							}
						}
					}
				}
				
			}	


		} while (nL < nNN-1);
	}

	return true;
}




namespace Wml {

template class LinearSystemExt<float>;
template class LinearSystemExt<double>;

}