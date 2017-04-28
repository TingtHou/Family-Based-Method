//////////////////////////////////////////////////////////////////////////
////////  The class is used to calculate the inverse matrix
///////	  LU, QR and SVD are all in this class
///////   QR has some problem
///////   I intensely recommend to use SVD mod, like :
///////							InvertMatrixSVD(OriMatrix,InvMatrix)
///////
////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Matrix_calc.h"
#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

// REMEMBER to update "lu.hpp" header includes from boost-CVS 
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 
#include <boost/numeric/ublas/io.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/svd.hpp>
namespace ublas = boost::numeric::ublas;


Matrix_calc::Matrix_calc()
{
}

bool Matrix_calc::InvertMatrixLU(const ublas::matrix<double> &input, ublas::matrix<double> &inverse)
{
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input 
	matrix<double> A(input);
	// create a permutation matrix for the LU-factorization 
	pmatrix pm(A.size1());

	// perform LU-factorization 
	int res = lu_factorize(A, pm);
	if (res != 0) return false;


	// create identity matrix of "inverse" 
	(inverse).assign(identity_matrix<double>(A.size1()));

	// back substitute to get the inverse 
	lu_substitute(A, pm, inverse);

	return true;
}

bool Matrix_calc::InvertMatrixSVD(mtl::dense2D<double>& input, mtl::dense2D<double>& inverse)
{
	using namespace mtl;

	dense2D<double> Ori(input);

	dense2D<double> S(Ori.num_rows(), Ori.num_cols()),
		V(Ori.num_rows(), Ori.num_cols()),
		D(Ori.num_rows(), Ori.num_cols());
	set_to_zero(S);
	set_to_zero(V);
	set_to_zero(D);
	boost::tie(S, V, D) = svd(Ori,0);

	dense2D<double> V_INV(V);
	for (int i = 0; i < V_INV.num_rows(); i++)
	{
		if (V[i][i]>1.0e-10)
		{
			V_INV[i][i] = 1 / V[i][i];
		}
		else
		{
			V_INV[i][i] = 0;
		}
		
	}

	set_to_zero(inverse);
	inverse = S*V_INV*trans(D);


	return true;
}
////////////////////////////////////////incorrect///////////////////////////////////////////////////////////
////bool Matrix_calc::InvertMatrixQR(mtl::dense2D<double>& input, mtl::dense2D<double>& inverse)
////{
////	using namespace mtl;
////
////	dense2D<double> Ori(input);
////
////	std::cout << Ori << std::endl;
////	dense2D<double> Q(Ori.num_rows(), Ori.num_cols()),
////				R_inv(Ori.num_rows(), Ori.num_cols()),
////		R(Ori.num_rows(), Ori.num_cols());
////	set_to_zero(Q);
////	set_to_zero(R_inv);
////	set_to_zero(R);
////	boost::tie(Q, R) = mat::qr_givens(Ori);
////	std::cout << Q << std::endl;
////	
////	std::cout << R << std::endl;
////	std::cout << Q*trans(R) << std::endl;
////	R_inv = inv(trans(R)*R)*trans(R);
////	inverse = R_inv*trans(Q);
////	std::cout << inverse*Ori << std::endl;
////	return true;
////}
//////////////////////////////////////////////////////////////////////////////////////////////////////
#endif //INVERT_MATRIX_HPP