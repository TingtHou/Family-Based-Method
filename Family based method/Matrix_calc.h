#pragma once
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/utility.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/svd.hpp>
namespace ublas = boost::numeric::ublas;
class Matrix_calc
{
public:
	Matrix_calc();
	static bool InvertMatrixLU(const ublas::matrix<double> &input, ublas::matrix<double> &inverse);
	static bool InvertMatrixSVD(mtl::dense2D<double>& input, mtl::dense2D<double>& inverse);
	//static bool Matrix_calc::InvertMatrixQR(mtl::dense2D<double>& input, mtl::dense2D<double>& inverse);
};
