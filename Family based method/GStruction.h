#pragma once
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
class GStruction
{
public:
	GStruction();
	int GenID;
	double P_value;
	double LR;
	int it;
	dense_vector<double> b_value;
};

struct InputData
{
	mtl::dense_vector<double> phe;
	mtl::dense2D<double> Xinfo;
	mtl::dense2D<double> X_LRT;
	mtl::dense2D<double> triadsGene;
};

struct GaPMatrix
{
	InputData InData;
	
	GStruction Tmarker;
};

;

struct Procedure
{
	double nowProce;
	double TotalP;
};

struct GSPar {
	dense_vector<double> min;
	dense_vector<double> max;
	double partion = 10;
	double step = 3;
	InputData iData;

};