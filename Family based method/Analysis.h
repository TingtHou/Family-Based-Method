#pragma once
#include "Simulate.h"
#include "GStruction.h"
#include <boost/numeric/mtl/mtl.hpp>
using namespace mtl;


class Analysis
{
public:
	Analysis(Simulate SimResult);
	int ThreadRun();
	dense_vector<GaPMatrix*> Result;
protected:
	
	Simulate SR;
	void GetYandX(int Gene, GaPMatrix* ma);
	dense_vector<double> TranslateX(Subject SubInd, int loci);
	dense_vector<double> TranslateXC(Subject SubInd, int loci);
};

