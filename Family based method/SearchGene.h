#pragma once
#include "Simulate.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/atomic.hpp>
#include "RandomDouble.h"
#include "Analysis.h"
#include "GStruction.h"
using namespace mtl;

class SearchGene
{
public:
	SearchGene(); 
//	static GaPMatrix *pG1S;
	static void StartAnalysis(GaPMatrix *G1Struct, Procedure *bar);
	static int StartAnalysis(GaPMatrix *G1Struct);
	static void StartAnalysis(GaPMatrix *G1Struct, Procedure *bar, bool Grid);
protected:
	
	static double CalcPvalue(double LRT);
	static double ScoreTest(GaPMatrix *G1Struct);
	static double ScoreTest(GaPMatrix *G1Struct, bool grid);
};

