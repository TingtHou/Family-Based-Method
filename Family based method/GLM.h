#pragma once

#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/svd.hpp>
#include "GStruction.h"
using namespace mtl;

#define COL_VECTOR mtl::vec::parameters<mtl::tag::col_major>
#define ROW_VECTOR mtl::vec::parameters<mtl::tag::row_major>
class GLM
{
public:
	GLM(InputData& inData, bool J,bool LRT);
	GLM(InputData & inData, bool LRT);
	dense2D<double> getVarianceMatrix();
	dense_vector<double, COL_VECTOR> getB();
	double getDeviace();
	double getTheta();
	int it_final = 0;
protected:
	dense2D<double> V;
	InputData Data;
	dense_vector<double,COL_VECTOR> b;
	dense_vector<double, COL_VECTOR> Y;
	dense2D<double> X;
	dense_vector<double, COL_VECTOR> Score;
	int nind = 0;
	int np = 0;
	//double likelihood = 0;
	bool Joca = false;
	bool h = false;
	dense2D<int> triadsGene;
	double Deviance = 0;
	double theta = 0.5;
//	void GLM::uBlastoMTL(ublas::matrix<double>& max, dense2D<double>& translate);
//	dense_vector<double,COL_VECTOR> uBlastoMTL(ublas::vector<double>& v);
	void CalcD(dense_vector<double, COL_VECTOR> Ee, double &like);
	virtual void fitLM() = 0;
	void initalGLM();
	void UpdateX(double NewTheta);
	
};