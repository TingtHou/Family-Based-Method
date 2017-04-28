#pragma once
#include "GLM.h"
#include "Subject.h"
#include "Matrix_calc.h"
using namespace mtl;
class Logit_Jacobian :public  GLM
{
public:

	Logit_Jacobian(InputData& inData, bool J,bool h);
	~Logit_Jacobian();
	void fitLM();
//	void fitLMH();
	
protected:
	
	dense2D<double> CalcJacobian();
	dense2D<double> JacoLRT();
	dense2D<double> JacoAlter();
// 	dense2D<double> Logit_Jacobian::CalcJacobianH();
	dense_vector<double> calcY(bool LRT);
	dense2D<double> Hessian;
	void CalHess();
//	dense_vector<double> calcYH();
//	void Likelihood();
};

