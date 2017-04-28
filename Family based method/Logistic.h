#pragma once
#include "GLM.h"
#include "Matrix_calc.h"
#include <math.h>


class Logistic :public  GLM
{
public:
	Logistic(InputData& inData, bool J, bool h);
	~Logistic();
	void fitLM();

protected:
	
	 
	
	
};

