#pragma once
#include "GLM.h"
class Grid_Newton :
	public GLM
{
public:
	Grid_Newton(InputData& inData,  bool h);
	void fitLM();
	void SetPartion(int parten);
	~Grid_Newton();
protected:
	double Newton();
	void FitAlter();
	void FitNull();
	void Initial();
	int parten = 10;
	double mintheta =0;
	double maxtheta = 1;
	dense_vector<double> b_NT;
	double BestMLE = -INFINITY;
	dense_vector<double> Best_b_NT;
	double Best_Theta = 0;
	int bestIt = 0;
	
};

