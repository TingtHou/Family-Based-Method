////////////////////////////////////////////////////////////////////////////////////////
/////////////      This class is grid search
///////////
////////////
//////////////////////////////////////////////////////////////////////


#include "stdafx.h"
#include "GridSearch.h"
#include <algorithm>
#include <list>

GridSearch::GridSearch()
{
}

GridSearch::GridSearch(GSPar* par, bool LRT)
{

	//Set ranges of each parameter  
	min = par->min;
	max = par->max;
	//set precision
	step = par->step;
	// set partion in each iteration  
	partion = par->partion;
	
	Data = par->iData;
	if (LRT)
	{
		E_Result = dense_vector<double>(3);
	}
	else
	{
		E_Result = dense_vector<double>(4);
	}
	
	
	set_to_zero(E_Result);

	StartEstimate(LRT);
	Likelihood(E_Result,LRT,likelihood);

}

mtl::dense_vector<double> GridSearch::getEstimate()
{
	return E_Result;
}

double GridSearch::getLikelihood()
{
	return likelihood;
}
//start to search
void GridSearch::StartEstimate(bool LRT)
{
	assert(size(max) == size(min));
	int l = size(max);
	dense_vector<double> currentSteps(l);
	dense_vector<double> currentMin(l);
	dense_vector<double> currentMax(l);
	dense_vector<double> bestvalue(l);
	dense_vector<double> value(l);
	dense_vector<int> state(l);
	dense_vector<int> maxState(l);

	double bestScore = -INFINITY,Score=0;
	for (int i=0;i<size(maxState);i++)
	{
		maxState[i] = partion;
	}
	for (int i=0;i<size(currentMax);i++)
	{
		currentMax[i] = max[i];
		currentMin[i] = min[i];
	}
	for (int s=0;s<step;s++)
	{

		
		for (int i = 0; i < l; i++)
		{
			currentSteps[i] = (currentMax[i] - currentMin[i]) / partion; //calculate each step value
		}
		for (int i = 0; i < size(state); i++)
		{
			state[i] = 0; //indicate whether steps arrive mix
		}
		while (true)
		{

			for (int i = 0; i < size(value); i++)
			{
				value[i] = currentMin[i] + currentSteps[i] * state[i];
			}
			Likelihood(value, LRT, Score);
			if (Score > bestScore)
			{
				bestScore = Score;
				for (int i = 0; i < size(value); i++)
				{
					bestvalue[i] = value[i];
					
				}
				
			}
			
			if (state == maxState)
			{
				break;
			}
			// this loop is used to generate next search array.
			for (int j = l - 1; j >= 0; j--)
			{
				if (state[j] != partion)
				{
					state[j] = state[j] + 1;
					for (int k = j + 1; k < l; k++)
					{
						state[k] = 0;
					}
					break;
				}
			}
		}
		//all combination has been searched, set the new max value and min value for next step
		for (int i = 0; i < l; i++) 
		{
			currentMin[i] = std::max(
				min[i],
				bestvalue[i] - currentSteps[i]);
			currentMax[i] = std::min(
				max[i],
				bestvalue[i] + currentSteps[i]);

		}
	}
	E_Result = bestvalue;

}
//curValue contains ui,alpha, beta, theta
double GridSearch::CalcSocreSquare(dense_vector<double,COL_VECTOR> &curValue,bool LRT)
{
	//
	dense_vector<double> n;
	dense_vector<double, COL_VECTOR> LinearB;
	dense2D<double> X;
	if (LRT==true)
	{
		LinearB = dense_vector<double, COL_VECTOR>(3);
		LinearB[0] = curValue[0];//ui
		LinearB[1] = curValue[1];//alpha
		LinearB[2] = curValue[2];//beta
		X = Data.X_LRT;
		
	}
	else
	{
		LinearB=dense_vector<double, COL_VECTOR>(6);
		LinearB[0] = curValue[0];//ui
		LinearB[1] = curValue[1];//alpha

		LinearB[2] = curValue[2];//beta
		LinearB[3] = curValue[1] * curValue[3];//alpha*theta
		LinearB[4] = curValue[2] * curValue[3];//beta**theta)
		LinearB[5] = pow(curValue[3], 2) * curValue[2];//square(theta)*beta
		X = Data.Xinfo;
	}
	

	n = X*LinearB;

	
	dense_vector<double> pi(size(n));
	for (int i=0;i<size(n);i++)
	{
		pi[i] = exp(n[i]) / (1 + exp(n[i]));
	}
	dense_vector<double, COL_VECTOR> score(size(n));
	for (int i = 0; i < size(n); i++)
	{
		score[i] = Data.phe[i] - pi[i];
	}
	double S2 = 0;
	for (int i=0;i<size(score);i++)
	{
		S2 += score[i] * score[i];
	}
	return sqrt(S2);
}

///////////////////////////////////////////
void GridSearch::Likelihood(dense_vector<double, COL_VECTOR> &curValue, bool LRT,double &likeli)
{
	likeli = 0;
	dense_vector<double, COL_VECTOR>  p_e;
	dense_vector<double, COL_VECTOR> LinearB;
	dense_vector<double, COL_VECTOR> n;
	dense2D<double> X;
	if (LRT == true)
	{
		LinearB = dense_vector<double, COL_VECTOR>(3);
		LinearB[0] = curValue[0];//ui
		LinearB[1] = curValue[1];//alpha
		LinearB[2] = curValue[2];//beta
		X = Data.X_LRT;
	}
	else
	{
		LinearB = dense_vector<double, COL_VECTOR>(6);
		LinearB[0] = curValue[0];//ui
		LinearB[1] = curValue[1];//alpha

		LinearB[2] = curValue[2];//beta
		LinearB[3] = curValue[1] * curValue[3];//alpha*theta
		LinearB[4] = curValue[2] * curValue[3];//beta**theta)
		LinearB[5] = pow(curValue[3], 2) * curValue[2];//square(theta)*beta
		X = Data.Xinfo;
		
	}


	n = X*LinearB;
	
	p_e = dense_vector<double, COL_VECTOR>(size(n));
	set_to_zero(p_e);
	for (int i = 0; i<size(p_e); i++)
	{
		p_e[i] = 1 / (1 + exp(-n[i]));
	}
	
	for (int i = 0; i<size(p_e); i++)
	{
		double f = Data.phe[i] * log(p_e[i]) + (1 - Data.phe[i])*log(1-p_e[i]); //the probability density function of the Bernoulli distribution, f=p^y*(1-p)^(1-y) 
																//	double f = Y[i]*tmp[i]-log(1+exp(tmp[i]));
		assert(!std::isinf(f));

		likeli += f;
		//likelihood += f;
		assert(!std::isinf(likeli));
	}
	
}