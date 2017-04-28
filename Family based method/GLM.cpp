//////////////////////////////////////////////////////////////////////////
////////  The class is based class of the GLM 
///////	  !!!!!!!This Class can not be used directly, PLEASE USING THE Logistic Class
///////         
///////
///////
///////
////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <fstream>
#include "GLM.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "RandomDouble.h"
#include "GridSearch.h"
#include <boost/thread/thread.hpp>

////////////////////////////////////////////////////////////////////////////////////
///////Initial the GLM class
//////inData is the structure InputData
/////J is an indicator whether use linear regression or non-linear regression
/////   true indicates using non-linear model
////LRT is an signal which is used to indicate whether it is under the null hypothesis. True means it is under the null hypothesis
GLM::GLM(InputData& inData, bool J, bool LRT)
{


	Data = inData;
	Y = inData.phe;
	nind = size(Y);
	
	triadsGene = inData.triadsGene;
	Joca = J;
	h = LRT;
////NOTES: When non-linear regression model is chosen, the X matrix can't be used. A new X matrix should be updated. 
	if (h)
	{
		X = inData.X_LRT;
	}
	else
	{
		X = inData.Xinfo;
	}
	np = X.num_cols();
	initalGLM();
	
}

GLM::GLM(InputData & inData,bool LRT)
{
	Data = inData;
	Y = inData.phe;
	nind = size(Y);

	triadsGene = inData.triadsGene;
	h = LRT;
	////NOTES: When non-linear regression model is chosen, the X matrix can't be used. A new X matrix should be updated. 
}






//Calculating the log-likelihood 
void GLM::CalcD(dense_vector<double,COL_VECTOR> Ee,double &like)
{
	dense_vector<double, COL_VECTOR> p_e,tmp;
	like = 0;
	double test = 0;
	
	tmp= X*Ee;
	p_e = dense_vector<double, COL_VECTOR>(size(tmp));
	set_to_zero(p_e);

	for (int i=0;i<size(tmp);i++)
	{
		double done= 1 / (1 + exp(-tmp[i]));
		
		p_e[i] = done;
	}
	
	for (int i=0;i<size(p_e);i++)
	{
	
		test = Y[i] * log(p_e[i]) + (1 - Y[i])*log(1 - p_e[i]);
	
		like += test;
	}

	//Deviance = Deviance * 2;
	return ;
}

//initial the parameters
void GLM::initalGLM()
{
	V = dense2D<double>(nind, nind);
	RandomDouble rd;
	
	if (Joca)
	{
		if (!h)
		{
			b = dense_vector<double, COL_VECTOR>(4);
		}
		else
		{
			b = dense_vector<double, COL_VECTOR>(3);
		}
		
		set_to_zero(b);
		for (int i = 0; i < size(b); i++)
		{
			using namespace boost;
			RandomDouble rd;
			double initial = rd.generation01();
			b[i] = initial;
		}

		
	//	std::cout << b << std::endl;
	}
	else
	{
		b = dense_vector<double, COL_VECTOR>(np);
		set_to_zero(b);
		for (int i = 0; i < size(b); i++)
		{
			using namespace boost;
			RandomDouble rd;
			double initial = rd.generation01();
			b[i] = initial;
		}
	}
		
/*

	if (Joca==true)
	{
		if (h==false)
		{
			b = dense_vector<double, COL_VECTOR>(4);
		}
		else
		{
			b = dense_vector<double, COL_VECTOR>(3);
		}
		
		for (int i = 0; i < size(b); i++)
		{
			using namespace boost;
			RandomDouble rd;
			double initial = rd.generation();
			b[i] = initial;
		}
		if (b[3]>0.5)
		{
			b[3] = b[3] / 2;
		}
	}
	else
	{
		b = dense_vector<double, COL_VECTOR>(np);
		for (int i = 0; i < size(b); i++)
		{
			using namespace boost;
			RandomDouble rd;
			double initial = rd.generation();
			b[i] = initial;
		
		}
		
	}
	*/
}
////According to the θ estimated by Grid search, a new covariate matrix has to been updated.
void GLM::UpdateX(double NewTheta)
{

	
	dense2D<double> NewX(nind,3);
	set_to_zero(NewX);
	double pi = 1 - 2 * NewTheta;
	for (int i=0;i<nind;i++)
	{
		NewX[i][0] = 1;
		int sumP = triadsGene[i][0] + triadsGene[i][1];
		int sumA = triadsGene[i][0] + triadsGene[i][1] + triadsGene[i][2];
	
		if ((triadsGene[i][0] == triadsGene[i][1]) && (triadsGene[i][1] == 1))
		{
			switch (triadsGene[i][2])
			{
			case 0:
			{
				NewX[i][1] = -2 * pi;
				NewX[i][2] = -pi*pi;

			}
			break;
			case 1:
			{
				NewX[i][1] = 0;
				NewX[i][2] = pi*pi;

			}
			break;
			case 2:
			{
				NewX[i][1] = 2 * pi;
				NewX[i][2] = -pi*pi;
			}
			break;
			}
		}
		else
		{

			switch (sumA)
			{
			case 0:
			{
				NewX[i][1] = -2;
				NewX[i][2] = -1;
			}
			break;
			case 1:
			{
				NewX[i][1] = -(1+pi);
				NewX[i][2] = -pi;

			}
			break;
			case 2:
			{
				NewX[i][1] = -(1 - pi);
				NewX[i][2] = pi;
			}
			break;
			case 3:
			{
				NewX[i][1] = 0;
				NewX[i][2] = 1;
			}
			break;
			case 4:
			{
				NewX[i][1] = (1 - pi);
				NewX[i][2] = pi;

			}
			break;
			case 5:
			{
				NewX[i][1] = 1 + pi;
				NewX[i][2] = -pi;
			}
			break;
			case 6:
			{
				NewX[i][1] = 2;
				NewX[i][2] = -1;
			}
			break;
			}
		}
		
	
	}
	
	set_to_zero(X);
	X = NewX;
	np = X.num_cols();
}

double GLM::getTheta()
{
	return theta;
}

dense2D<double> GLM::getVarianceMatrix()
{
	
	return V;
}

dense_vector<double,COL_VECTOR> GLM::getB()
{
	return b;
}

double GLM::getDeviace()
{
	return Deviance;
}

