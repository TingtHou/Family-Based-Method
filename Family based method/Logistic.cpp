///////////////////////////////////////////////////////////////////////////////////////////
////////  The class is used to estimate the variable using the Newton-Raphson method
///////    this class can be used like this:
///////              Logistic *a = new Logistic(Y,X),  Y is a vector and X is a matrix
///////                                                Y,X are all using the form of boost
///////              a->fitLM(), starting estimation
///////				 a->getB(),  get the estimated values
////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include "Logistic.h"
#include <fstream>
#include "RandomDouble.h"
using namespace mtl;

Logistic::Logistic(InputData& inData,bool J,bool h) : GLM(inData,J,h)
{

}

Logistic::~Logistic()
{

}
//this function is used to start the Newton-Raphson estimating
void Logistic::fitLM()
{

	
	if (nind==0||np==0)
	{
		return;
	}
	bool converge = false;
	int it = 0;
	
	dense2D<double> InfoMatrix_inv;
	
	while ( !converge && it<200 )
	{
	
		dense_vector<double> t;
		set_to_zero(t);
		
		t = X*b;
		if (_isnan(t[0]))
		{
			std::cout << "NAN" << std::endl;
		}
		
		dense_vector<double> u(size(t));
		set_to_zero(u);
		set_to_zero(V);
		
		//// Determine p and V
		for (int i = 0; i < nind; i++)
		{
			
			double tmp = 1 / (1 + exp(-t[i]));
			u[i]=tmp;
			V[i][i]=tmp*(1-tmp);
		}
		
		// Update coefficients
		// b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
		dense2D<double> Xt;
	
		Xt= trans(X);
		dense2D<double> XtV;
		XtV=Xt*V;
		dense2D<double> InfoMatrix;
		
		InfoMatrix=XtV*X;
		InfoMatrix_inv=dense2D<double>(InfoMatrix.num_rows(), InfoMatrix.num_cols());
		set_to_zero(InfoMatrix_inv);
		Matrix_calc::InvertMatrixSVD(InfoMatrix, InfoMatrix_inv); //calculating the inverse matrix of the informatiom matrix using SVD; 
	
		set_to_zero(Score);
		Score = Xt*(Y - u);

		
		dense_vector<double> newB;
		set_to_zero(newB);
		newB= InfoMatrix_inv*Score;
	
		// Update coefficients, and check for 
		// convergence
		
		b += newB;
	
	//	std::cout << " b : " << b << std::endl;
		double Delta = 0;
		
		for (int j = 0; j < size(b); j++)
		{
			Delta += abs(newB[j]);
		
		}
	
		if (Delta < 1e-10) //if b is converged 
		{
			converge = true;

		}

		//whether the estimated value is out the range, if true, re-initial the b value and restart the Newton-Raphson step 
		if (converge!=true)
		{
			//bool outRange = false;
			for (int i=0;i<size(b);i++)
			{
				if (std::abs(b[i])>10)
				{
				//	outRange = true;
					for (int i = 0; i < size(b); i++)
					{
						using namespace boost;
						RandomDouble rd;
						double initial = rd.generation01();
						b[i] = initial;
						it = 0;
					}
					break;
				}

			}
		}
		// Next iteration
		it++;
		
	
	
	}
	it_final = it;
//	std::cout << b << std::endl;
	CalcD(b,Deviance);

}


