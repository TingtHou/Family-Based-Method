#include "Grid_Newton.h"
#include "RandomDouble.h"
#include "Matrix_calc.h"

Grid_Newton::Grid_Newton(InputData& inData, bool h) : GLM(inData, h)
{
	b_NT = dense_vector<double>(3);
	Best_b_NT = dense_vector<double>(3);
	b = dense_vector<double>(4);
	X = dense2D<double>(nind, 3);
	set_to_zero(X);
	set_to_zero(b_NT);
	set_to_zero(Best_b_NT);
	set_to_zero(b);
}

void Grid_Newton::fitLM()
{
	if (h)
	{
		FitNull();
	}
	else
	{
		FitAlter();
	}
}




double Grid_Newton::Newton()
{
	bestIt = 0;
	if (nind == 0 || np == 0)
	{
		return -1;
	}
	bool converge = false;
	int it = 0;

	dense2D<double> InfoMatrix_inv;
	double Oldlikelihood=0, newlikelihood=0;
	while (!converge && it < 200)
	{
		
		dense_vector<double> t;
		set_to_zero(t);

		t = X*b_NT;
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
			u[i] = tmp;
			V[i][i] = tmp*(1 - tmp);
		}

		// Update coefficients
		// b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
		dense2D<double> Xt;

		Xt = trans(X);
		dense2D<double> XtV;
		XtV = Xt*V;
		dense2D<double> InfoMatrix;

		InfoMatrix = XtV*X;
		InfoMatrix_inv = dense2D<double>(InfoMatrix.num_rows(), InfoMatrix.num_cols());
		set_to_zero(InfoMatrix_inv);
		Matrix_calc::InvertMatrixSVD(InfoMatrix, InfoMatrix_inv); //calculating the inverse matrix of the informatiom matrix using SVD; 

		set_to_zero(Score);
		Score = Xt*(Y - u);


		dense_vector<double> newB;
		set_to_zero(newB);
		newB = InfoMatrix_inv*Score;

		// Update coefficients, and check for 
		// convergence

		b_NT += newB;
		
	//	std::cout << " b : " << b_NT << std::endl;
		double Delta = 0;

		for (int j = 0; j < size(b_NT); j++)
		{
			Delta += abs(newB[j]);

		}

		if (Delta < 1e-10) //if b is converged 
		{
			converge = true;

		}
	
		//whether the estimated value is out the range, if true, re-initial the b value and restart the Newton-Raphson step 
		if (converge != true)
		{
			//bool outRange = false;
			for (int i = 0; i < size(b_NT); i++)
			{
				if (std::abs(b_NT[i]) > 10)
				{
					//	outRange = true;
					Initial();
					break;
				}

			}
		}
		// Next iteration
		it++;



	}
	bestIt = it;
	double DV = 0;
	CalcD(b_NT, DV);
	return DV;
}


Grid_Newton::~Grid_Newton()
{
	
}

void Grid_Newton::FitAlter()
{
	

	for (int perm = 0; perm < 10; perm++)
	{
		theta = mintheta;
		BestMLE = -INFINITY;
		int totalPart = 10 * parten;
		double step = (maxtheta - mintheta) / (double)totalPart;
		for (int i=0;i<=totalPart;i++)
		{
		
			UpdateX(theta);
			Initial();
			double nowMLE = 0;
			nowMLE = Newton();
/*
			std::cout << std::setprecision(10) << theta << "\t" << std::setprecision(10) << nowMLE << std::endl;
			std::cout << b_NT << std::endl;*/
			if (nowMLE > BestMLE)
			{
				BestMLE = nowMLE;
				Best_b_NT = b_NT;
				Best_Theta = theta;
				it_final = bestIt;
			}
			theta += step;
		}
		
		theta = Best_Theta;
		maxtheta = std::min(maxtheta, theta + step);
		mintheta = std::max(mintheta, theta - step);
	}

	std::cout << " b : " << Best_b_NT << "\t" << Best_Theta << std::endl;
	
	b[3] = theta;
	for (int i = 0; i < size(Best_b_NT); i++)
	{
		b[i] = Best_b_NT[i];
	}
	CalcD(Best_b_NT, Deviance);


}

void Grid_Newton::FitNull()
{
	theta = 0.5;
	UpdateX(theta);
	Initial();
	double nowMLE = Newton();
	BestMLE = -INFINITY;
	if (nowMLE > BestMLE)
	{
		BestMLE = nowMLE;
		Best_b_NT = b_NT;
		Best_Theta = theta;
	}
	std::cout << " b : " << Best_b_NT << std::endl;
	CalcD(Best_b_NT,Deviance);

}

void Grid_Newton::Initial()
{
	V = dense2D<double>(nind, nind);
	set_to_zero(V);
	b_NT[0] = 1;
	b_NT[1] = 0;
	b_NT[2] = 0;
	/*RandomDouble rd;
	for (int i=0;i<size(b_NT);i++)
	{
		b_NT[i] = rd.generation01();
	}*/
//	
}

void Grid_Newton::SetPartion(int parten)
{
	this->parten = parten;
}

