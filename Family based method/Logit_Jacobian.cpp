#include "stdafx.h"
#include "Logit_Jacobian.h"

Logit_Jacobian::Logit_Jacobian(InputData & inData, bool J, bool h):GLM(inData,J,h)
{

}


Logit_Jacobian::~Logit_Jacobian()
{
}

void Logit_Jacobian::fitLM()
{
	

	if (nind == 0 || np == 0)
	{
		return;
	}
	bool converge = false;
	int it = 0;

	dense2D<double> InfoMatrix_inv;
	double thresold = 1e-6;
	double landa = 1;
	dense2D<double> I(size(b), size(b));
	set_to_zero(I);
	for (int i=0;i<I.num_cols();i++)
	{
		I[i][i] = 1;
	}
	bool sign = false;
 	while (!converge && it < 200)
	{
	

		dense_vector<double> t;
		set_to_zero(t);
		t = calcY(h);
	
		if (_isnan(t[0]))
		{
			std::cout << "NAN" << std::endl;
			std::cout << b << std::endl;
		}

		dense_vector<double> u(size(t));
		set_to_zero(u);
		set_to_zero(V);
		
	
		//// Determine p and V
		dense2D<double> W(nind, nind);
		set_to_zero(W);
		//dense2D<double> D(nind, nind);
		//set_to_zero(D);
		for (int i = 0; i < nind; i++)
		{

			double tmp = 1 / (1 + exp(-t[i]));
			u[i] = tmp;
			V[i][i] = tmp*(1 - tmp);
			W[i][i] = 1 / (tmp*(1 - tmp));
	
		}
		

		dense2D<double> J;
		
		J = CalcJacobian();
		
		
		// Update coefficients
		// b <- b +  solve( t(J) %*% J ) %*% t(J) %*% ( y - p ) 

		dense2D<double> Jt;
	
		Jt = trans(J);

/*
		dense2D<double> JtV;
		JtV = Jt*V;


		dense2D<double> JtVJ;
		set_to_zero(JtVJ);
		JtVJ = JtV*J;
	//	std::cout << JtJ << std::endl;
		Score = Jt*(Y - u);
		dense_vector<double> newB;

		dense2D<double> JtVJ_INV = dense2D<double>(JtVJ.num_rows(), JtVJ.num_cols());
		set_to_zero(JtVJ_INV);
		Matrix_calc::InvertMatrixSVD(JtVJ, JtVJ_INV);
		newB = JtVJ_INV*Score;*/

		dense2D<double> JtJ;
		JtJ = Jt*J;


		
		//	std::cout << JtJ << std::endl;
		Score = Jt*(Y - u);
		dense_vector<double> newB;
		dense2D<double> TMP(I.num_rows(), I.num_cols());
		set_to_zero(TMP);
		for (int id = 0; id < I.num_cols(); id++)
		{
			TMP[id][id] = landa*I[id][id];
		}
		dense2D<double> LMA;
		LMA = JtJ +TMP;
		dense2D<double> LMA_INV = dense2D<double>(LMA.num_rows(), LMA.num_cols());
		set_to_zero(LMA_INV);
		Matrix_calc::InvertMatrixSVD(LMA, LMA_INV);
		newB = LMA_INV*Score;

		// Update coefficients, and check for 
		// convergence


		double P_LM=0;
		for (int i = 0; i < size(newB); i++)
		{
			P_LM += newB[i] * newB[i];
		}
		P_LM = sqrt(P_LM);
		if (P_LM>thresold)
		{
			landa = landa / 10;
			//sign = false;
		}
		else
		{
			landa = landa * 10;
		}

		





		double D = 0;
		for (int j = 0; j < size(b); j++)
		{
			D += abs(newB[j]);
		
		}
		b += newB;
		
		if (D < 1e-10) //if b is converged 
		{
			converge = true;

		}
		if (converge != true)
		{
			//bool outRange = false;
			for (int i = 0; i < size(b); i++)
			{
				if (std::abs(b[i]) > 5)
				{
					//	outRange = true;
					initalGLM();
					break;
				}
				if (!h)
				{
					if (b[3]>1||b[3]<0)
					{
						initalGLM();
						break;
					}
					
				}
				
			}

		}

		
		//whether the estimated value is out the range, if true, re-initial the b value and restart the Newton-Raphson step 
		
		// Next iteration
		it++;

		
	}
	it_final = it;
	dense_vector<double, COL_VECTOR> var(6);
	double u = b[0];
	double alpha = b[1];
	double beta = b[2];
	double theta = 0.5;
	set_to_zero(var);
	if (h)
	{
		
		CalcD(b, Deviance);

	}
	else
	{
		theta = b[3];
		var[0] = u;
		var[1] = alpha;
		var[2] = beta;
		var[3] = theta*alpha;
		var[4] = theta*beta;
		var[5] = theta*theta*beta;

		CalcD(var,Deviance);
	}
	


}




dense2D<double> Logit_Jacobian::CalcJacobian()
{
	dense2D<double> JocabianMatrix;
	if (h)
	{
		JocabianMatrix=JacoLRT();
	}
	else
	{
		JocabianMatrix = JacoAlter();
	}
	return JocabianMatrix;
}

dense2D<double> Logit_Jacobian::JacoLRT()
{
	dense2D<double> J(nind, 3);
	set_to_zero(J);
	J = X;
	return J;
}

dense2D<double> Logit_Jacobian::JacoAlter()
{
	dense2D<double> J(nind, 4);
	set_to_zero(J);
	dense_vector<double, COL_VECTOR> var(6);
	double u = b[0];
	double alpha = b[1];
	double beta = b[2];
	double theta = b[3];

	for (int i = 0; i < nind; i++)
	{
		
		J[i][0] = 1;
		int sumP = triadsGene[i][0] + triadsGene[i][1];
		int sumA = triadsGene[i][0] + triadsGene[i][1] + triadsGene[i][2];

		if ((triadsGene[i][0] == triadsGene[i][1]) && (triadsGene[i][1] == 1))
		{


			switch (triadsGene[i][2])
			{
			case 0:
			{
				J[i][1] = -2+4*theta;
				J[i][2] = -1+4*theta-4*theta*theta;

				J[i][3] = 4*alpha+4*beta-8*theta*beta;

			}
			break;
			case 1:
			{
				J[i][1] = 0;
				J[i][2] = 1-4*theta+4*theta*theta;
				J[i][3] = -4*beta+8*theta*beta;

			}
			break;
			case 2:
			{
				J[i][1] = 2 - 4 * theta;
				J[i][2] = -1+4*theta-4*theta*theta;

				J[i][3] = -4*alpha+4*beta-8*theta*beta;

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
				J[i][1] = -2;
				J[i][2] = -1;
				J[i][3] = 0;


			}
			break;
			case 1:
			{
				J[i][1] = 2 * (theta - 1);
				J[i][2] = 2 * theta - 1;

				J[i][3] = 2 * (alpha + beta);

			}
			break;
			case 2:
			{
				J[i][1] = -2 * theta;
				J[i][2] = 1 - 2 * theta;

				J[i][3] = -2 * (alpha + beta);

			}
			break;
			case 3:
			{
				J[i][1] = 0;
				J[i][2] = 1;

				J[i][3] = 0;

			}
			break;
			case 4:
			{
				J[i][1] = 2 * theta;
				J[i][2] = 1 - 2 * theta;

				J[i][3] = 2 * (alpha - beta);

			}
			break;
			case 5:
			{
				J[i][1] = 2 * (1 - theta);
				J[i][2] = -1 + 2 * theta;

				J[i][3] = -2 * (alpha - beta);

			}
			break;
			case 6:
			{
				J[i][1] = 2;
				J[i][2] = -1;

				J[i][3] = 0;

			}
			break;
			}
		}
		
	}

	return J;
}


dense_vector<double> Logit_Jacobian::calcY(bool LRT)
{

	double u = b[0];
	double alpha = b[1];
	double beta = b[2];
	dense_vector<double, COL_VECTOR> var(6);
	double theta = 0.5;
	set_to_zero(var);
	dense_vector<double> phe;
	if (LRT)
	{
		phe = X*b;
	
	}
	else
	{
		theta = b[3];
		var[0] = u;
		var[1] = alpha;
		var[2] = beta;
		var[3] = theta*alpha;
		var[4] = theta*beta;
		var[5] = theta*theta*beta;

		phe = X*var;
	}
	
	return phe;
}

void Logit_Jacobian::CalHess()
{
	Hessian = dense2D<double>(4, 4);

}



