//////////////////////////////////////////////////////////////////////////
////////  The class is used to analysis the simulation data
///////
///////
///////
///////
///////
////////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "SearchGene.h"
#include "Logistic.h"
#include <iostream>
#include <iomanip>
#include "boost/thread/mutex.hpp"  
#include "GridSearch.h"
#include "Logit_Jacobian.h"
#include "Grid_Newton.h"
boost::atomic_int TotalCp(0);
SearchGene::SearchGene()
{
	
}

//estimating the parameter and the good-fit of the model using Likelihood Ratio Test
//the function will be very very very very slow because all calculation is in one thread
void SearchGene::StartAnalysis(GaPMatrix *G1Struct,Procedure *bar)
{
	//std::cout << G1Struct->phe << std::endl;
	//	pG1S = G1Struct;
	int it = 0;
	double LR = 0;
	do{
		Logit_Jacobian *lgGLM = new Logit_Jacobian(G1Struct->InData, true, false);
		lgGLM->fitLM();
		//get estimated value 
		G1Struct->Tmarker.b_value = lgGLM->getB();
		//calculating the Likelihood Ratio
		double likelihoodofAlter = lgGLM->getDeviace();// get the likelihood under the alternate hypothesis
														  //	std::cout << G1Struct->Tmarker.b_value;
		double likelihoodofNull = ScoreTest(G1Struct);//get the likelihood under the null hypothesis
		
		delete lgGLM;

		LR = (likelihoodofAlter - likelihoodofNull);//calculating the likelihood ratio, D=-2*(ln(likelihood fot null model)-ln(likelihood fot alternate model))
		it++;
	} while (it < 10 && LR < 0);
		
	
		G1Struct->Tmarker.LR = LR;
		if (LR > 0)
		{
			G1Struct->Tmarker.P_value = CalcPvalue(LR);
		}
		else
		{
			G1Struct->Tmarker.P_value = NULL;
		}
	//	G1Struct->Tmarker.P_value = CalcPvalue(LR);// calculating the P-value under the chi-squared distribution with 3 freedom
		TotalCp = TotalCp + 1;
/*
		boost::mutex lock;
		lock.lock();
		std::cout << "\r\r"<<TotalCp;
		lock.unlock();
		lock.destroy();*/

		boost::mutex lock;
		lock.lock();
		bar->nowProce++;
		double percent = bar->nowProce / bar->TotalP * 100;
		int barline = percent;
		std::cout << "\r\r[ ";
		for (int i = 0; i < barline - 1; i++)
		{
			std::cout << "=";
		}
		std::cout << "= " << std::setw(103 - barline) << "]" << std::setw(5) << std::setprecision(4)<< percent << "%";
		lock.unlock();
		lock.destroy();

/*
		GeneID++;
		///////////////////////////Show a process bar/////////////////////////////////////////
		double percent = (double)GeneID / (double)SR.nSNPs * 100;
		int bar = percent;
		std::cout << "\r\r[ ";
		for (int i = 0; i < bar-1; i++)
		{
			std::cout << "=";
		}
		std::cout << "= " <<std::setw(103-bar)<<"]" << std::setw(3)<< percent << "%";
	
	}
	std::cout << "\n";*/
}

//this function is used as grid search
void SearchGene::StartAnalysis(GaPMatrix *G1Struct, Procedure *bar,bool Grid)
{
	//std::cout << G1Struct->phe << std::endl;
	//	pG1S = G1Struct;
	int it = 1;
	double LR = 0;

	do 
	{
		Grid_Newton *lgGLM = new Grid_Newton(G1Struct->InData, false);
		lgGLM->SetPartion(it);
		lgGLM->fitLM();
		//get estimated value 
		G1Struct->Tmarker.b_value = lgGLM->getB();
		G1Struct->Tmarker.it = lgGLM->it_final;
		//calculating the Likelihood Ratio
		double likelihoodofAlter = lgGLM->getDeviace();// get the likelihood under the alternate hypothesis
													   //	std::cout << G1Struct->Tmarker.b_value;
		double likelihoodofNull = ScoreTest(G1Struct, true);//get the likelihood under the null hypothesis

		delete lgGLM;

		LR = 2*(likelihoodofAlter - likelihoodofNull);//calculating the likelihood ratio, D=-2*(ln(likelihood fot null model)-ln(likelihood fot alternate model))
		it++;
	} while (LR<-0.000001&&it<6);
	


	if (LR<0&&LR>-0.000001)
	{
		LR = 0;
		G1Struct->Tmarker.it = 1;
	}
	std::cout << G1Struct->Tmarker.b_value <<"\t"<< G1Struct->Tmarker.it<<std::endl;
	G1Struct->Tmarker.LR = LR;
	G1Struct->Tmarker.P_value = CalcPvalue(LR);// calculating the P-value under the chi-squared distribution with 3 freedom
	TotalCp = TotalCp + 1;


}

int SearchGene::StartAnalysis(GaPMatrix *G1Struct)
{
	
	int count = 0;
	double LR = -10;
	double likelihoodofNull=0;
	dense_vector<double> nowE;

	do 
	{
		Logit_Jacobian *lgGLM = new Logit_Jacobian(G1Struct->InData, true, false);
		lgGLM->fitLM();
		//get estimated value 
		nowE = lgGLM->getB();
		std::cout << nowE << std::endl;
		//calculating the Likelihood Ratio
		double likelihoodofAlter = lgGLM->getDeviace();// get the likelihood under the alternate hypothesis

		likelihoodofNull = ScoreTest(G1Struct);//get the likelihood under the null hypothesis
		G1Struct->Tmarker.it = lgGLM->it_final;
		delete lgGLM;
		LR = 2*(likelihoodofAlter - likelihoodofNull);//calculating the likelihood ratio, D=-2*(ln(likelihood fot null model)-ln(likelihood fot alternate model))

		
		count++;

	} while (LR<-0.000001&&count<3);
	G1Struct->Tmarker.LR = LR;
	if (LR<0 && LR>-0.000001)
	{
		LR = abs(LR);
		G1Struct->Tmarker.it = 1;
	}
	if (LR<-0.000001)
	{
		return 0;
	}
	G1Struct->Tmarker.P_value = CalcPvalue(LR);
	G1Struct->Tmarker.b_value = nowE;
/*
	for (int j=0;j<size(nowE);j++)
	{
		G1Struct->Tmarker.b_value[j] = nowE[j];
	}
	G1Struct->Tmarker.b_value[3] = T;*/
	// calculating the P-value under the chi-squared distribution with 3 freedom
	return 1;
}
//This function is used to estimating the variables under the H0, in which the 0=0.5.
//In this function, the pi*alpha, pi*beta, pi^2*beta are considered as 0, for the six parameters model is only considered here.
//after the function has been done, a likelihood value under null hypothesis will been return.

double SearchGene::ScoreTest(GaPMatrix *G1Struct)
{
	 
	
	Logit_Jacobian *lgGLM = new Logit_Jacobian(G1Struct->InData, true,true);
	lgGLM->fitLM();

	double d = lgGLM->getDeviace();
	delete lgGLM;
	return d;
}


double SearchGene::ScoreTest(GaPMatrix *G1Struct,bool gre)
{

	
	Grid_Newton *lgGLM = new Grid_Newton(G1Struct->InData, true);
	lgGLM->fitLM();

	double d = lgGLM->getDeviace();
	delete lgGLM;
	return d;
}

//calculating the P value under the chi squared distribution with  1 freedom
double SearchGene::CalcPvalue(double LRT)
{
	boost::math::chi_squared dist(1);
	double p = boost::math::cdf(boost::math::complement(dist, LRT));
	return p;
}

