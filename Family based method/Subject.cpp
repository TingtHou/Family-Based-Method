//////////////////////////////////////////////////////////////////////////
////////  The class is used to simulate data of one individual
///////   this class is not using only.
///////   it has been contained in the Simulate.cpp
///////   
///////
////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include "Subject.h"
#include <cmath>
#include <fstream>
#include <random>
#include <windows.h>
//std::fstream f;


void Subject::GenerationOfGEN()
{
	GenerationOfCausalGEN();

	for (int i = 0; i < nSNPs; i++)
	{
		ParentSNPGene[i] = dense_vector<double>(4);
		dense_vector<double> Parent1snp = GenerationParentGEN(i);
		ParentGen[i]= Parent1snp;
	//	f << ParentSNPGene[i][0] << CausalofM[0] << "\t" << ParentSNPGene[i][1] << CausalofM[1] << "\t"<<ParentSNPGene[i][2] << CausalofP[0] << "\t" << ParentSNPGene[i][3] << CausalofP[1] << "\t";
		double child = GenerationChildGEN(Parent1snp,i);
		ChildGen[i]= child;
		
	
	}
	

}

void Subject::GenerationOfCausalGEN()
{
	CausalofM = dense_vector<double>(2);
	
	double alleleM1 = -1;
	double alleleM2 = -1;
	RandomDouble rd;
	if (rd.generation01() < MAFofQ)
	{
		alleleM1 = 0;
	}
	else
	{
		alleleM1 = 1;
	}
	if (rd.generation01() < MAFofQ)
	{
		alleleM2 = 0;
	}
	else
	{
		alleleM2 = 1;
	}
	assert((alleleM1 != -1) & (alleleM2 != -1));
	CausalofM[0] = alleleM1;
	CausalofM[1] = alleleM2;

	CausalofP = dense_vector<double>(2);

	double alleleP1 = -1;
	double alleleP2 = -1;
	if (rd.generation01() < MAFofQ)
	{
		alleleP1 = 0;
	}
	else
	{
		alleleP1 = 1;
	}
	if (rd.generation01() < MAFofQ)
	{
		alleleP2 = 0;
	}
	else
	{
		alleleP2 = 1;
	}
	assert((alleleP1 != -1) & (alleleP2 != -1));
	CausalofP[0] = alleleP1;
	CausalofP[1] = alleleP2;

}

void Subject::GenerationOfPHE()
{
	double n;
	int num_Upper = CausalofChild[0] + CausalofChild[1];

	
	n = m + ((num_Upper-1)?(num_Upper-1)*GenEffect:num_Upper*d);

//	n = m + (num_Upper-1)*GenEffect;
/*
	dense_vector<double, ROW_VECTOR> par = Paremter();
	n = par*var;*/
	double p = 1 / (1 + exp(-n));


	RandomDouble rd;
	double crd = rd.generation01();
	if (crd<p)
	{
		phe = 1;
	}
	else
	{
		phe = 0;
	}
}
//the Generated markers has been stored in the vector Parent, Parent[0] is maternal marker,  Parent[1] is paternal marker
dense_vector<double> Subject::GenerationParentGEN(int i)
{
	
	
	
	dense_vector<double> Parent(2);
	double alleleM1 = -1;
	double alleleM2 = -1;
	double alleleP1 = -1;
	double alleleP2 = -1;
	RandomDouble rd;
	double P_MQ = (1 - MAFofM)*(1 - MAFofQ) + Delta;
	double P_mQ = MAFofM*(1 - MAFofQ) - Delta;
	double P_Mq = (1 - MAFofM)*MAFofQ - Delta;
	double P_mq = MAFofM*MAFofQ + Delta;
	double P_MM1 = -1;
	double P_MM2 = -1;
	/////////////////////////////maternal/////////////////////////////////
	if (CausalofM[0]==0)
	{
		P_MM1 = P_Mq / MAFofQ;
	//	P_m1 = P_mq / MAFofQ;
	}
	else
	{
		P_MM1 = P_MQ / (1 - MAFofQ);
//		P_m1 = P_mQ / (1 - MAFofQ);
	}
	if (CausalofM[1] == 0)
	{
		P_MM2 = P_Mq / MAFofQ;
		//	P_m1 = P_mq / MAFofQ;
	}
	else
	{
		P_MM2 = P_MQ / (1 - MAFofQ);
		//		P_m1 = P_mQ / (1 - MAFofQ);
	}

	if (rd.generation01()<P_MM1)
	{
		alleleM1 =1;
	}
	else
	{
		alleleM1 = 0;
	}
	if (rd.generation01() < P_MM2)
	{
		alleleM2 = 1;
	}
	else
	{
		alleleM2 = 0;
	}
	ParentSNPGene[i][0] = alleleM1;
	ParentSNPGene[i][1] = alleleM2;
	Parent[0] = alleleM1 + alleleM2;
///////////////////////////paternal////////////////////////////
	double P_MP1 = -1;
	double P_MP2 = -1;
	if (CausalofP[0] == 0)
	{
		P_MP1 = P_Mq / MAFofQ;
	//	P_m2 = P_mq / MAFofQ;
	}
	else
	{
		P_MP1 = P_MQ / (1 - MAFofQ);
//		P_m2 = P_mQ / (1 - MAFofQ);
	}

	if (rd.generation01() < P_MP1)
	{
		alleleP1 = 1;
	}
	else
	{
		alleleP1 = 0;
	}
	if (CausalofP[1] == 0)
	{
		P_MP2 = P_Mq / MAFofQ;
		//	P_m2 = P_mq / MAFofQ;
	}
	else
	{
		P_MP2 = P_MQ / (1 - MAFofQ);
		//		P_m2 = P_mQ / (1 - MAFofQ);
	}

	
	if (rd.generation01() <P_MP2)
	{
		alleleP2 = 1;
	}
	else
	{
		alleleP2 = 0;
	}
	assert((alleleP1 != -1) & (alleleP2 != -1));
	ParentSNPGene[i][2] = alleleP1;
	ParentSNPGene[i][3] = alleleP2;
	Parent[1]= alleleP1 + alleleP2;
	
	
	return Parent;
}
///According to the parental genotype, genotypes of the kids have been generated under the HWE and Mendel's genetic law;
double Subject::GenerationChildGEN(dense_vector<double> ParentAllele, int i)
{
	CausalofChild = dense_vector<double>(2);
	double ChildGen = -1;
	double thershold = -1;

	thershold = Theta;
	RandomDouble rd;
	int ChromIdM=-1;
	double v = rd.generation01();
	if (v<0.5)
	{
		ChromIdM = 0;
	}
	else
	{
		ChromIdM = 1;
	}
	int ChromIdP = -1;
	v = rd.generation01();
	if (v<0.5)
	{
		ChromIdP = 2;
	}
	else
	{
		ChromIdP = 3;
	}
	//f << ChromIdM+1 << "\t" << ChromIdP-1 << "\t";
	//////////////maternal//////////////////////////////////
	//	std::cout << ParentSNPGene[i][0]<<CausalofM[0]<<"\t"<<ParentSNPGene[i][1]<<CausalofM[1] << std::endl;
	double MfromMa = ParentSNPGene[i][ChromIdM];
	double AlleleC1 = MfromMa;
	
	if (rd.generation01() < Theta)
	{
		CausalofChild[0] = CausalofM[ChromIdM ? 0 : 1];

	}
	else
	{
		CausalofChild[0] = CausalofM[ChromIdM];
	}

	

		////////////////paternal///////////////////////
	//	std::cout << ParentSNPGene[i][2] << CausalofP[0] << "\t" << ParentSNPGene[i][3] << CausalofP[1] << std::endl;
	double MfromDad = ParentSNPGene[i][ChromIdP];
	double AlleleC2 = MfromDad;
	if (rd.generation01() < Theta)
	{
		CausalofChild[1] = CausalofP[ChromIdP - 2 ? 0 : 1];

	}
	else
	{
		CausalofChild[1] = CausalofP[ChromIdP - 2];
	}
	ChildGen = AlleleC1 + AlleleC2;
	
	
	
//////////////////////////////////////////////////////////////////////////////
	
//	f << AlleleC1 << CausalofChild[0] << "\t" << AlleleC2 << CausalofChild[1] << std::endl;
//	f.flush();
// 	f.close();
	
	
	return ChildGen;
}



dense_vector<double,ROW_VECTOR> Subject::Paremter()
{
	dense_vector<double,ROW_VECTOR> para(6);
	para[0]= 1;
	int sumP = ParentGen[loci][0] + ParentGen[loci][1];
	int sumA = ParentGen[loci][0] + ParentGen[loci][1] + ChildGen[loci];
	para[1] = 2 * (ChildGen[loci] - 1);
	para[2] = pow(-1, ChildGen[loci] + 1);
	if ((ParentGen[loci][0] == ParentGen[loci][1]) && (ParentGen[loci][0] == 1))
	{
		
	
		para[3] = 4 * (1 - ChildGen[loci]);
		para[4] = 4 * pow(-1, ChildGen[loci]);
		para[5] = 4 * pow(-1, ChildGen[loci] + 1);

	}
	else
	{
		para[3] = sumA % 3 ? 2 * pow(-1, sumA % 3 + 1) : 0;
		para[4] = sumA % 3 ? (sumA % 2 ? 2 : -2) : 0;
		para[5] = 0;
		
	}
//	std::cout << ParentGen[loci][0] << "\t" << ParentGen[loci][1] << "\t" << ChildGen[loci] << "\t:\t" << para << std::endl;
	return para;
}

Subject::Subject(SimPar par)
{
	ninds = par.GetNumofInd();
	nSNPs = par.GetNumofSNP();
	MAFofM = par.GetMAFofM();
	MAFofQ = par.GetMAFofQ();
	SNPs = par.GetLinkSNP();
	Delta = par.GetDelta();
	Theta = par.GetTheta();
	GenEffect = par.GetGenEffect();
	loci = par.GetLinkSNP();
	d = par.getd();
	m = par.getm();
	ChildGen=dense_vector<double>(nSNPs);
	ParentGen= dense_vector<dense_vector<double>>(nSNPs);
	ParentSNPGene = dense_vector<dense_vector<double>>(nSNPs);
	var = par.GetVar();
}

Subject::Subject()
{
}

void Subject::Generation()
{
//	f.open("causal.txt", std::ios::app);
	GenerationOfGEN();
	GenerationOfPHE();
}