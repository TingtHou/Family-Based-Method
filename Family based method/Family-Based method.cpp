// Family-Based method.cpp: 主项目文件。

#include "stdafx.h"
#include "Logistic.h"
#include "GLM.h"
#include <math.h>
#include "SimPar.h"
#include "RandomDouble.h"
#include "Simulate.h"
#include "Analysis.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric//ublas/io.hpp>
#include <fstream>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include "ReadFile.h"
using namespace mtl;
void Trans2PED(Simulate sm);

int main()
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////Read file to analysis////////////////////////////
/*
	ReadFile rf("sim.txt");
	std::cout << "Read File Over" << std::endl;
	Analysis out(rf.getSR());
	std::ofstream result;
	result.open("result.txt");
	std::cout << "Analysis starting" << std::endl;
	out.ThreadRun();
	dense_vector<GaPMatrix*> p = out.Result;


	for (int i = 0; i < size(p); i++)
	{
		result << "SNP" << i << "\t";
		result << p[i]->Tmarker.LR << "\t";
		result << p[i]->Tmarker.P_value << "\t";
		for (int j=0;j<size(p[i]->Tmarker.b_value);j++)
		{
			result << p[i]->Tmarker.b_value[j] << "\t";
		
		}
		result << p[i]->Tmarker.it << std::endl;
		result.flush();
	}

	result.close();
	std::cout << "Finished" << std::endl;
	getchar();
	getchar();
*/


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////simulation////////////////////////////////////////////////////////////////////////
/*


	int indNum = 10000;
	int snpNum = 1;
	int loci = 0;


/ *
	std::cout << "Please input individual number: ";
	std::cin >> indNum;
	std::cout << "\n";
	std::cout << "Please input SNP number: ";
	std::cin >> snpNum;
	std::cout << "\n";
	std::cout << "Please input link loci: ";
	std::cin >> loci;
	std::cout << "\n";
	while (loci > snpNum)
	{
		std::cout << "link loci must be smaller than  SNP number ";
		std::cout << "\n";
		std::cout << "Please input link loci: ";
		std::cin >> loci;
		std::cout << "\n";
	}	* /
	clock_t starttime = clock();
	SimPar a;
	a.SetLinkSNP(loci);
	a.SetDelta(0.2);
	a.SetTheta(0.1);
	a.SetMAFofM(0.3);
	a.SetMAFofQ(0.4);
	a.SetNumofInd(indNum);
	a.SetNumofSNP(snpNum);
	a.Setd(0);
	a.Setm(0.5);
	a.SetGenEffect(1);
	a.CalcParemeter();
	
	Simulate test(a);
	test.Generation();

	test.printout();

	std::ofstream result;

	result.open("result.txt");
	std::cout << "the simulation has been done" << std::endl;
	result << a.toString();
	result << "SNP ID\tP_Value\tu\ta\tb\tpia\tpib\tpi2b\n";
	result.flush();
	Analysis out(test);
	std::cout << "Analysis starting" << std::endl;
	out.ThreadRun();
	dense_vector<GaPMatrix*> p = out.Result;


	for (int i=0;i<size(p);i++)
	{
		result << "SNP" << i << "\t";
		result << p[i]->Tmarker.LR << "\t";
		result << p[i]->Tmarker.P_value << "\t";
		result << p[i]->Tmarker.b_value << "\t";
		result << p[i]->Tmarker.it << "\t";
		result << "\n";
		result.flush();

	}
	
	result.close();
	clock_t endtime = clock();
	std::cout << "Finished" << std::endl;
	std::cout << "Elapsed time : " << (endtime - starttime) / CLOCKS_PER_SEC << " S" << std::endl;
	getchar();
	getchar();
	*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////Type I error////////////////////////////////////////////////////////////////////////////
	clock_t starttime = clock();
	SimPar a;
	a.SetLinkSNP(0);
	a.SetDelta(0);
	a.SetTheta(0.5);
	a.SetMAFofM(0.40);
	a.SetMAFofQ(0.30);
	a.SetNumofInd(1000);
	a.SetNumofSNP(1);
	a.SetGenEffect(2);
	a.Setd(0);
	a.Setm(0);
	a.CalcParemeter();
	std::ofstream result;
	result.open("result.txt");
	result << "P_Value" << std::endl;
	/*result << a.toString();
	result << "SNP ID\tP_Value\tu\ta\tb\tpia\tpib\tpi2b\n";*/
	result.flush();
	for (int times=0;times<500;)
	{
		Simulate test(a);
		test.Generation();
		test.printout();
		Analysis out(test);
		std::cout << "Analysis starting" << times + 1 << std::endl;
		int Cont=out.ThreadRun();
		if (!Cont)
		{
			continue;
		}
		dense_vector<GaPMatrix*> p = out.Result;

		for (int i = 0; i < size(p); i++ )
		{
		
			if (p[i]->Tmarker.P_value==-1)
			{
				break;
			}
			else 
			{
				result << p[i]->Tmarker.P_value << std::endl;
				result.flush();
				times++;
		/*		for (int i = 0; i < size(p); i++)
				{
					result << "SNP" << i << "\t";
					result << p[i]->Tmarker.LR << "\t";
					result << p[i]->Tmarker.P_value << "\t";
					result << p[i]->Tmarker.b_value << "\t";
					result << p[i]->Tmarker.it << "\t";
					result << "\n";
					result.flush();

				}
				times++;*/
			}
			
		}

	}
	
	result.close();
	clock_t endtime = clock();
	std::cout << "Finished" << std::endl;
	std::cout << "Elapsed time : " << (endtime - starttime) / CLOCKS_PER_SEC << " S" << std::endl;
	

	return 0;
}

void Trans2PED(Simulate sm)
{
	std::ofstream f("SM.lemped");
	int k = 1;
	for (int i=0;i<sm.nInds;i++)
	{
		
		f << i + 1 << "\t" << k++ << "\t" << k << "\t" << k+1 << "\t" << 0 << "\t" << sm.sim[i].phe;
		f.flush();
		for (int j = 0; j < sm.nSNPs; j++)
		{
			f << "\t" << sm.sim[i].ChildGen[j];
		}
		f << "\n";
		f << i + 1 << "\t" << k++ << "\t" << 0 << "\t" << 0 << "\t" << 1 << "\t" << 0;//male
		for (int j = 0; j < sm.nSNPs; j++)
		{
			f << "\t" << sm.sim[i].ParentGen[j][0];
		}
		f << "\n";
		f << i + 1 << "\t" << k++ << "\t" << 0 << "\t" << 0 << "\t" << 2 << "\t" << 0;//female
		for (int j = 0; j < sm.nSNPs; j++)
		{
			f << "\t" << sm.sim[i].ParentGen[j][1];
		}
		f << "\n";
		f.flush();
	}
	f.close();
	std::cout << "Done" << std::endl;
}