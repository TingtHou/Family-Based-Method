//////////////////////////////////////////////////////////////////////////
////////  The class is used to simulate data
///////   Using it like this:
///////  Simulate a(par);
///////  a.Generation();
///////
///////
////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"
#include "Simulate.h"
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include "RandomDouble.h"
#include <fstream>
Simulate::Simulate(SimPar parameter)
{
	par = parameter;
	sim=dense_vector<Subject>(par.GetNumofInd());
	nSNPs = par.GetNumofSNP();
	nInds=par.GetNumofInd();
}

Simulate::Simulate()
{
}
//This function is used to generate the data including phenotype and genotypes of  par.GetNumofInd()  size individuals
void Simulate::Generation()
{
	
	int indnum = par.GetNumofInd();
	for (int i=0;i<nInds;i++)
	{
		Subject newone(par);
		newone.Generation();
		sim[i] = newone;
	}
	/*
	int controlnum, casenum;
	casenum=controlnum= 0;
	int design = par.GetNumofInd()/2;
	int i = 0;
	while (casenum<design||controlnum<design)
	{
		Subject newone(par);
		newone.Generation();
		if (newone.phe == 1)
		{
			if (casenum < design)
			{
				sim[i++] = newone;
				casenum++;
			}
			
		}
		else
		{
			if (controlnum < design)
			{
				sim[i++] = newone;
				controlnum++;
			}
			
		}
	}*/
//	std::random_shuffle(sim.begin(), sim.end());
}

void Simulate::printout()
{

	std::ofstream  f;
	f.open("sim.txt");
	f << "Phenotype\t";
	for (int i=0;i<nSNPs;i++)
	{
		f << "SNP" << i << "\t";
	}
	f << "\n";
	for (int i=0;i<nInds;i++)
	{
		f << sim[i].phe << "\t";
		for (int j=0;j<nSNPs;j++)
		{
			f << sim[i].ParentGen[j][0] << " " << sim[i].ParentGen[j][1] << " " << sim[i].ChildGen[j] << "\t";
		}
		f << "\n";
	}
	f.close();
}


