////////////////////////////////////////////////////////////////////////////////////////////////////
////////  The class is used to set some parameters to simulation
///////   this must to set all parameters, which have been declared in the SimPar.h
///////   then using the function CalcParemeter(), it will caluate the 6 variables in our model 
///////
///////
///////
//////////////////////////////////////////////////////////////////////////////////////////////////


#include "SimPar.h"
#include <sstream>


SimPar::SimPar()
{
	
}

void SimPar::Setd(double d)
{
	this->d = abs(d) < Precision_Zero ? 0 : d;
}

void SimPar::Setm(double m)
{
	this->m = abs(m) < Precision_Zero ? 0: m;
}

void SimPar::SetMAFofM(double MAF)
{
	
	MAFofM =MAF;
}

void SimPar::SetMAFofQ(double MAF)
{
	MAFofQ = MAF;
}

void SimPar::SetNumofSNP(int SNPNUM)
{
	nSNPs = SNPNUM;
}

void SimPar::SetNumofInd(int INDNUM)
{
	ninds = INDNUM;
}

void SimPar::SetTheta(double T)
{
	Theta = T;

}

void SimPar::SetDelta(double D)
{
	Delte = D;
}

void SimPar::SetLinkSNP(int snps)
{
	SNPs = snps;
	/*using namespace std;
	
	list<boost::iterator_range<string::iterator>> l;

	boost::split(SNPs, snps, boost::is_any_of(",.-+"));
	*/
	/*
	for (string x:SNPs)
	{
		
	//	SNPs.insert_element(i++, x);
		cout << x << "\n" << endl;
	}*/
}

void SimPar::SetGenEffect(double GeneEffect)
{
	GenEffect = GeneEffect;
}

double SimPar::GetMAFofM()
{
	return MAFofM;
}

double SimPar::GetMAFofQ()
{
	return MAFofQ;
}

int SimPar::GetNumofSNP()
{
	return nSNPs;
}

int SimPar::GetNumofInd()
{
	return ninds;
}

double SimPar::GetTheta()
{
	return Theta;
}

double SimPar::GetDelta()
{
	return Delte;
}

double SimPar::GetGenEffect()
{
	return GenEffect;
}

int SimPar::GetLinkSNP()
{
	return SNPs;
}

dense_vector<double,COL_VECTOR> SimPar::GetVar()
{
	return var;
}

std::string SimPar::toString()
{
	
	std::ostringstream outline;
	outline << "The number of the individuals : " << ninds<<"\n";
	outline << "The number of the SNPs : " << nSNPs << "\n";
	outline << "The linkage loci : " << SNPs << "\n";
	outline<<"Gene Effect : "<< GenEffect<< "\n";
	outline << "The 4 Parameters : u " << u << "\talpha : " << alpha << "\tbeta : " << beta << "\ttheta : " << Theta << "\n";
	outline << "The 6 Parameters : u " << u << "\talpha : " << alpha << "\tbeta : " << beta << "\ttheta*alpha : " << var[3] << "\ttheta*beta : " << var[4] << "\ttheta^2*beta : " << var[5] << "\n";
	
	return outline.str();
}
void SimPar::getu()
{
	double PQ = 1 - 2 * MAFofQ;
	double PMTH = (1-2*MAFofM)*Delte / ((1 - MAFofM)*MAFofM);
	double test1 = (PQ - PMTH)*(GenEffect);
	double test2 = (1 - (PQ - PMTH)*(PQ - PMTH)*(PQ - PMTH)*(PQ - PMTH))*d / 2;
	u = m+ test1 + test2;
	std::cout << u << std::endl;
//	u = (PQ - PMTH)* GenEffect;
}

void SimPar::CalcParemeter()
{
	getu();
	getbeta();
	getalpha();
//	pi = 1 - 2 * Theta;
	var=dense_vector<double,COL_VECTOR>(6);
	var[0]= u;
	var[1]= alpha;
	var[2]= beta;
	var[3]= Theta*alpha;
	var[4]= Theta*beta;
	var[5]= pow(Theta, 2)*beta;

}

void SimPar::getalpha()
{
	double PQ = 1-2 * MAFofQ;
	double PMTH = (1-2 * MAFofM)*Delte / ((1 - MAFofM)*MAFofM);
	alpha = Delte*((GenEffect) - (PQ - PMTH)*d) / (2 * (1 - MAFofM)*MAFofM);
//	alpha = Delte*GenEffect  / (2 * (1 - MAFofM)*MAFofM);
}

void SimPar::getbeta()
{
	beta = pow((Delte / ((1 - MAFofM)*MAFofM)), 2)*d / 2;
//	beta = 0;
}
