#pragma once
#include "stdafx.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#define COL_VECTOR mtl::vec::parameters<mtl::tag::col_major>
#define ROW_VECTOR mtl::vec::parameters<mtl::tag::row_major>
using namespace mtl;
class SimPar
{

protected:
	double MAFofM;//the MAF of the link SNP
	double MAFofQ;//the MAF of the disease allele
	int nSNPs; //the number of the SNPs
	int ninds; //the number of the individuals
	double Theta;//the combination rate between the link SNP and the disease alleles
	double Delte;//the association rate between the link SNP and the disease alleles
	double GenEffect;//genetic effect under the additive model
	int SNPs; //the SNPs which affect the disease
	double d;
	double m;
public:
	SimPar();
	void Setd(double d) ;
	void Setm(double m);
	double getd() { return d; };
	double getm() { return m; };
	void SetMAFofM(double MAF);
	void SetMAFofQ(double MAF);
	void SetNumofSNP(int SNPNUM);
	void SetNumofInd(int INDNUM);
	void SetTheta(double T);
	void SetDelta(double D);
	void SetLinkSNP(int snps);
	void SetGenEffect(double GeneEffect);
	double GetMAFofM();
	double GetMAFofQ();
	int GetNumofSNP();
	int GetNumofInd();
	double GetTheta();
	double GetDelta();
	double GetGenEffect();
	int GetLinkSNP();
	dense_vector<double,COL_VECTOR> GetVar();
	std::string toString();
	void getu();
	void CalcParemeter();
	void getalpha();
	void getbeta();
	double u;
	double alpha;
	double beta;
	double pi;
	dense_vector<double,COL_VECTOR> var;

};

