#pragma once
#include <boost/numeric/mtl/mtl.hpp>
#include "RandomDouble.h"
#include "SimPar.h"
using namespace mtl;
#define COL_VECTOR mtl::vec::parameters<mtl::tag::col_major>
#define ROW_VECTOR mtl::vec::parameters<mtl::tag::row_major>
class Subject
{
public:
	Subject(SimPar par);
	Subject();
	void Generation();//Generating  the genotypes and phenotype of a individual;
public:
	dense_vector<dense_vector<double>> ParentGen; //genotypes of the parents, whose sizes equals to ParentGen size
	dense_vector<double> ChildGen;//genotypes of the child of the parent
	double phe; //the phenotype of the child, 1 means having disease, others not
protected:
	
	void GenerationOfGEN();//Generating  the genotypes
	void GenerationOfCausalGEN();//Generating  the genotypes
	void GenerationOfPHE();//Generating  the phenotype
	dense_vector<double> GenerationParentGEN(int i);//Generating  the genotypes of the Parents
	double GenerationChildGEN(dense_vector<double> ParentAllele,int i);//Generating  the genotypes of the child
	dense_vector<double,ROW_VECTOR> Paremter();//Translate the marker type to the coefficient of the model.
protected:
	double MAFofM;//the MAF of the link SNP
	double MAFofQ;//the MAF of the disease allele
	int nSNPs; //the number of the SNPs
	int ninds; //the number of the individuals
	double Theta;//the combination rate between the link SNP and the disease alleles
	double Delta;//the association rate between the link SNP and the disease alleles
	dense_vector<std::string> SNPs; //the SNPs which affect the disease
	double GenEffect;//genetic effect under the additive model
	int loci;//the location of the linage snp;
	double m = 0;//average gene effect
	double d = 0;//domain gene effect
	dense_vector<double> CausalofChild;
	dense_vector<double> CausalofM;
	dense_vector<double> CausalofP;
	dense_vector<dense_vector<double>> ParentSNPGene;
	dense_vector<double,COL_VECTOR> var;// the variable of the model
	
};

