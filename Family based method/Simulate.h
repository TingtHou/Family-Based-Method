#pragma once
#include "SimPar.h"
#include "Subject.h"
#include <boost/numeric/mtl/mtl.hpp>
using namespace mtl;
class Simulate
{
public:
	Simulate(SimPar parameter);
	Simulate();
protected:
    SimPar par;



public:
	void Generation();

	void printout();
public:
	dense_vector<Subject> sim;
	int nSNPs;
	int nInds;
};

