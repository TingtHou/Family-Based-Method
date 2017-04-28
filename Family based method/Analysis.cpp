/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////  This class is used to initial the thread pool and start the multi-thread calculation 
//////////   After Analysis class initialization, just use ThreadRun to start.
//////////
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Analysis.h"
#include "SearchGene.h"

#include <boost/thread/mutex.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/bind.hpp>
#include "ThreadPool.h"

Analysis::Analysis(Simulate SimResult)
{
	SR = SimResult;
	Result=dense_vector<GaPMatrix*>(SR.nSNPs);
}

int Analysis::ThreadRun()
{
	using nbsdx::concurrent::ThreadPool;
//	ThreadPool<10> tp; //set max thread

	Procedure *bar=new Procedure; //the a procedure bar
	bar->nowProce = 0;
	bar->TotalP = SR.nSNPs;
//	dense_vector<std::future<int> > Procedure(SR.nSNPs);
	int SignalLRT=1;
	for (int i=0;i<SR.nSNPs;i++)//create thread to analysis
	{
		Result[i] = new GaPMatrix;
		GetYandX(i, Result[i]);
	//	SearchGene T1Search(OneT);
	//	SearchGene::StartAnalysis(Result[i]);
	//	tp.AddJob(boost::bind(SearchGene::StartAnalysis, Result[i],bar,true));
		SearchGene::StartAnalysis(Result[i],bar,true);
	}
	//tp.wait();//waiting all thread done
	//tp.JoinAll();
	return SignalLRT;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////Store the covariate matrix, phenotype vector and gene data into GaPMatrix///////////////
void Analysis::GetYandX(int Gene,GaPMatrix* ma)
{
	
	ma->InData.phe = dense_vector<double>(SR.nInds);
	ma->InData.Xinfo = dense2D<double>(SR.nInds, 6);
	ma->InData.triadsGene = dense2D<double>(SR.nInds, 3);
	ma->InData.X_LRT = dense2D<double>(SR.nInds, 3);
	ma->Tmarker.GenID = Gene;
	dense2D<double> X(SR.nInds, 6);
	dense_vector<double> Y(SR.nInds);
	//initial the X matrix, and Y vector
	for (int IndID = 0; IndID < SR.nInds; IndID++)
	{
		dense_vector<double> indX = TranslateX(SR.sim[IndID], Gene);
	
		for (int i = 0; i < size(indX); i++)
		{
			ma->InData.Xinfo[IndID][i] = indX[i];
	
		}
		ma->InData.phe[IndID] = (double)SR.sim[IndID].phe;
		//store gene data into triadsGene which will be used in regression
		ma->InData.triadsGene[IndID][0] = SR.sim[IndID].ParentGen[Gene][0];
		ma->InData.triadsGene[IndID][1] = SR.sim[IndID].ParentGen[Gene][1];
		ma->InData.triadsGene[IndID][2] = SR.sim[IndID].ChildGen[Gene];

	}
	//store covariate matrix which is under null hypothesis θ=0.5
	for (int IndID = 0; IndID < SR.nInds; IndID++)
	{
		
		dense_vector<double> indX_LRT = TranslateXC(SR.sim[IndID], Gene);
		for (int i = 0; i < size(indX_LRT); i++)
		{
		
			ma->InData.X_LRT[IndID][i] = indX_LRT[i];
		}
	
	}
	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Get matrix of X under the alternative hypothesis θ≠0.5
dense_vector<double> Analysis::TranslateX(Subject SubInd, int loci)
{
	dense_vector<double> para(6);
	//dense_vector<int, ROW_VECTOR> para(6);
	set_to_zero(para);
	para[0] = 1;
	int sumP = SubInd.ParentGen[loci][0] + SubInd.ParentGen[loci][1];
	int sumA = SubInd.ParentGen[loci][0] + SubInd.ParentGen[loci][1] + SubInd.ChildGen[loci];
	para[1] = 2 * (SubInd.ChildGen[loci] - 1);
	para[2] = pow(-1, SubInd.ChildGen[loci] + 1);
	if ((SubInd.ParentGen[loci][0] == SubInd.ParentGen[loci][1]) && (SubInd.ParentGen[loci][0] == 1))
	{


		para[3] = 4 * (1 - SubInd.ChildGen[loci]);
		para[4] = 4 * pow(-1, SubInd.ChildGen[loci]);
		para[5] = 4 * pow(-1, SubInd.ChildGen[loci] + 1);

	}
	else
	{
		para[3] = sumA % 3 ? 2 * pow(-1, (sumA % 3) + 1) : 0;
		para[4] = sumA % 3 ? (sumA % 2 ? 2 : -2) : 0;
		para[5] = 0;

	}

	return para;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Get matrix of X under the alternative hypothesis θ=0.5
dense_vector<double> Analysis::TranslateXC(Subject SubInd, int loci)
{

	dense_vector<double> para(3);
	set_to_zero(para);
	para[0] = 1;
	int sumP = SubInd.ParentGen[loci][0] + SubInd.ParentGen[loci][1];
	int sumA = SubInd.ParentGen[loci][0] + SubInd.ParentGen[loci][1] + SubInd.ChildGen[loci];
	para[1] = sumP - 2;
	switch (sumP)
	{
	case 4:
		para[2] = -1;
		break;
	case 3:
		para[2] = 0;
		break;
	case 2:
		if (SubInd.ParentGen[loci][0] == SubInd.ParentGen[loci][1])
		{
			para[2] = 0;
		}
		else
		{
			para[2] = 1;
		}
		break;
	case 1:
		para[2] = 0;
		break;
	case 0:
		para[2] = -1;
		break;
	}
	return para;
}
