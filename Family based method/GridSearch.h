
#include "GStruction.h"
#include <boost/numeric/mtl/mtl.hpp>
using namespace mtl;
#define COL_VECTOR mtl::vec::parameters<mtl::tag::col_major>
#define ROW_VECTOR mtl::vec::parameters<mtl::tag::row_major>
class GridSearch
{
public:
	GridSearch();
	GridSearch(GSPar* par,bool LRT);
	dense_vector<double> getEstimate();
	double getLikelihood();
private:
	dense_vector<double> max;
	dense_vector<double> min;
	double likelihood = 0;
	double step = 3;
	double partion = 10;
	InputData Data;
	dense_vector<double,COL_VECTOR> E_Result;
	void StartEstimate( bool LRT);
	double CalcSocreSquare(dense_vector<double> &curValue,bool LRT);
	void GridSearch::Likelihood(dense_vector<double, COL_VECTOR> &curValue, bool LRT, double &likeli);
};

