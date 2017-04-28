//////////////////////////////////////////////////////////////////////////
////////  The class is used to generate a double random number
///////	  it is based on the boost
///////   you can use it like this:
///////           RandomDouble x;
///////	          x.generation(), then it returns a double random numder
///////
////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "RandomDouble.h"
#include <time.h>
#include <windows.h>
using namespace boost;

RandomDouble::RandomDouble()
{

	LARGE_INTEGER nFrequency;
	QueryPerformanceFrequency(&nFrequency);
	LARGE_INTEGER nStartCounter;
	QueryPerformanceCounter(&nStartCounter);

	engine = new boost::mt19937((unsigned)nStartCounter.LowPart);
	
	
}

double RandomDouble::generation01()
{
	boost::uniform_01<> *u01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<> > *die;
	u01 = new boost::uniform_01<>();
	die = new boost::variate_generator<boost::mt19937&, boost::uniform_01<> >(*engine, *u01);
	return (*die)();
}

double RandomDouble::generation02()
{
	
	boost::uniform_real<> *nm;
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > *die;
	nm = new boost::uniform_real<>(-1,1);
	die = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(*engine, *nm);
	return (*die)();
}