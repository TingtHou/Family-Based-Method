#pragma once
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <time.h>
class RandomDouble
{
public:
	RandomDouble();
	boost::mt19937 *engine;

public:
	double generation01();
	double generation02();
};

