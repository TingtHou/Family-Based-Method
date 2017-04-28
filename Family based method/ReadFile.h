#pragma once
#include <string>
#include <fstream>
#include "Simulate.h"
#include "Subject.h"
#include "boost/algorithm/string/split.hpp"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>
using namespace std;

class ReadFile
{
public:
	ReadFile(string fileaddress);
	~ReadFile();
	Simulate getSR();
protected:
	Simulate SR;
	int linecount = 0;
	int columnscount = 0;
	string address;
	void CalcRowColCount();
	void Start2Read();
};

