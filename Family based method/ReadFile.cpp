#include "ReadFile.h"
#include <iostream>
#include <iomanip>
using namespace mtl;

ReadFile::ReadFile(string fileaddress)
{
	address = fileaddress;
	CalcRowColCount();
	Start2Read();
}

ReadFile::~ReadFile()
{
}

Simulate ReadFile::getSR()
{
	return SR;
}
void ReadFile::CalcRowColCount()
{
	ifstream in(address);
	
	if (!in.is_open())
	{
		cout << "ERROR, Can not open the file" << endl;
		return;
	}
	string tmp;
	getline(in, tmp);
	
	dense2D<string> tmplist;
	using namespace boost;
	vector<string> listnm;
	split(listnm, tmp,boost::is_any_of("\t"), boost::token_compress_on);
	columnscount = listnm.size()-1;
	while (!in.eof())
	{
		getline(in, tmp);
		linecount++;
	}
	in.close();
}

void ReadFile::Start2Read()
{
	SR.nInds = linecount;
	SR.nSNPs = columnscount;

	SR.sim = dense_vector<Subject>(linecount);

	ifstream in(address);

	if (!in.is_open())
	{
		cout << "ERROR, Can not open the file" << endl;
		return;
	}
	string tmp;
	getline(in, tmp);
	int c = 0;

	while (!in.eof())
	{
		string lineone;
		getline(in, lineone);
		vector<string> listnm;
		split(listnm, lineone, boost::is_any_of("\t"), boost::token_compress_on);
		Subject ind;
		ind.ParentGen = dense_vector<dense_vector<double>>(columnscount);
		ind.ChildGen = dense_vector<double>(columnscount);
		ind.phe = stod(listnm.at(0));
		for (int id=1;id<listnm.size();id++)
		{
			
			vector<string> genlist;
			split(genlist, listnm.at(id), boost::is_any_of(" "), boost::token_compress_on);
			dense_vector<double> parentGen(2);
			parentGen[0] = stod(genlist.at(0));
			parentGen[1] = stod(genlist.at(1));
			ind.ParentGen[id - 1] = parentGen;
			ind.ChildGen[id - 1] = stod(genlist.at(2));

		}
		SR.sim[c++] = ind;
	}
	in.close();
}
