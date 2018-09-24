#ifndef ADDITIONALFUNC_H
#define ADDITIONALFUNC_H

#include <fstream>
#include <string>
#include <vector>
//#include <list>
//#include <set>
//#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>



// initializing the array of prime numbers to be used in generation of unique smiles;
int prime(int i);

int factorial (int a);

int moleculesize(std::string a);

int patternsize(std::string a);

void throwsmileserror(std::string S);

void catchsmileserror(std::pair<std::string,int> error);

bool TokenizeIntoWords(std::string& S, std::string::size_type position, std::vector<std::string>& tokens);

double StringToDouble(std::string S);

int StringToInt(std::string S);

std::string DoubleToString (double d);

std::string IntToStr (int i);

void WriteInfoFile(std::string str, std::ofstream& file);

/*static 	vector<double> TemperatureRangePointsForGA()
{
	vector<double> TempRange;
	TempRange.push_back(298.);
	TempRange.push_back(400.);
	TempRange.push_back(500.);
	TempRange.push_back(600.);
	TempRange.push_back(800.);
	TempRange.push_back(1000.);
	TempRange.push_back(1500.);
	return TempRange;
}*/

#endif
