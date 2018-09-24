#ifndef ADDITIONALFUNC_H
#define ADDITIONALFUNC_H

#include <iostream>
#include <fstream>
#include <string>
//#include <cstring>
#include <sstream>
#include <vector>
//#include <list>
#include <set>
//#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
#include "additionalfunc.h"


using namespace std;



// initializing the array of prime numbers to be used in generation of unique smiles;
int prime(int i)
{
	int primenumbers[100]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
				73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,
				173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,
				271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,
				383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
				491,499,503,509,521,523,541};

	return primenumbers[i];
}

int factorial (int a)
{
  if (a > 1)
   return (a * factorial (a-1));
  else
   return (1);
}

int moleculesize(string a)
{
	int count=0;
	bool InsideCurlyBrackets = false;
	bool InsideSquareBrackets = false;
	for (unsigned int i=0; i<a.length(); i++)
	{
		string c;
		c=a[i];
		
		if (c=="[")
		{
			count++;
			InsideSquareBrackets = true;
			
		}
		if (c=="]")
		{
			InsideSquareBrackets = false;
			
		}

		if (c=="{")
		{
			if (!InsideSquareBrackets)count++;
			InsideCurlyBrackets = true;
			
		}
		if (c=="}")
		{
			InsideCurlyBrackets = false;
			
		}
		if ((isalpha(a[i])|| (c=="$") || (c=="&") )&& !(InsideCurlyBrackets) && !(InsideSquareBrackets))
		{
			count++;
		}
	}
	return count;
}

int patternsize(string a)
{
	int count =0;
	int flag=0;
	bool InsideCurlyBrackets = false;
	int InsideSquareBrackets = 0;
	bool InsideSingleQuotes = false;
	unsigned int i=0;
	while (i<a.length())
	{
		string c;
		c=a[i];
				
		if(c=="[")
		{
			InsideSquareBrackets++;
			//cout<<"inside sq br"<<endl;
		}
		if (c=="]")
		{		
			InsideSquareBrackets--;
			//cout<<"outside sq br"<<endl;
		}
		if (c=="'")
		{
			if (!InsideSingleQuotes && InsideSquareBrackets==0)
			{
				//cout<<"counting inside single quotes"<<endl;
				count++;
				InsideSingleQuotes = true;
			//	cout<<"inside single"<<endl;
			}
			else 
			{InsideSingleQuotes = false; }//cout<<"outside single"<<endl;}
		}
		if (c=="{")
		{
			if (InsideSquareBrackets==0 && !InsideSingleQuotes)count++;
			//cout<<"Inside cu br"<<endl;
			
			while (c!="}")
			{
				i++;
				c=a[i];
				
			}
			
			//cout<<"outside cu br"<<endl;
		}
		if ((isalpha(a[i]) && (c!="H") && (c!="R"))  && (InsideSquareBrackets==0) && (!InsideSingleQuotes))
		{
			//cout<<"yes"<<endl;
			count++;
		}
		if (count==0 && c=="H" && !InsideSingleQuotes && InsideSquareBrackets==0)count++;

		if (((c=="$") || (c=="&"))&& (InsideSquareBrackets==0) && (!InsideSingleQuotes))
		{
			count++;
			//cout<<"yes 2"<<endl;
		}

		i++;

	}
	return count;
}


void throwsmileserror(string S)
{
	if (S.compare(0,1,")")==0)throw pair<string,int>(S,0);
	int opencount=0;
	int closecount=0;
	set<string> AcceptedSymbols;
	AcceptedSymbols.insert("C");
	AcceptedSymbols.insert("c");
	AcceptedSymbols.insert("O");
	AcceptedSymbols.insert("o");
	AcceptedSymbols.insert("N");
	AcceptedSymbols.insert("n");
	AcceptedSymbols.insert("P");
	AcceptedSymbols.insert("p");
	AcceptedSymbols.insert("S");
	AcceptedSymbols.insert("s");
	AcceptedSymbols.insert("H");
	AcceptedSymbols.insert("D");
	AcceptedSymbols.insert("{");
	AcceptedSymbols.insert("}");
	AcceptedSymbols.insert("(");
	AcceptedSymbols.insert(")");
	AcceptedSymbols.insert("[");
	AcceptedSymbols.insert("]");
	AcceptedSymbols.insert(".");
	AcceptedSymbols.insert(":");
	AcceptedSymbols.insert("\"");
	AcceptedSymbols.insert("'");
	AcceptedSymbols.insert("*");
	AcceptedSymbols.insert("+");
	AcceptedSymbols.insert("=");
	AcceptedSymbols.insert("-");
	AcceptedSymbols.insert("_");
	AcceptedSymbols.insert("#");



	bool InsideSquareBraces = false;
	bool InsideCurlyBraces = false;
	for (unsigned int i=0;i<S.length();i++)
	{
		if (S.compare(i,1,"(")==0)
		{
			opencount++;
		}
		if (S.compare(i,1,")")==0)
		{
			closecount++;
			if (S.compare(i-1,1,"(")==0) throw pair<string,int>(S,1);
		}
		if (S.compare(i,1,"[")==0)
		{
			InsideSquareBraces = true;
		}
		if (S.compare(i,1,"]")==0)
		{
			if (!InsideSquareBraces) throw pair<string,int>(S,7);
			InsideSquareBraces = false;
			if (S.compare(i-1,1,"[")==0) throw pair<string,int>(S,2);
		}
		if (S.compare(i,1,"{")==0) InsideCurlyBraces = true;
		if (S.compare(i,1,"}")==0)
		{
			if (!InsideCurlyBraces) throw pair<string,int>(S,8);
			InsideCurlyBraces = false;
			if (S.compare(i-1,1,"{")==0) throw pair<string,int>(S,3);
		}
		string c;
		c=S[i];
		if (AcceptedSymbols.count(c)==0 && !InsideCurlyBraces && !isdigit(S[i]))throw pair<string,string>(S,c);
	}
	if (opencount!=closecount)throw pair<string,int>(S,4);
	if (InsideCurlyBraces)throw pair<string,int>(S,5);
	if (InsideSquareBraces) throw pair<string,int>(S,6);
}

void catchsmileserror(pair<string,int> error)
{
	if (error.second==0)cout<<"SMILES string "<<error.first<<" cannot start with a ')'"<<endl;
	if (error.second==1)cout<<"SMILES "<<error.first<<" cannot have nothing between braces'()'"<<endl;
	if (error.second==2)cout<<"SMILES "<<error.first<<" cannot have nothing between square braces'[]'"<<endl;
	if (error.second==3)cout<<"SMILES "<<error.first<<" cannot have nothing between curly braces'{}'"<<endl;
	if (error.second==4)cout<<"SMILES "<<error.first<<" has an open and close braces mismatch"<<endl;
	if (error.second==5)cout<<"SMILES "<<error.first<<" has a missing } brackets"<<endl;
	if (error.second==6)cout<<"SMILES "<<error.first<<" has a missing ]"<<endl;
	if (error.second==7)cout<<"SMILES "<<error.first<<" has a ] when it doesnt have a corresponding ["<<endl;
	if (error.second==8)cout<<"SMILES "<<error.first<<" has a } when it doesnt have a corresponding {"<<endl;
}



bool TokenizeIntoWords(string& S, string::size_type position, vector<string>& tokens)
{
	
	
	string::size_type position2;
	position= S.find_first_not_of(" ",position);
	if (position==string::npos)
		return true;
	else
	{
		position2 = S.find_first_of(" ",position);
		if (position2!=string::npos)
		{
		
			tokens.push_back(S.substr(position, position2-position));
			return TokenizeIntoWords(S,position2, tokens);
		}
		else 
		{	
			tokens.push_back(S.substr(position, S.length()-position));
			return true;
		}
	}
}


double StringToDouble(string S)
{
	stringstream ss(S);
	double DoubleValue;
	ss>>DoubleValue;

	if (ss.fail())
	{
		cout<<"could not convert "<<S<<"into double!"<<endl;
		throw 1;
	}

	return DoubleValue;
}

int StringToInt(string S)
{
	stringstream ss(S);
	int IntValue;
	ss>>IntValue;

	if (ss.fail())
	{
		cout<<"could not convert "<<S<<"into double!"<<endl;
		throw 1;
	}

	return IntValue;
}
string DoubleToString (double d)
{
	stringstream ss;
	ss << d;
	string str = ss.str();
	return str;
}

string IntToStr (int i)
{
	stringstream ss;
	ss << i;
	return  ss.str();
}

void WriteInfoFile(string str, ofstream& file)
{
	file<<str;
}

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
