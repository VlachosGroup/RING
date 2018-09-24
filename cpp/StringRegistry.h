#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <deque>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

class StringRegistry {

private:
	static map<string, string*> registry;
	

public:
	static bool InRegistry(const string &);
	static void InsertIntoRegistry(const string &);
	static string* getStringPointer(const string &);
	static void RemoveFromRegistry(const string &);

};

class CompositeAtomsRegistry
{
	protected:
		static vector<string> CompositeAtomsList;
	public:
		static void InsertIntoList(string);
		static int getIndexOfAtom(string);
};


