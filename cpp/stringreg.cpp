#include <string>
#include <vector>
#include <utility>
#include <map>

#include "stringreg.h"

using std::vector; using std::map; using std::string; using std::pair;

map<string, string*> StringRegistry::registry;
vector<string> CompositeAtomsRegistry::CompositeAtomsList;

bool StringRegistry::InRegistry(const std::string &s)
{
	if (registry.count(s)>0)
        return true;
	else 
        return false;
}

void StringRegistry::InsertIntoRegistry(const std::string & s)
{
	if (registry.count(s)==0)
	{
		string * sptr = new string (s);
		registry.insert(pair<string, string*>(*sptr,sptr));
	}
}

string* StringRegistry::getStringPointer(const std::string & s)
{
	if (registry.count(s)== 0)
	{
		InsertIntoRegistry(s);
	}
	string * result;
	result = registry[s];
	return result;
}

void StringRegistry::RemoveFromRegistry(const std::string & s)
{
	if (registry.count(s)>0)
	{
		string * molptr = registry[s];
		delete molptr;
		registry.erase(registry.find(s));
	}
}

void CompositeAtomsRegistry::InsertIntoList(string S)
{
	CompositeAtomsList.push_back(S);
}

int CompositeAtomsRegistry::getIndexOfAtom(string S)
{
	for (int i=0;i<CompositeAtomsList.size();i++)
		if (CompositeAtomsList[i]==S)
			return i;
	return -1;
}
