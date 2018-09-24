#ifndef STRINGREG_H
#define STRINGREG_H

#include <string>
#include <vector>
#include <map>

class StringRegistry {
private:
	static std::map<std::string, std::string*> registry;
public:
	static bool InRegistry(const std::string &);
	static void InsertIntoRegistry(const std::string &);
	static std::string* getStringPointer(const std::string &);
	static void RemoveFromRegistry(const std::string &);
};

class CompositeAtomsRegistry
{
	protected:
		static std::vector<std::string> CompositeAtomsList;
	public:
		static void InsertIntoList(std::string);
		static int getIndexOfAtom(std::string);
};

#endif
