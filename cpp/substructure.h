#ifndef SUBSTRUCTURE_H
#define SUBSTRUCTURE_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <vector>
//#include <list>
//#include <set>
//#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;

#include "common.h"
//#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
#include "atomcontainer.h"
//#include "substructure.h"

class Substructure: public Atomcontainer//class Substructure derived from Atomcontainer
{
	protected:
		std::string fragmentstring; //contains the SMARTS like description of pattern
		std::vector <env_set>AtomEnv;//set of env_sets - each atom could have a set of atoms or groups and environmental constraints.
		std::vector <std::vector<int> > atom_flag;
		std::vector <Triplet> ringbondDefs;
	public:
		Substructure (std::string,int);//constructor
		void printstring();//prints out the fragment std::string
		std::string getstring();//gets the fragmentstd::string
		int isringbondcheck(int,int) const;//checks if a bond is in a ring: -1 implies no check reqd, 0 implies forbid ring bond and 1 implies have ring bond
		friend class Patternmatch;//Patternmatch is a friend class
		friend class rxn_net_gen;	
		friend class Reactiontype;
};

#endif
