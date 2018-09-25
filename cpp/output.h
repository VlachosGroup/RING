#ifndef OUTPUT_H
#define OUTPUT_H

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
#include <algorithm>
//using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "common.h"
#include "clonable.h"
#include "element.h"
#include "atom.h"
#include "atomcontainer.h"
#include "automorphs.h"
#include "substructure.h"
#include "patternmatch.h"
#include "kinetics.h"
#include "molecule.h"
#include "reaction.h"
#include "lumping.h"
#include "groupadditivity.h"
#include "generated_rxn.h"
#include "rng.h"
#include "qsar.h"
#include "pathways.h"
#include "mechanisms.h"

enum AstNodeType
{
	Undefined, OperatorAND, OperatorOR,OperatorNOT,Boolean
};

enum TokenType
{
	Error, AND, OR, NOT, EndOfText, Openbraces, Closedbraces, Size, Charge, IsAromatic, IsCyclic, MinRingSize, MaxRingSize, Fragment
};

struct paircomp {
	bool operator() (const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int,unsigned int>& rhs)
	{
		if (lhs.first<rhs.first) return true;
		else if (rhs.first==lhs.first)
			return lhs.second<rhs.second;
		else return false;
	}
};


#endif
