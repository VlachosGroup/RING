#ifndef CLASS_HEADER_H
#define CLASS_HEADER_H

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
using namespace std;

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
	bool operator() (const pair<unsigned int, unsigned int>& lhs, const pair<unsigned int,unsigned int>& rhs)
	{
		if (lhs.first<rhs.first) return true;
		else if (rhs.first==lhs.first)
			return lhs.second<rhs.second;
		else return false;
	}
};

class CHEMKinFiles
{
	protected:
		rxn_net_gen* Network;
		set<int> TempsForCpCalc;
		vector<KineticParamPtr>* kinetics;
		double siteDensity;
		double Rg;
		map<string,double> AllSites;
		multimap<string,string> DenticityInfoMap;
		set<string> SitesRequiringDenticity;

		
		int CreateCHEMKINfileAndWriteSpecies();
		vector<string> CHEMKINRxnStrings;
		vector<string> CHEMKINPreExpInfo;
		vector<string> CHEMKINIndexInfo;
		vector<string> CHEMKINActEInfo;
		set<int> RxnsWithStickingCoeff;

		string generateCHEMKINRxnStrings(int);
		void generateCHEMKINkineticsStrings(double, double, double, bool);
		int WriteRxnInfoInCHEMKINfile();
		map<string*,string, classcomp> CHEMKINSpeciesName;
		map<string*,string, classcomp> CHEMKINBulkName;
		set<string> AllElements;
		void FindAllElements();
		bool isGasPhaseRxn(int);
		bool SetKinetics();
		multimap<string, pair<string*,int> > SiteOccupancyMap;//this just checks for composite heterogeneous atoms
		multimap<string, pair<string*,int> > HeterogeneousSitesOccupancyMap;//this checks for the actually defined heterogeneous sites -note the small but important difference - a heterogeneous site will require a heterogeneous composite site atom but can include more. for e.g. Zeo could be a heterogeneous composite site atom, but the actual free heterogeneous site is [{Zeo}H]
		bool PopulateSiteOccupancyMap(string*, string);
		void PopulateHeterogeneousSiteOccupanyMap(string*);
		set<string*> surfaceSpecies;
		int calculateNetSiteCount(string,int);

	public:
		CHEMKinFiles(rxn_net_gen*, set<int>&,vector<KineticParamPtr>&, double, map<string,double>&,multimap<string,string>&);
		void generateFiles();
};

class GAMSFiles
{
	protected:
		rxn_net_gen* Network;
		char* filename;
		vector<KineticParamPtr>* kineticsFns;
		set<string> AllSites;
		multimap<string,string> DenticityInfoMap;
		set<string> SitesRequiringDenticity;
		set<string*> surfaceSpecies;
		bool PopulateSiteOccupancyMap(string*, string);
		multimap<string, pair<string*,int> > SiteOccupancyMap;
		map<string, double> SiteDensity; //assuming right now that only one type of site exists -- need to allow for more
	public:
		GAMSFiles(rxn_net_gen*, vector<KineticParamPtr>&, set<string>&, map<string, double>&, multimap<string,string>&, char*);
		bool generateFiles();
};

class ThermoGroupsIdentification
{
	
	public:
		ThermoGroupsIdentification(vector<pair<string, string> >&, vector<pair<string, string> >&, char*,char*);
};


#endif
