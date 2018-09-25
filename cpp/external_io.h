#ifndef EXTERNAL_IO_H
#define EXTERNLA_IO_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>

#include "kinetics.h"
#include "rng.h"

class CHEMKinFiles
{
	protected:
		rxn_net_gen* Network;
		std::set<int> TempsForCpCalc;
		std::vector<KineticParamPtr>* kinetics;
		double siteDensity;
		double Rg;
        std::map<std::string,double> AllSites;
		std::multimap<std::string,std::string> DenticityInfoMap;
		std::set<std::string> SitesRequiringDenticity;

		
		int CreateCHEMKINfileAndWriteSpecies();
		std::vector<std::string> CHEMKINRxnStrings;
		std::vector<std::string> CHEMKINPreExpInfo;
		std::vector<std::string> CHEMKINIndexInfo;
		std::vector<std::string> CHEMKINActEInfo;
		std::set<int> RxnsWithStickingCoeff;

		std::string generateCHEMKINRxnStrings(int);
		void generateCHEMKINkineticsStrings(double, double, double, bool);
		int WriteRxnInfoInCHEMKINfile();
        std::map<std::string*,std::string, classcomp> CHEMKINSpeciesName;
        std::map<std::string*,std::string, classcomp> CHEMKINBulkName;
		std::set<std::string> AllElements;
		void FindAllElements();
		bool isGasPhaseRxn(int);
		bool SetKinetics();
		std::multimap<std::string, std::pair<std::string*,int> > SiteOccupancyMap;//this just checks for composite heterogeneous atoms
		std::multimap<std::string, std::pair<std::string*,int> > HeterogeneousSitesOccupancyMap;//this checks for the actually defined heterogeneous sites -note the small but important difference - a heterogeneous site will require a heterogeneous composite site atom but can include more. for e.g. Zeo could be a heterogeneous composite site atom, but the actual free heterogeneous site is [{Zeo}H]
		bool PopulateSiteOccupancyMap(std::string*, std::string);
		void PopulateHeterogeneousSiteOccupanyMap(std::string*);
		std::set<std::string*> surfaceSpecies;
		int calculateNetSiteCount(std::string,int);

	public:
		CHEMKinFiles(rxn_net_gen*, std::set<int>&,std::vector<KineticParamPtr>&, double, std::map<std::string,double>&,std::multimap<std::string,std::string>&);
		void generateFiles();
};

class GAMSFiles
{
	protected:
		rxn_net_gen* Network;
		char* filename;
		std::vector<KineticParamPtr>* kineticsFns;
		std::set<std::string> AllSites;
		std::multimap<std::string,std::string> DenticityInfoMap;
		std::set<std::string> SitesRequiringDenticity;
		std::set<std::string*> surfaceSpecies;
		bool PopulateSiteOccupancyMap(std::string*, std::string);
		std::multimap<std::string, std::pair<std::string*,int> > SiteOccupancyMap;
        std::map<std::string, double> SiteDensity; //assuming right now that only one type of site exists -- need to allow for more
	public:
		GAMSFiles(rxn_net_gen*, std::vector<KineticParamPtr>&, std::set<std::string>&, std::map<std::string, double>&, std::multimap<std::string,std::string>&, char*);
		bool generateFiles();
};

class ThermoGroupsIdentification
{
	
	public:
		ThermoGroupsIdentification(std::vector<std::pair<std::string, std::string> >&, std::vector<std::pair<std::string, std::string> >&, char*,char*);
};


#endif
