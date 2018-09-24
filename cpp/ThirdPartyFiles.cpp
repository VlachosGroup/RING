#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <vector>
#include <cctype>
#include <queue>
#include <algorithm>
#include <set>
#include <cmath>
using namespace std;



#include "Classheader.h"
#include "AdditionalFunctions.h"
#include "StringRegistry.h"



CHEMKinFiles::CHEMKinFiles(rxn_net_gen* net, set<int>& cpTemps, vector<KineticParamPtr>& kineticFns, double GasConstant, map<string, double>& Sites, multimap<string,string>& Denticity)
{
	Network = net;
	TempsForCpCalc = cpTemps;
	kinetics = &kineticFns;
	Rg = GasConstant;
	AllSites = Sites;
	DenticityInfoMap = Denticity;

	for (multimap<string,string>::iterator it = DenticityInfoMap.begin();it!=DenticityInfoMap.end();it++)
		SitesRequiringDenticity.insert(it->second);//adding into this set all the site composite atoms that requires keeping track of denticity

	

}

void CHEMKinFiles::generateFiles()
{
	CreateCHEMKINfileAndWriteSpecies();

	if (!SetKinetics())
		cout<<"kinetics of reactions not set completely"<<endl;
	else
		WriteRxnInfoInCHEMKINfile();
}


int CHEMKinFiles::CreateCHEMKINfileAndWriteSpecies()
{
	ofstream file("CHEMKINsurf.txt");
	ofstream file2("CHEMKINspecies.txt");
	ofstream file3("CHEMKINgas.txt");

	FindAllElements();//find all the elements in the system
	file3<<"ELEMENTS"<<endl;
	file2<<"! This file contains CHEMKING specis names, its SMILES name, and thermo values: Enthalpy of formation (kcal/mol), Entropy (cal/mol-K), and cp (cal/mol-K) either at user specified temperature or at 298, 400, 500, 600, 800, and 1000 K."<<endl;

	//cout<<AllElements.size()<<endl;
	for (set<string>::iterator it = AllElements.begin();it!=AllElements.end();it++)
		file3<<"  "<<*it;
	file3<<endl;
	file3<<"END"<<endl;
	file3<<"SPECIES"<<endl;

	


	

	map<string,double>::iterator it;
	
	map<string*,int,classcomp>::iterator molIter;

	map<string, int> IdenticalCHEMKINstrings;

	set<string> CompositeHeterogeneousSites;

	for (int i =0;i<Network->CompositeSites.size();i++)
	{
		if (Network->CompositeSites[i].second == Heterogeneous)
			CompositeHeterogeneousSites.insert(Network->CompositeSites[i].first);
		
	}

	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{
		
		Molecule mol(*(molIter->first), moleculesize(*(molIter->first)));
		mol.unique_smiles();

		string CHEMKINstring ="";
		
		if (AllSites.count(mol.moleculestring())>0)
			CHEMKINstring = mol.GetMF();
		else
			CHEMKINstring = mol.getMFwithoutAtoms(CompositeHeterogeneousSites);
		
		if (IdenticalCHEMKINstrings.count(CHEMKINstring)>0)
		{
			IdenticalCHEMKINstrings[CHEMKINstring]+=1;
			CHEMKINstring+="-"+IntToStr(IdenticalCHEMKINstrings[CHEMKINstring]);

		}			
		else IdenticalCHEMKINstrings[CHEMKINstring]=1;
		
		CHEMKINSpeciesName[molIter->first] = CHEMKINstring;

		for (set<string>::iterator it = CompositeHeterogeneousSites.begin(); it!=CompositeHeterogeneousSites.end();it++)
		{
			if (("[" + *it + "]" == *(molIter->first))  || PopulateSiteOccupancyMap(molIter->first, *it))
				surfaceSpecies.insert(molIter->first);//if there exists site atoms of any kind, insert species into surfaceSpecies
			
		}

		if (surfaceSpecies.count(molIter->first)==0 && AllSites.count(*molIter->first)==0)
		{
			file3<<CHEMKINSpeciesName[molIter->first]<<endl;
		}
		if (surfaceSpecies.count(molIter->first)>0)
			PopulateHeterogeneousSiteOccupanyMap(molIter->first);
	}

	
	file3<<"END"<<endl;


	for (it =AllSites.begin();it!=AllSites.end();it++)
	{
	
		//string siteNameWithoutSqBrackets="";
		string* sitename =  StringRegistry::getStringPointer(it->first);


		
		if (HeterogeneousSitesOccupancyMap.count(it->first)>0)
		{
			file<<"SITE/SURFACE/      SDEN/"<<it->second<<"/"<<endl;
			multimap<string, pair<string*,int> >::iterator map_it;
			pair<multimap<string, pair<string*,int> >::iterator, multimap<string, pair<string*,int> >::iterator > pairIter;

			pairIter = HeterogeneousSitesOccupancyMap.equal_range(it->first);

			for (map_it=pairIter.first;map_it!=pairIter.second;map_it++)
			{

				CHEMKINSpeciesName[map_it->second.first]+="(S)";
				file<<CHEMKINSpeciesName[map_it->second.first]<<"/"<<map_it->second.second<<"/"<<endl;
			}


			CHEMKINBulkName[sitename] = CHEMKINSpeciesName[sitename]+"(B)";
			CHEMKINSpeciesName[sitename]+="(S)";
			file<<CHEMKINSpeciesName[sitename]<<endl;
			file<<"BULK "<<CHEMKINBulkName[sitename]<<"/     /"<<endl;
			file<<"END"<<endl;
		}		
	}
		



	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{

		file2<<CHEMKINSpeciesName[molIter->first]<<" --> "<<*molIter->first;

		Molecule mol(*(molIter->first), moleculesize(*(molIter->first)));
		mol.unique_smiles();

		//write thermo of each species into file 2 (enthalpy, entropy, and Cp at specified temps
		//TempsForCpCalc should be user specified. If not, calculate at 298, 400, 500, 600, 800, 1000K

		if (TempsForCpCalc.size()==0)
		{
			TempsForCpCalc.insert(298);
			TempsForCpCalc.insert(400);
			TempsForCpCalc.insert(500);
			TempsForCpCalc.insert(600);
			TempsForCpCalc.insert(800);
			TempsForCpCalc.insert(1000);
		}
		
		double dH=0.0;
		double dS = 0.0;
		bool isSiteMol = false;

		if (Network->containsSiteAtom(*molIter->first))isSiteMol = true;

		ThermoGA::calculateDeltaH(mol,isSiteMol,298.0,dH);
		ThermoGA::calculateDeltaS(mol,isSiteMol,298.0,dS);
		file2<<"  "<<dH/4.18;
		file2<<"  "<<dS/4.18;
		
		for (set<int>::iterator SetIt =TempsForCpCalc.begin();SetIt!=TempsForCpCalc.end();SetIt++)
		{
			double Cp = 0.0;

			ThermoGA::calculateCp(mol,isSiteMol,*SetIt,Cp);
			file2<<"  "<<Cp/4.18;
		}
		file2<<endl;		
	}




	
	file<<"REACTIONS	MWOFF	KCAL/MOLE    "<<endl;
	file<<endl;
	file<<"!-----------------------------------------------------------------------------	"<<endl;
	file<<"! Reactions begin now and have been grouped according to their families"<<endl;
	file<<"!-----------------------------------------------------------------------------	"<<endl;
	file<<endl;

	file3<<"REACTIONS"<<endl;

	file.close();
	file2.close();
	file3.close();

	return 0;
	
}

void CHEMKinFiles::PopulateHeterogeneousSiteOccupanyMap(string* mol)
{
	Molecule m(*mol, moleculesize(*mol));
	m.unique_smiles();
	for (map<string,double>::iterator it = AllSites.begin();it!=AllSites.end();it++)
	{
		pair<multimap<string,string>::iterator, multimap<string,string>::iterator> pairIter;

		pairIter = DenticityInfoMap.equal_range(it->first);

		int NumMatches = 0;
		for (multimap<string,string>::iterator m_it=pairIter.first;m_it!=pairIter.second;m_it++)
		{
			string pattern = m_it->second;
			NumMatches+= Patternmatch(m,Substructure(pattern,patternsize(pattern)),0).GetDistinctMatches();
		}
		if (NumMatches>0)
			HeterogeneousSitesOccupancyMap.insert(pair<string, pair<string*,int> >(it->first, pair<string*,int>(mol, NumMatches)));
			
	}
}



bool CHEMKinFiles::PopulateSiteOccupancyMap(string* mol, string site)
{
	if (!Network->CompositeSites.empty())
	{
		Molecule m(*mol, moleculesize(*mol));
		m.unique_smiles();
		int NumMatches = 0;
		 NumMatches = Patternmatch(m,Substructure(site,patternsize(site)),0).number_of_matches();

		if (NumMatches>0)
		{
			SiteOccupancyMap.insert(pair<string, pair<string*,int> >(site, pair<string*,int>(mol, NumMatches)));//TODO: alternatively consider taking the proper SMILES string instead of just "{M}"
			return true;
		}
	
		
	}
	return false;

}


int CHEMKinFiles::WriteRxnInfoInCHEMKINfile()
{
	ofstream file;
	ofstream file3;
	file.open("CHEMKINsurf.txt", ios::app);
	file3.open("CHEMKINgas.txt", ios::app);


	int maxRxnLength=0;
	for (int i =0;i<CHEMKINRxnStrings.size();i++)
	{
		int length = CHEMKINRxnStrings[i].length();
		if (maxRxnLength<=length)
			maxRxnLength = length;
	}
	int maxPreExpLength =0; 
	int maxIndexLength = 0; 
	int maxEaLength = 0;

	for (int i =0;i<CHEMKINPreExpInfo.size();i++)
	{
		int length = CHEMKINPreExpInfo[i].length();
		if (maxPreExpLength<=length)
			maxPreExpLength = length;
	}

	for (int i =0;i<CHEMKINIndexInfo.size();i++)
	{
		int length = CHEMKINIndexInfo[i].length();
		if (maxIndexLength<=length)
			maxIndexLength = length;
	}

	for (int i =0;i<CHEMKINActEInfo.size();i++)
	{
		int length = CHEMKINActEInfo[i].length();
		if (maxEaLength<=length)
			maxEaLength = length;
	}


	multimap<int, int>::iterator it;
	int prevRule = -1;

	int rxnCount =1;
	set<int> RxnsWrittenSoFar;

	for (it=Network->ReactionsMap.begin();it!=Network->ReactionsMap.end();it++)
	{
		
		if (it->first!=prevRule)
		{
			file<<endl;
			file<<"!----------------------------------------------------------"<<endl;
			file<<"            !Reactions of rule ";
			file<<Network->Rtlist[it->first].getRuleName()<<":"<<endl;
			file<<"!----------------------------------------------------------"<<endl;
			file<<endl;

			prevRule = it->first;
		}

		
		
		string str(maxRxnLength+5-CHEMKINRxnStrings[(*it).second-1].length(), ' ');
		string str2(maxPreExpLength+3-CHEMKINPreExpInfo[(*it).second-1].length(), ' ');
		string str3(maxIndexLength+3-CHEMKINIndexInfo[(*it).second-1].length(), ' ');
		string str4(maxEaLength+3-CHEMKINActEInfo[(*it).second-1].length(), ' ');

		//checking if reverse reaction already included or not! 
		bool ReverseAlreadyIncluded = false;
		for (set<int>::iterator setit=RxnsWrittenSoFar.begin();setit!=RxnsWrittenSoFar.end();setit++)
		{
			//cout<<"checking"<<endl;
			if (Network->AllReactions[it->second-1].isReverseReaction(Network->AllReactions[*setit]))
			{
				ReverseAlreadyIncluded = true;
				cout<<"yes it's reverse"<<endl;
				cout<<Network->AllReactions[it->second-1].reactionstring()<<endl;
				break;
			}
		}

		if (!ReverseAlreadyIncluded)
		{
			RxnsWrittenSoFar.insert(it->second-1);

			if (!isGasPhaseRxn(it->second-1))
			{
				WriteInfoFile(CHEMKINRxnStrings[(*it).second-1],file);
				WriteInfoFile(str,file);
				WriteInfoFile(CHEMKINPreExpInfo[(*it).second-1],file);
				WriteInfoFile(str2,file);
				WriteInfoFile(CHEMKINIndexInfo[(*it).second-1],file);
				WriteInfoFile(str3,file);
				WriteInfoFile(CHEMKINActEInfo[(*it).second-1],file);
				WriteInfoFile(str4,file);
				file<<"  !"<<rxnCount;
				file<<endl;
				if (RxnsWithStickingCoeff.count((*it).second-1)==1)
					file<<"STICK"<<endl;
			}
			else
			{
				WriteInfoFile(CHEMKINRxnStrings[(*it).second-1],file3);
				WriteInfoFile(str,file3);
				WriteInfoFile(CHEMKINPreExpInfo[(*it).second-1],file3);
				WriteInfoFile(str2,file3);
				WriteInfoFile(CHEMKINIndexInfo[(*it).second-1],file3);
				WriteInfoFile(str3,file3);
				WriteInfoFile(CHEMKINActEInfo[(*it).second-1],file3);
				WriteInfoFile(str4,file3);
				file3<<"  !"<<rxnCount;
				file<<endl;

			}
			rxnCount++;
		}

	}
	file<<"END"<<endl;
	file3<<"END"<<endl;

	return 0;


}


void CHEMKinFiles::FindAllElements()
{
	for (set<string*,classcomp>::iterator set_it=Network->InitialReactants.begin();set_it!=Network->InitialReactants.end();set_it++)
	{

		//cout<<"initial reactant: "<<*(*set_it)<<endl;
		set<string> elementlist = Molecule(*(*set_it), moleculesize(*(*set_it))).GetElements();
		//cout<<elementlist.size()<<endl;

		for (set<string>::iterator it = elementlist.begin();it!=elementlist.end();it++)
			AllElements.insert(*it);
	}
	//cout<<AllElements.size()<<endl;
}

bool CHEMKinFiles::isGasPhaseRxn(int rxn)
{
	//if any of hte reactants/products is a surface species or a site, then the reaction is a surface reaction. 
	for (int i =0;i<Network->AllReactions[rxn].number_pdcts();i++)
	{
		if (surfaceSpecies.count(Network->AllReactions[rxn].get_products(i))>0 || AllSites.count(*Network->AllReactions[rxn].get_products(i))>0)//it is a surface species or a site
			return false;
	}

	return true;
}

bool CHEMKinFiles::SetKinetics()
{

	ofstream BEPFile;
	BEPFile.open("BEPInfo.txt");
	map<pair<double,double>, int > BEPValues;
	multimap<int, int> BEPIndexAndRxnMap;


	for (int i=0;i<Network->AllReactions.size();i++)
	{
		
		double PreExp, ActE, TempIndex, kinValue;
		int rule = Network->AllReactions[i].get_rule();
		
		
		PreExp = 0.0; ActE=0.0; TempIndex = 0.0;kinValue = 0.0;
		bool calcK = false; bool usesBEP = false; bool usesLFER = false;
		double alpha = 0.0; double beta = 0.0; bool stick = false; bool calculateCollisionFreq = 0.0;

		double dH = Network->calculateThermoOfRxn(EnthalpyType,Network->AllReactions[i],Network->Temperature)*1000.0;
		double dS = Network->calculateThermoOfRxn(EntropyType,Network->AllReactions[i],Network->Temperature);
		double delN = Network->getDeltaNGasPhase(Network->AllReactions[i]);


		KineticsInfo kinInfo(*kinetics,Network->AllReactions[i],Rg,Network->Temperature, dH, dS, delN, siteDensity, Network->firstGasSpecies(i,0), Network->firstGasSpecies(i,1)); 

		kinInfo.setCollisionFreq(false);//PreExp will only be collission freq.

		if (!kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK, usesBEP, usesLFER, alpha, beta, stick))
			return false;

		
		CHEMKINRxnStrings.push_back(generateCHEMKINRxnStrings(i));

		if (usesBEP)
		{
			//cout<<"adding into BEP map"<<endl;
			pair<double,double> P(alpha,beta);
			if (BEPValues.count(P)==0)BEPValues[P]=BEPValues.size()+1;
			BEPIndexAndRxnMap.insert(pair<int,int> (BEPValues[P],i));
		}

		generateCHEMKINkineticsStrings(PreExp, TempIndex, ActE/4180.0, stick); 

				

	}
	for (map<pair<double,double>, int>::iterator it = BEPValues.begin();it!=BEPValues.end();it++)
		BEPFile<<it->second<<"  "<<it->first.first<<"  "<<it->first.second<<endl;

	for (multimap<int,int>::iterator it = BEPIndexAndRxnMap.begin();it!=BEPIndexAndRxnMap.end();it++)
		BEPFile<<CHEMKINRxnStrings[it->second]<<"  "<<it->first<<endl;


	BEPFile.close();
	return true;
}

string CHEMKinFiles::generateCHEMKINRxnStrings(int rxn)
{
	map<string,int> ReactantsStoich;
	map<string,int> ProductsStoich;
	set<string*, classcomp> SitesInRxn;


	for (int i =0;i<Network->AllReactions[rxn].number_reactants();i++)
	{
		string* molstrptr = Network->AllReactions[rxn].get_reactants(i);

		if (DenticityInfoMap.count(*molstrptr)==0)// if not a site
		{
		
			string CHEMKINname = CHEMKINSpeciesName[molstrptr];
			if (ReactantsStoich.count(CHEMKINname)==0)
				ReactantsStoich[CHEMKINname]=1;
			else ReactantsStoich[CHEMKINname]+=1;
		}
		else 
			SitesInRxn.insert(molstrptr);
		
			/*if (SitesInReactants.count(*molstrptr)==0)
			{
				SitesInReactants.insert(*molstrptr);
				string CHEMKINname = CHEMKINSpeciesName[molstrptr];

				//calculate site count;
				//add enough bulk species;
				ReactantsStoich[CHEMKINname] = -calculateNetSiteCount(*molstrptr,reactants, products);
				string bulkName = CHEMKINBulkName[molstrptr];
				ProductsStoich[bulkName] = -ReactantsStoich[CHEMKINname];

				
			}*/

	}

	for (int i =0;i<Network->AllReactions[rxn].number_pdcts();i++)
	{
		string* molstrptr = Network->AllReactions[rxn].get_products(i);
		
		if (DenticityInfoMap.count(*molstrptr)==0)// if not a site
		{
			string CHEMKINname = CHEMKINSpeciesName[molstrptr];
			if (ProductsStoich.count(CHEMKINname)==0)
				ProductsStoich[CHEMKINname]=1;
			else ProductsStoich[CHEMKINname]+=1;

			
		}
		else 
			SitesInRxn.insert(molstrptr);
		
			/*if (SitesInReactants.count(*molstrptr)==0)
			{
				SitesInProducts.insert(*molstrptr);
				
				string CHEMKINname = CHEMKINSpeciesName[molstrptr];

				//calculate site count;
				//add enough bulk species;
				ProductsStoich[CHEMKINname] = calculateNetSiteCount(*molstrptr,reactants, products);
				string bulkName = CHEMKINBulkName[molstrptr];
				ReactantsStoich[bulkName] = -ProductsStoich[CHEMKINname];
			}*/
		
			
	}

	for (set<string*, classcomp>::iterator it = SitesInRxn.begin();it!=SitesInRxn.end();it++)
	{
		int netCount = calculateNetSiteCount(*(*it), rxn);

		if (netCount<0)
		{
			string CHEMKINname = CHEMKINSpeciesName[*it];
			ProductsStoich[CHEMKINname] = -netCount;
			string bulkName = CHEMKINBulkName[*it];
			ReactantsStoich[bulkName] = ProductsStoich[CHEMKINname];
		}

		if (netCount >0)
		{
			string CHEMKINname = CHEMKINSpeciesName[*it];
			ReactantsStoich[CHEMKINname] = netCount;
			string bulkName = CHEMKINBulkName[*it];
			ProductsStoich[bulkName] = ReactantsStoich[CHEMKINname];
		}
	}

	string rxnstr="";

	map<string,int>::iterator it;

	for (it = ReactantsStoich.begin();it!=ReactantsStoich.end();it++)
	{
		if (it!=ReactantsStoich.begin())
			rxnstr+="+";
		if (it->second==1)
			rxnstr+=it->first;
		else rxnstr=rxnstr+IntToStr(abs(it->second))+it->first;
	}
	rxnstr+="=";
	for (it = ProductsStoich.begin();it!=ProductsStoich.end();it++)
	{
		if (it!=ProductsStoich.begin())
			rxnstr+="+";
		if (it->second==1)
			rxnstr+=it->first;
		else rxnstr=rxnstr+IntToStr(abs(it->second))+it->first;
	}

	return rxnstr;	
}

int CHEMKinFiles::calculateNetSiteCount(std::string site, int rxn)
{
	int siteCount=0;
	
	for (int i =0;i<Network->AllReactions[rxn].number_reactants();i++)
	{
		string* molstrptr = Network->AllReactions[rxn].get_reactants(i);

		
		if (DenticityInfoMap.count(site)>0 && site!=*molstrptr)//checking if site is indeed a surface site and if molstring is not site --because I want to count the number of sites "site" occupies (a little confusing, i know :D)
		{
			
			multimap<string, pair<string*,int> >::iterator it;
			pair <multimap<string, pair<string*,int> >::iterator, multimap<string, pair<string*,int> >::iterator> pairIter;
			pairIter = HeterogeneousSitesOccupancyMap.equal_range(site);

			for (it = pairIter.first;it!=pairIter.second;it++)
			{
				
				if (*(it->second.first)==*molstrptr)//if the species occupying the site matches reactants[i]
				{
					
					siteCount-=it->second.second;//count sites (subtract because reactants)
				}
			}
		}
		
	}

	for (int i =0;i<Network->AllReactions[rxn].number_pdcts();i++)
	{
		string* molstrptr = Network->AllReactions[rxn].get_products(i);

		if (DenticityInfoMap.count(site)>0 && site!=*molstrptr)//checking if site is indeed a surface site and if molstring is not site --because I want to count sites of 
		{
			multimap<string, pair<string*,int> >::iterator it;
			pair <multimap<string, pair<string*,int> >::iterator, multimap<string, pair<string*,int> >::iterator> pairIter;
			pairIter = HeterogeneousSitesOccupancyMap.equal_range(site);//(site.substr(1,site.length()-2));

			for (it = pairIter.first;it!=pairIter.second;it++)
			{
				if (*(it->second.first)==*molstrptr)//if the species occupying the site matches reactants[i]
					siteCount+=it->second.second;//count sites
			}
		}
	}
	
	return siteCount;
}

					
				




void CHEMKinFiles::generateCHEMKINkineticsStrings(double A, double n, double Ea, bool stick)
{
	
	if (stick)RxnsWithStickingCoeff.insert(CHEMKINPreExpInfo.size());
	CHEMKINPreExpInfo.push_back(DoubleToString(A));
	CHEMKINIndexInfo.push_back(DoubleToString(n));
	CHEMKINActEInfo.push_back(DoubleToString(Ea));
	
}
	

GAMSFiles::GAMSFiles(rxn_net_gen * net, std::vector<KineticParamPtr> & kinFns, std::set<string> & sites, map<string, double>& Density,  std::multimap<string,string> & denticity, char* file)
{
	Network = net;
	filename = file;
	kineticsFns = &kinFns;
	AllSites = sites;
	DenticityInfoMap = denticity;
	SiteDensity = Density;

	for (multimap<string,string>::iterator it = DenticityInfoMap.begin();it!=DenticityInfoMap.end();it++)
		SitesRequiringDenticity.insert(it->first);//adding into this set all the free heterogeneous sites that require keeping track of denticity
}

bool GAMSFiles::generateFiles()
{
	ofstream file(filename);

	map<string, pair<string, string> > TypesOfAtomicBonding;
	map<string, double> GammaValues;

	file<<"*------------------"<<endl;
	file<<"* Sets"<<endl;
	file<<"*------------------"<<endl;
	file<<"Sets"<<endl;

	file<<"i_g(i)	gaseous molecules /"<<endl;

	//find surface and gas phase molecules

	map<string*,int,classcomp>::iterator molIter;
	
	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{
		
		Molecule mol(*(molIter->first), moleculesize(*(molIter->first)));
		mol.unique_smiles();
		
		for (int i =0;i<Network->CompositeSites.size();i++)
		{
			if (Network->CompositeSites[i].second == Heterogeneous)
			{
				if (("[" + Network->CompositeSites[i].first + "]" == *(molIter->first))  || PopulateSiteOccupancyMap(molIter->first, Network->CompositeSites[i].first))
					surfaceSpecies.insert(molIter->first);//if there exists site atoms of any kind, insert species into surfaceSpecies
			}
		}

		if (surfaceSpecies.count(molIter->first)==0 && AllSites.count(*molIter->first)==0)
			file<<"\'"<< Network->SMILESGAMSSpeciesMap[molIter->first]<<"\'"<<endl;//write into file gas phase species!
	}
	file<<"/"<<endl;


	//CM4
	TypesOfAtomicBonding["CM4"] = pair<string, string> ("CM", "C[|-{M}|==4]");
	GammaValues["CM4"] = 1.0;
	//CM3
	TypesOfAtomicBonding["CM3"] = pair<string, string> ("CM", "C[|-{M}|==3]");
	GammaValues["CM3"] = 0.75;
	//CM2
	TypesOfAtomicBonding["CM2_1"] = pair<string, string> ("CM", "C[|-{M}|==2][!|=O|>=1]");
	GammaValues["CM2_1"] = 0.5;
	//CM2 for CO
	TypesOfAtomicBonding["CM2_2"] = pair<string, string> ("CM", "C[|-{M}|==2][|=O|>=1]");
	GammaValues["CM2_2"] = 0.5;
	//CM
	TypesOfAtomicBonding["CM1"] = pair<string, string> ("CM", "C[|-{M}|==1]");
	GammaValues["CM1"] = 0.25;
	//OM2
	TypesOfAtomicBonding["OM2"] = pair<string, string> ("OM", "O[|-{M}|==2]");
	GammaValues["OM2"] = 1.0;
	//OM1
	TypesOfAtomicBonding["OM1"] = pair<string, string> ("OM", "O[|-{M}|==1]");
	GammaValues["OM1"] = 0.5;
	//HM
	TypesOfAtomicBonding["HM"] = pair<string, string> ("HM", "{M}[|-H|==1]");
	GammaValues["HM"] = 1.0;
	//nOH
	TypesOfAtomicBonding["nOH"] = pair<string, string> ("OM", "O[|_{M}|==1]");
	GammaValues["nOH"] = 0.1;

	file<<"a	types of atomic site bonding /"<<endl;

	map<string, pair<string, string> >::iterator map_it;

	for (map_it = TypesOfAtomicBonding.begin();map_it!=TypesOfAtomicBonding.end();map_it++)
		file<<"\'"<<map_it->first<<"\'"<<endl;
	file<<"/"<<endl;

	file<<"aa	categories of atomic site bonding /"<<endl;
	file<<"\'CM\'"<<endl;
	file<<"\'OM\'"<<endl;
	file<<"\'HM\'"<<endl;
	file<<"/"<<endl;


	file<<"set a_corr(a,aa) corresponding type of atomic site bonding /"<<endl;

	for (map_it = TypesOfAtomicBonding.begin();map_it!=TypesOfAtomicBonding.end();map_it++)
		file<<"\'"<<map_it->first<<"\' . \'"<<map_it->second.first<<"\'"<<endl;

	file<<"/"<<endl;	

	file<<";"<<endl;


	file<<"*------------------"<<endl;
	file<<"* Parameters"<<endl;
	file<<"*------------------"<<endl;
	file<<"Parameters"<<endl;

	

	map<int,pair<double,double> > RxnBEPParamaters;
	map<int, pair<double,double> > RxnPreExpMap;

	for (int i=0;i<Network->AllReactions.size();i++)
	{
		
		double PreExp, ActE, TempIndex, kinValue;
		int rule = Network->AllReactions[i].get_rule();
		
		
		PreExp = 0.0; ActE=0.0; TempIndex = 0.0;kinValue = 0.0;
		bool calcK = false; bool usesBEP = false; bool usesLFER = false;
		double alpha = 0.0; double beta = 0.0; bool stick = false;
		double dH = Network->calculateThermoOfRxn(EnthalpyType,Network->AllReactions[i],Network->Temperature)*1000.0;
		double dS = Network->calculateThermoOfRxn(EntropyType,Network->AllReactions[i],Network->Temperature);
		double delN = Network->getDeltaNGasPhase(Network->AllReactions[i]);

		//TODO: siteDen is density of sites - right now assuming only one type of site exists. Need to change that soon! 
		double siteDen = 0.1;

		if (SiteDensity.size()>0)
			siteDen = SiteDensity.begin()->second;



		KineticsInfo kinInfo(*kineticsFns,Network->AllReactions[i],1/Network->Temperature,Network->Temperature, dH, dS, delN, siteDen, Network->firstGasSpecies(i,0), Network->firstGasSpecies(i,1)); //Rg is 1/Temp because everything is in pressure units and I want RT to be 1.0!!
		if (!kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK, usesBEP, usesLFER, alpha, beta, stick))
			return false;

		int Symmetry = Network->AllReactions[i].get_frequency();

		if (!usesBEP) beta = ActE;
		

		if (calcK)
		{
			if (usesBEP) alpha = 1.0 - alpha;
			Symmetry = 1; //exp(ds/R) is multiplied to preExp factor. That accounts for the correct symmetry number to be used
		}
		RxnPreExpMap.insert(pair<int, pair<double,double> > (i, pair<double,double> (PreExp, Symmetry)));

		RxnBEPParamaters.insert(pair<int, pair<double,double> > (i,pair<double,double> (alpha,beta)));
			
	}

	file<<"alpha(j)	coefficient for enthalpy in calculating Ea /"<<endl;
	for (int i=0;i<Network->AllReactions.size();i++)
		file<<"\'rxn"<<i+1<<"\' "<<RxnBEPParamaters[i].first<<endl;
	file<<"/"<<endl;

	file<<"beta(j)		constant in calculating Ea /"<<endl;

	for (int i=0;i<Network->AllReactions.size();i++)
		file<<"\'rxn"<<i+1<<"\' "<<RxnBEPParamaters[i].second<<endl;
	file<<"/"<<endl;

	map<string, double>::iterator gamma_it;

	file<<"gamma(a) /"<<endl;
	for (gamma_it = GammaValues.begin();gamma_it!=GammaValues.end();gamma_it++)
		file<<"\'"<<gamma_it->first<<"\' "<<gamma_it->second<<endl;
	file<<"/"<<endl;

	file<<"prefactor(j)	pre-exponential factor A for reaction j /"<<endl;

	for (int i=0;i<Network->AllReactions.size();i++)
		file<<"\'rxn"<<i+1<<"\' "<<RxnPreExpMap[i].first<<endl;
	file<<"/"<<endl;

	file<<"symfactor(j)	factor for reaction j that accounts for symmetry/"<<endl;


	for (int i=0;i<Network->AllReactions.size();i++)
		file<<"\'rxn"<<i+1<<"\' "<<RxnPreExpMap[i].second<<endl;
	file<<"/"<<endl;

	file<<"freq(i,a)	frequency of atomic site bonding a for molecule i /"<<endl;

	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{
		if (surfaceSpecies.count(molIter->first)>=0 || AllSites.count(*molIter->first)>=0)
		{
			Molecule mol(*(molIter->first), moleculesize(*(molIter->first)));
			mol.unique_smiles();
			for (map_it = TypesOfAtomicBonding.begin();map_it!=TypesOfAtomicBonding.end();map_it++)
			{
				Substructure Sub(map_it->second.second,patternsize(map_it->second.second));

				Patternmatch P(mol,Sub,0);
				int matches = P.GetDistinctMatches();

				if (matches>0)
					file<<"\'"<< Network->SMILESGAMSSpeciesMap[molIter->first]<<"\' . \'"<<map_it->first<<"\' "<<matches<<endl;
			}
		}
	}
	file<<"/"<<endl;

	file<<"siteOccupancy(i)    sites occupied by each species i /"<<endl;

	//right now there is only one type of site. In general, there could be many, and maybe then siteOccupancy will be two dimensional (i,b) where b is set of all site types
	//for each species calculate the total number of sites it occupies. for others (gas phase) sites are zero! 

	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{
		if (surfaceSpecies.count(molIter->first)>0)
		{
			int NumMatches = 0;
			pair<multimap<string,string>::iterator, multimap<string,string>::iterator > pairIter;

			Molecule m(*molIter->first, moleculesize(*molIter->first));
			m.unique_smiles();

			for (set<string>::iterator set_it = SitesRequiringDenticity.begin();set_it!=SitesRequiringDenticity.end();set_it++)
			{
				

				pairIter = DenticityInfoMap.equal_range(*set_it);

				for (multimap<string,string>::iterator it = pairIter.first;it!=pairIter.second;it++)
				{
					string pattern = it->second;
					NumMatches += Patternmatch(m,Substructure(pattern,patternsize(pattern)),0).GetDistinctMatches();
				}
			}
			if (NumMatches>0)
				file<<"\'"<< Network->SMILESGAMSSpeciesMap[molIter->first]<<"\' "<<NumMatches<<endl;
			else
				file<<"\'"<< Network->SMILESGAMSSpeciesMap[molIter->first]<<"\' 0"<<endl;

			
		}
		else 
			file<<"\'"<< Network->SMILESGAMSSpeciesMap[molIter->first]<<"\' 0"<<endl;
	}




	file<<"/"<<endl;

	file<<"Q_Pt(aa)	binding energy contribution of atomic site bonding aa on platinum /"<<endl;
	file<<"\'CM\' 689.7"<<endl;
	file<<"\'OM\' 375.37"<<endl;
	file<<"\'HM\' 234.5"<<endl;
	file<<"/"<<endl;

	return true;
}



bool GAMSFiles::PopulateSiteOccupancyMap(string* mol, string site)
{
	if (!Network->CompositeSites.empty())
	{
		Molecule m(*mol, moleculesize(*mol));
		int NumMatches = 0;
		NumMatches = Patternmatch(m,Substructure(site,patternsize(site)),0).number_of_matches();
		if (NumMatches>0)
		{
			SiteOccupancyMap.insert(pair<string, pair<string*,int> >(site, pair<string*,int>(mol, NumMatches)));
			return true;
		}
	
		
	}
	return false;

}			


ThermoGroupsIdentification::ThermoGroupsIdentification(vector<pair<string,string> >& GroupFragments, vector<pair<string,string> >& GroupCorrections, char* speciesInput, char* outputfile)
{
	ifstream speciesfile(speciesInput);
	ofstream resultsFile(outputfile);

	if (speciesfile.is_open())
    {
		string linestr;
		int linecounter = 1;
		while ( speciesfile.good() )
		{
		   getline(speciesfile, linestr);
		   cout<<"Checking species in line "<<linecounter<<endl;
		   linestr = linestr.substr(0,  linestr.find_first_of("\t\r\n\v\f ",0));
		   linecounter++;

		    Molecule mol(linestr, moleculesize(linestr));
		    mol.unique_smiles();

			resultsFile<<mol.moleculestring();
			

			

			for (int i=0;i<GroupFragments.size();i++)
			{
				string fragment = GroupFragments[i].second;

				Substructure Sub(fragment, patternsize(fragment));
				Patternmatch P(mol,Sub,0);
				set<int> atomsCoveredByFirstFragAtom = P.AtomsCoveredBySubstr(0); // the set of atoms matching the first atom of the pattern described by the substructure S
			

				resultsFile<<"   "<<GroupFragments[i].first<<" => "<<atomsCoveredByFirstFragAtom.size();
			}

			for (int i=0;i<GroupCorrections.size();i++)
			{
				string fragment = GroupCorrections[i].second;

				Substructure Sub(fragment, patternsize(fragment));
				Patternmatch P(mol,Sub,0);
				
				resultsFile<<"   "<<GroupCorrections[i].first<<" => "<<P.GetDistinctMatches();
			}
			resultsFile<<endl;
		}
	}
	else
	{
		cout<<"Species file not found"<<endl;
	}

}

		
