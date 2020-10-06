#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <cctype>
#include <deque>
#include <bitset>
#include <iterator>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Classheader.h"
#include "AdditionalFunctions.h"
#include "StringRegistry.h"


Pathways::Pathways(rxn_net_gen* net, PathwayConstraints& P, const char* name)
{
	Network = net;
	Constr = P;
	DistinctNature = false;
	filename = name;
	canCalculateActE = false;
}

void Pathways::setDistinctNature(bool a)
{
	DistinctNature = a;
}


void Pathways::GeneratePathways()
{
	
	vector<string*> FinalProducts;

	multimap<string*,int>::iterator it;

	for (it=Network->MolProductMap.begin();it!=Network->MolProductMap.end();it++)
	{
		if (Network->MolReactantMap.count((*it).first)==0)
		{
			FinalProducts.push_back((*it).first);
		}
	}
	

	if (FinalProducts.size()==0)cout<<"no molecule exists that seems to be a final product"<<endl;
	else
	{
		for (int i=0;i<FinalProducts.size();i++)
		{
			multimap<string*, vector<int> > Pathways;
			int SpeciesRank = Network->AllMolecules[it->first];
			int RelLength = Constr.getMaxRelativeLength();
			int ActualMaxLength = 0;
			if (RelLength ==-1)ActualMaxLength = Constr.getMaxLength();
			else ActualMaxLength = SpeciesRank + RelLength;
			Pathways = FindAllPathways(FinalProducts[i], ActualMaxLength);	
			ofstream myfile;
			multimap<int,string*>ValuesMap;
			if (i==0)
				myfile.open(filename);
			else 
				myfile.open(filename,ios::app);
			multimap<string*,vector<int> >::iterator map_it;
			for (map_it=Pathways.begin();map_it!=Pathways.end();map_it++)
			{
				vector<generated_rxn> Rxn;
				vector<int>::reverse_iterator it2;
				for (it2=(*map_it).second.rbegin();it2<(*map_it).second.rend();++it2)
				{
					Rxn.push_back(Network->AllReactions[(*it2)-1]);
				}
				RxnPathway RP(&Rxn);
				
				
				if (Constr.IsSatisfied(RP))
				{
					pair<int, map<int,int> > PValue;
					PValue = CalculatePathwayValue((*map_it).second);
					if ((!similarRulePresent(&ValuesMap, PValue.first, map_it->first) && DistinctNature) || (!DistinctNature))
					{
					
						for (it2=(*map_it).second.rbegin();it2<(*map_it).second.rend();++it2)
						{
							myfile<<Network->AllReactions[(*it2)-1].reactionstring()<<endl;
						}
						myfile<<endl;
						ValuesMap.insert(pair<int,string*>(PValue.first,map_it->first));
					}
				}	
			}
			myfile.close();
		}
	}

	
}
void Pathways::QueryPlusGenerate(ConstrPtr Cons)
{
	multimap<string*, int>::iterator it;
	int Counter = 0;
	string previous ="";
	ofstream pathwayfile;

	ofstream uniquePathwayRxns;
	uniquePathwayRxns.open("uniquePathwayRxnSet.txt");
				
	
	string molecule_strfile(filename);
	cout<<"pathways to be stored in "<<molecule_strfile<<endl;
	int molcount = 0;
	pathwayfile.open(filename);
	for (it=Network->MolProductMap.begin();it!=Network->MolProductMap.end();it++)
	{
		
		vector<Molecule> M;
		
		bool flag = true;
		
		if (previous==*((*it).first))flag = false;
		else
		{
			previous = *((*it).first);
			Molecule mol(*((*it).first), moleculesize(*((*it).first)));
			
			
			flag = (Cons)(mol);
			
		}
		if (flag)
		{
			
			molcount++;
			Counter++;
			multimap<string*, vector<int> > Pathways;
			
			pathwayfile<<"molecule "<<*(it->first)<<" satisfies constraints"<<endl;
			cout<<"molecule "<<*(it->first)<<" satisfies constraints"<<endl;
			pathwayfile<<"shortest path to the molecule is "<<Network->AllMolecules[(*it).first]<<endl;
			cout<<"shortest path to the molecule is "<<Network->AllMolecules[(*it).first]<<endl;

			int SpeciesRank = Network->AllMolecules[it->first];
			int RelLength = Constr.getMaxRelativeLength();
			int ActualMaxLength = 0;
			if (RelLength ==-1)ActualMaxLength = Constr.getMaxLength();
			else ActualMaxLength = SpeciesRank + RelLength;
			
			if (SpeciesRank <=ActualMaxLength)
			{
				cout<<"finding pathways"<<endl;
				
				Pathways = FindAllPathways((*it).first, ActualMaxLength);
				
				map<pair<int,string*>,int > ValuesMap;
				map<int,map<int,int> > ValueInfoMap;
				
				pathwayfile<<"Pathways for "<<*(*it).first<<" "<<endl;
				uniquePathwayRxns<<"Unique reaction set for "<<*(it->first)<<" is as below"<<endl;
				
				pathwayfile<<endl;
				multimap<string*,vector<int> >::iterator map_it;
				int valid_path_count=0;
				int distinct_path_count = 0;
				map<int,int> LengthCounterMap;//counts how many pathways of a given length are found
				
				set<int> OnlyUniqueRxnsSet;//this stores the unique rxns - most pathways have rxns that have occured already - so the set of unique rxns will be much smaller than the total number of reactions in all pathways
				for (map_it=Pathways.begin();map_it!=Pathways.end();map_it++)
				{
					vector<generated_rxn> Rxn;
					vector<int> RxnIndices;
					vector<int>::reverse_iterator it2;
					for (it2=(*map_it).second.rbegin();it2<(*map_it).second.rend();++it2)
					{
						Rxn.push_back(Network->AllReactions[(*it2)-1]);
						RxnIndices.push_back((*it2)-1);
						
					}

					vector<double> ActEnergy;
					

									
					RxnPathway RP(&Rxn);
					if (canCalculateActE)
					{
						calculateActEnergy(RxnIndices, ActEnergy);
						RP.setRxnActE(&ActEnergy);
					}
					if (Constr.IsSatisfied(RP))
					{
					
						valid_path_count++;
						pair<int,map<int,int> > PValue;
			 			PValue = CalculatePathwayValue((*map_it).second);
						pair<int,string*> ValueMolPair(PValue.first,map_it->first);
						if ((ValuesMap.count(ValueMolPair)==0 && DistinctNature) || (!DistinctNature))
						{
							distinct_path_count++;
							int act_index =0;
							for (it2=(*map_it).second.rbegin();it2<(*map_it).second.rend();++it2)
							{
								pathwayfile<<Network->AllReactions[(*it2)-1].reactionstring();
								if(canCalculateActE)pathwayfile<<"  "<<ActEnergy[act_index]<<endl;
								if (OnlyUniqueRxnsSet.find((*it2))==OnlyUniqueRxnsSet.end())
								{
									OnlyUniqueRxnsSet.insert((*it2));
									uniquePathwayRxns<<Network->AllReactions[(*it2)-1].reactionstring()<<endl;
									
								}
								act_index++;
								pathwayfile<<endl;
							}
							pathwayfile<<endl;
							//writing the pathway information at the end
							pathwayfile<<"the above pathway is of length "<<Rxn.size()<<" and contains"<<endl;
							if (LengthCounterMap.count(Rxn.size())==0)
								LengthCounterMap[Rxn.size()]=1;
							else LengthCounterMap[Rxn.size()]+=1;


							map<int,int>::iterator RuleMap_it;
							for (RuleMap_it=PValue.second.begin();RuleMap_it!=PValue.second.end();RuleMap_it++)
							{
								pathwayfile<<"rule "<<Network->Rtlist[RuleMap_it->first].getRuleName()<<" occurs "<<RuleMap_it->second<<" times"<<endl;
							}
							pathwayfile<<endl;
							
							if (ValuesMap.count(ValueMolPair)==0)ValuesMap[ValueMolPair]=0;
							ValueInfoMap[PValue.first]=PValue.second;
						}
						
						if (ValuesMap.count(ValueMolPair)>0)
							ValuesMap[ValueMolPair]+=1;
					
				
					}
				}
				pathwayfile<<"total valid pathways is "<<valid_path_count<<endl;
				pathwayfile<<"total distinct pathways is "<<distinct_path_count<<endl;
				
				map<pair<int,string*>,int>::iterator ValuesMap_it;
				for (ValuesMap_it=ValuesMap.begin();ValuesMap_it!=ValuesMap.end();ValuesMap_it++)
				{
					pathwayfile<<ValuesMap_it->second<<" pathways from "<<*(ValuesMap_it->first.second)<<" containing "<<endl;
					map<int,int>::iterator RuleMap_it;
					for (RuleMap_it=ValueInfoMap[ValuesMap_it->first.first].begin();RuleMap_it!=ValueInfoMap[ValuesMap_it->first.first].end();RuleMap_it++)
					{
						pathwayfile<<RuleMap_it->second<< "of rule "<<Network->Rtlist[RuleMap_it->first].getRuleName()<<endl;
					}
					pathwayfile<<endl;
				}
				
				pathwayfile<<endl;

				cout<<"valid pathways "<<valid_path_count<<", distinct pathways "<<distinct_path_count<<endl;
				map<int,int>::iterator count_map;

				for (count_map=LengthCounterMap.begin();count_map!=LengthCounterMap.end();count_map++)
				{
					pathwayfile<<"Number of pathways of length "<<count_map->first<<" is: "<<count_map->second<<endl;
				}


				uniquePathwayRxns<<endl;
				
				uniquePathwayRxns<<"Number of reactions is "<<OnlyUniqueRxnsSet.size()<<endl;
				uniquePathwayRxns<<endl;
				uniquePathwayRxns<<endl;

			}
			else
			{
				pathwayfile<<"could not find pathways to the molecule because shortest path length is "<<Network->AllMolecules[(*it).first]<<" and cannot be reached in "<<Constr.getMaxLength()<<" steps"<<endl;
				pathwayfile<<"pathway identification for this molecule has been skipped"<<endl;
				pathwayfile<<endl;

				cout<<"could not find pathways to the molecule because shortest path length is "<<Network->AllMolecules[(*it).first]<<" and cannot be reached in "<<Constr.getMaxLength()<<" steps"<<endl;
				cout<<"pathway identification for this molecule has been skipped"<<endl;
			}
			
		}
		
		
	}
	pathwayfile<<molcount<<" molecules satisfied the molecular constraints specified"<<endl;
	cout<<molcount<<" molecules satisfied the molecular constraints specified"<<endl;
	cout<<""<<endl;
	pathwayfile.close();
	
}

void Pathways::AddKineticFunctions(std::vector<KineticParamPtr> & KineticFunc)
{
	KineticFunctions = &KineticFunc;
	canCalculateActE = true;
}

bool Pathways::calculateActEnergy(vector<int>& rxns, vector<double>& actE)
{
	actE.clear();

	for (int i=0;i<rxns.size();i++)
	{		
		double PreExp, ActE, TempIndex, kinValue;
		int rule = Network->AllReactions[rxns[i]].get_rule();
		
		//PreExp = PreExpFactor->at(rule);
		//ActE = ActEnergy->at(rule);
		//TempIndex = N->at(rule);
		PreExp = 0.0; ActE=0.0; TempIndex = 0.0; kinValue = 0.0;
		bool calcK = false; bool usesBEP = false; bool usesLFER = false; double alpha = 0.0; double beta = 0.0; bool stick = false;

		double dH = Network->calculateThermoOfRxn(EnthalpyType,Network->AllReactions[rxns[i]],Network->Temperature)*1000.0;
		double dS = Network->calculateThermoOfRxn(EntropyType,Network->AllReactions[rxns[i]],Network->Temperature);
		double delN = Network->getDeltaNGasPhase(Network->AllReactions[rxns[i]]);

		KineticsInfo kinInfo(*KineticFunctions,Network->AllReactions[rxns[i]],0.0821,Network->Temperature, dH, dS, delN, 0.0, Network->firstGasSpecies(i,0), Network->firstGasSpecies(i,1)); 


		if (!kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK,usesBEP,usesLFER, alpha, beta, stick))
			return false;
		actE.push_back(ActE);
	}
	return true;

}


multimap<string*, vector<int> > Pathways::FindAllPathways(string* S, int ActualMaxLength)
{
	vector<int> Reactions;
	
	multimap<string*, vector<int> > GenPathways;
	vector<string*> PrimaryMolecules;
	vector<int> IndivRxnCounter;
	map<string*,int> RxnCounterMap;//stores how many formation reactions of the molecule has been used in the pathway thus far
	PrimaryMolecules.push_back(S);
	RxnCounterMap.insert(pair<string*,int>(S,0));
	IndivRxnCounter.push_back(0);
	
	while (PrimaryMolecules.size()>0 )
	{
		
		string* P = PrimaryMolecules[PrimaryMolecules.size()-1];
		
		map<string*,int>::iterator it2;
		it2=RxnCounterMap.find(P);//it will definitely find one
		int counter = (*it2).second;

		if (Network->InitialReactants.count(P)>0)
		{
			
			int totalpathwaycost = CalculateCost(Reactions);
			
			if (totalpathwaycost<=Constr.getMaxCost() && totalpathwaycost>=Constr.getMinCost() && Reactions.size()>=Constr.getMinLength())
				GenPathways.insert(pair<string*,vector<int> >(P,Reactions));
			

			PrimaryMolecules.pop_back();
			IndivRxnCounter.pop_back();
			if (Reactions.size()>0)Reactions.pop_back();
			
		}
		else if (Network->MolProductMap.count(P)>counter && Reactions.size()<=ActualMaxLength-1)
		{
			
			pair<multimap<string*,int>::iterator,multimap<string*,int>::iterator> ret;
			ret =Network->MolProductMap.equal_range(P);
			advance (ret.first,counter);
			string* NewMolecule = CloserMolecule((*ret.first).second,P);
			
			(*it2).second++;
			IndivRxnCounter[IndivRxnCounter.size()-1]++;
			if (Network->AllMolecules[NewMolecule]<=(ActualMaxLength - Reactions.size()))
			{
				PrimaryMolecules.push_back(NewMolecule);
				
				if (!HasCycles(PrimaryMolecules) && !isReverseAlreadyPresent(Reactions,(*ret.first).second))
				{
					Reactions.push_back((*ret.first).second);
					
					IndivRxnCounter.push_back(0);
					
					if (RxnCounterMap.count(NewMolecule)==0)
						RxnCounterMap.insert(pair<string*,int>(NewMolecule,0));
					
					it2=RxnCounterMap.find(NewMolecule);
				}
				else PrimaryMolecules.pop_back();
			}
			
		}
		else
		{
			
			PrimaryMolecules.pop_back();
			(*it2).second-=IndivRxnCounter.back();
			IndivRxnCounter.pop_back();
			if (Reactions.size()>0)Reactions.pop_back();
		}
	}
	return GenPathways;
}

string* Pathways::CloserMolecule(int i, std::string * S)
{
	return Network->AllReactions[i-1].get_reactants(Network->AllReactions[i-1].getParentMolecule(S));
}

void Pathways::SetConstraints(PathwayConstraints PC)
{
	Constr = PC;
}


int Pathways::CalculateCost(vector<int> R)
{
	int TotalCost = 0;

	for (int i=0;i<R.size();i++)
	{
		TotalCost+=Network->Rtlist[Network->AllReactions[R[i]-1].get_rule()].getCost();
	}
	
	return (TotalCost);
}

bool Pathways::isReverseAlreadyPresent(vector<int> R, int rxn)
{
	bool isPresent = false;

	for (int i=0;i<R.size();i++)
	{
		if (Network->AllReactions[R[i]-1].isReverseReaction(Network->AllReactions[rxn-1]))
		{
			isPresent = true;
			break;
		}
	}
	return isPresent;
}

pair<int,map<int,int> > Pathways::CalculatePathwayValue(vector<int> R)
{
	int Value = 1;
	
	map<int,int> RuleFreqMap;

	for (int i=0;i<R.size();i++)
	{
		int rule = Network->AllReactions[R[i]-1].get_rule();
		Value = Value*prime(rule);
		if (RuleFreqMap.count(rule)>0)
			RuleFreqMap[rule]+=1;
		else RuleFreqMap[rule]=1;
	}

	return pair<int,map<int,int> >(Value, RuleFreqMap);
}

bool Pathways::HasCycles(std::vector<string*> v)
{
	bool result = false;
	for (int i=0;i<v.size();i++)
	{
		for (int j=0;j<i;j++)
		{
			if ((*v[i])==(*v[j]))
			{
				result = true;
				break;
			}
		}
	}
	return result;

}

bool Pathways::similarRulePresent(std::multimap<int,string*>* ValueMap, int PValue, std::string * firstMolecule)
{
	if (ValueMap->count(PValue)==0)
		return false;
	else 
	{
		pair<multimap<int, string*>::iterator,multimap<int,string*>::iterator> ret;

		ret = ValueMap->equal_range(PValue);

		multimap<int,string*>::iterator it;
		
		bool returnvalue = false;
		for (it=ret.first;it!=ret.second;it++)
		{
			if ((*(it->second))==*firstMolecule)
			{
				returnvalue = true;
				break;
			}
		}
		return returnvalue;
	}
}

PathwayConstraints::PathwayConstraints()
{
	MaxLength = 20;//long enough, but not too large! 
	MinLength = 0;
	MaxCost = 10000000;//arbitrarily large! 
	MinCost = 0;
	MaxRelativeLength = -1;
	
}

void PathwayConstraints::SetMaxLength(int m)
{
	MaxLength = m;
}

int PathwayConstraints::getMaxLength()
{
	return MaxLength;
}

int PathwayConstraints::getMinCost()
{
	return MinCost;
}
int PathwayConstraints::getMinLength()
{
	return MinLength;
}

void PathwayConstraints::SetMaxRelativeLength(int m)
{
	MaxRelativeLength = m;
}

int PathwayConstraints::getMaxRelativeLength()
{
	return MaxRelativeLength;
}

void PathwayConstraints::SetMaxCost(int m)
{
	MaxCost = m;
}
void PathwayConstraints::SetMinLengthAndCost(int a, int b)
{
	MinLength = a;
	MinCost = b;
}
int PathwayConstraints::getMaxCost()
{
	return MaxCost;
}

bool PathwayConstraints::IsSatisfied(RxnPathway & RP)
{
	return (PathConstraints)(RP);
}


RxnInfo::RxnInfo(generated_rxn * rxn, double Enth, double Entr, double FrEn, double ActE)
{
	rxn = rxn;
	Enthalpy = Enth;
	Entropy = Entr;
	FreeEnergy = FrEn;
	ActEnergy = ActE;
}

double RxnInfo::ActE()
{
	return ActEnergy;
}

RxnPathway::RxnPathway(std::vector<generated_rxn> * Rxn)
{
	rxn_index=Rxn;
}

void RxnPathway::setRxnActE(vector<double> * ActE)
{
	RxnActE = ActE;
}


int RxnPathway::numRule(int Rule)
{
	int Counter=0;
	for (int i = 0;i<(*rxn_index).size();i++)
	{
		if ((*rxn_index)[i].get_rule()==Rule)Counter++;
	}
	return Counter;
}
int RxnPathway::numRuleWithParticipants(int Rule, ConstrPtr c, int p)
{
	int counter = 0;
	for (int i=0;i<(*rxn_index).size();i++)
	{
		if ((*rxn_index)[i].get_rule()==Rule)
		{
			
			int number;
			if (p==0)number = (*rxn_index)[i].number_reactants();
			else number = (*rxn_index)[i].number_pdcts();

			for (int j=0;j<number;j++)
			{
				string mol="";
				if (p==0)mol = *(*rxn_index)[i].get_reactants(j);
				else mol = *(*rxn_index)[i].get_products(j);
				
				if ((c)(Molecule(mol,moleculesize(mol))))
				{
					counter++;
					
					break;
				}
				
			}
		}
	}
	return counter;
}


int RxnPathway::numRuleWithParticipants(int Rule, CombinedConstrPtr c, int p)
{
	int counter = 0;
	
	for (int i=0;i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].get_rule()==Rule)
		{
			if (p==0 && (*rxn_index)[i].number_reactants()==2)
			{
				
				Molecule mol1(*(*rxn_index)[i].get_reactants(0), moleculesize(*(*rxn_index)[i].get_reactants(0)));
				Molecule mol2(*(*rxn_index)[i].get_reactants(1), moleculesize(*(*rxn_index)[i].get_reactants(1)));
				if ((c)(mol1,mol2))
				{
					counter++;
					break;
				}
			}
			else
			{
				bool satisfied=false;
				for (int j=0;j<(*rxn_index)[i].number_pdcts()-1;j++)
				{
					for (int k=j+1;k<(*rxn_index)[i].number_pdcts();j++)
					{

						Molecule mol1(*(*rxn_index)[i].get_products(j), moleculesize(*(*rxn_index)[i].get_products(j)));
						Molecule mol2(*(*rxn_index)[i].get_products(k), moleculesize(*(*rxn_index)[i].get_products(k)));
						if ((c)(mol1,mol2))
						{
							counter++;
							satisfied = true;
							break;
						}
					}
					if (satisfied)break;
				}
				
			}
		}
	}
	return counter;
}

int RxnPathway::numMolecule(ConstrPtr C, int p)
{
	int counter =0;
	if (p==0)
	{
		for (int i=0;i<(*rxn_index).size();i++)
		{
			int number = (*rxn_index)[i].number_reactants();
			for (int j=0;j< number ; j++)
			{
				string mol = *(*rxn_index)[i].get_reactants(j);
				if ((C)(Molecule(mol,moleculesize(mol))))
				{
					counter++;
					break;
				}
			}
		}
	}
	else 
	{
		for (int i=0;i<(*rxn_index).size();i++)
		{
			int number = (*rxn_index)[i].number_pdcts();
			for (int j=0;j< number ; j++)
			{
				string mol = *(*rxn_index)[i].get_products(j);
				if ((C)(Molecule(mol,moleculesize(mol))))
				{
					counter++;
					break;
				}
			}
		}
	}
	return counter;
}

int RxnPathway::numRuleIntramolecular(int r)
{
	int counter = 0;
	for (int i=0;i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].get_rule()==r)
		{
			if ((*rxn_index)[i].number_reactants()==1)//assuming that the bimolecularity  check has already been done! 
				counter++;
		}
	}
	return counter;
}

int RxnPathway::numIntramolecular()
{
	int counter = 0;
	for (int i=0;i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].isReactionIntramolecular())
			counter++;
	}
	return counter;
}

int RxnPathway::numSatisfying(int r, RxnConstraintPtr rcp)
{
	int counter = 0;
	for (int i =0; i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].get_rule()==r)
		{
			RxnInfo Rinfo(&rxn_index->at(i), 0.0, 0.0, 0.0,RxnActE->at(i));
			if ((rcp)(Rinfo))
				counter++;
		}
	}
	return counter;
}
		


void PathwayConstraints::AddPathwayConstrPtr(PathwayConstrPtr P)
{
	PathConstraints = P;
}


//...............................................
//MoleculeRxnQuery begins

MoleculeQuery::MoleculeQuery(rxn_net_gen* Net, ConstrPtr Cp, const char* file, bool generateGAMSfile, bool calculateProperties)
{
	Network = Net;
	filename = file;

	map<string*, int, classcomp>::iterator it;
	ofstream queryfile;
	string molecule_strfile(filename);
	cout<<"molecule queries to be stored in "<<molecule_strfile<<endl;
	queryfile.open(filename);
	ofstream CMCFile, CPFile, HLBFile, SurfTensFile;
	CMCFile.open("CMCValues.txt");
	CPFile.open("CPValues.txt");
	HLBFile.open("HLBValues.txt");
	SurfTensFile.open("SurfTens.txt");
	int count=0;
	queryfile<<"the following molecules satisfy all molecular constraints imposed"<<endl;
	for (it=Network->AllMolecules.begin();it!=Network->AllMolecules.end();it++)
	{
		Molecule mol(*((*it).first), moleculesize(*((*it).first)));
		if((Cp)(mol))
		{
			queryfile<<mol.getsize()<<"  "<<(*(it->first))<<endl;
			PassedMolecules.push_back(it->first);

			if (calculateProperties)
			{

				string GAMSName = Network->SMILESGAMSSpeciesMap[it->first];

				CMCFile<<GAMSName<<"  ";
				CPFile<<GAMSName<<"  ";
				HLBFile<<GAMSName<<"  ";
				SurfTensFile<<GAMSName<<"  ";

			
				for (int EO=4; EO<=12;EO++)
				{
					
					Surfqsar Sf(&mol, EO);
					CMCFile<<pow(10.0, Sf.calculateCMC())*(float(mol.MolecularWeight())+44.0*float(EO))*1000;
					CPFile<<Sf.calculateCP();
					HLBFile<<Sf.calculateHLB();
					SurfTensFile<<Sf.calculateSurfTens();
					if (EO!=12)
					{
						CMCFile<<"  ";
						CPFile<<"  ";
						HLBFile<<"  ";
						SurfTensFile<<"  ";
					}
				}
				CMCFile<<endl;
				CPFile<<endl;
				HLBFile<<endl;
				SurfTensFile<<endl;
			}
		}
	}
	queryfile.close();
	if (generateGAMSfile) GenerateParameterFileForGAMS();
}

void MoleculeQuery::GenerateParameterFileForGAMS()
{
	ofstream pFile;
	pFile.open("ParametersForGAMSProducts.txt");
	pFile<<endl;
	pFile<<"Sets"<<endl;
	pFile<<endl;

	pFile<<"i_p(i)   species that are important products /"<<endl;
	for (int i=0;i<PassedMolecules.size();i++)
		pFile<<"\'"<<Network->SMILESGAMSSpeciesMap[PassedMolecules[i]]<<"\'"<<endl;
	pFile<<"/"<<endl;

	pFile<<endl;

	/*pFile<<"j_p      product outputs /"<<endl;
	pFile<<"rxnp1*rxnp"<<PassedMolecules.size()<<endl;
	pFile<<"/"<<endl;

	pFile<<endl;*/

	/*pFile<<"A_p(i,j_p)       product stoichiometric coefficent matrix /"<<endl;
	for (int i=0;i<PassedMolecules.size();i++)
		pFile<<"\'"<<*PassedMolecules[i]<<"\' . rxnp"<<i+1<<"  -1"<<endl;
	pFile<<"/"<<endl;*/
		
		
	pFile.close();

}

//......................................................................................
//ReactionQuery begins
ReactionQuery::ReactionQuery(rxn_net_gen * Net, PathwayConstrPtr PcP,const char * file, int rule)
{
	Network = Net;
	filename = file;
	ofstream queryfile;
	string molecule_strfile(filename);
	cout<<"Reaction queries to be stored in "<<molecule_strfile<<endl;
	queryfile.open(filename);
	int count=0;
	queryfile<<"the following reactions satisfy all constraints imposed - given are reaction string, reaction number (as listed in reactions_list.txt, and heat of reaction (if available)"<<endl;
	vector<generated_rxn> OneRxn;
	if (rule <0)
	{
		for (int i=0;i<Network->AllReactions.size();i++)
		{
			if (checkRxn(Network->AllReactions[i], PcP))
			{
				queryfile<<Network->AllReactions[i].reactionstring()<<"  "<<i+1;
				if (Network->shouldCalcThermo)
				{
					double dHrxn = 0.0;
					dHrxn = Network->calculateThermoOfRxn(EnthalpyType, Network->AllReactions[i], Network->Temperature);
					queryfile<<"  "<<dHrxn;
				}
				queryfile<<endl;
			}
						
		}
	}
	else
	{
		multimap<int,int>::iterator it;
		pair<multimap<int,int>::iterator, multimap<int,int>::iterator> pair_it;
		pair_it = Network->ReactionsMap.equal_range(rule);

		for (it=pair_it.first;it!=pair_it.second;it++)
		{
			if (checkRxn(Network->AllReactions[(*it).second-1], PcP))
			{
				queryfile<<Network->AllReactions[(*it).second-1].reactionstring()<<"  "<<(*it).second;
				if (Network->shouldCalcThermo)
				{
					double dHrxn = 0.0;
					dHrxn = Network->calculateThermoOfRxn(EnthalpyType, Network->AllReactions[(*it).second-1], Network->Temperature);
					queryfile<<"  "<<dHrxn;
				}
				queryfile<<endl;
			}
		}		
	}
	queryfile.close();

}

bool ReactionQuery::checkRxn(generated_rxn Rxn, PathwayConstrPtr& PcP)
{
	vector<generated_rxn> OneRxnPath;
	OneRxnPath.push_back(Rxn);
	RxnPathway RP(&OneRxnPath);
	return (PcP(RP));
}
	



	