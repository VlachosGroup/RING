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
#include <algorithm>
#include <cmath>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include "Classheader.h"
#include "AdditionalFunctions.h"
#include "StringRegistry.h"

void rxn_net_gen::GenerateNetwork()
{
	if (!setInitialReactants()) throw 1;	
	
	vector<vector<pair<Molecule*,Patternmatch> > > ProcessedRuleVector;
	/*this vector basically stores the patterns of the molecules of the processed molecule list that are potential partners of bimolecular reaction types.
	Note that each of these rules have two reactants and hence take up two elements of the outermost vector. For e.g. if there are 10 bimolecular rules, 
	ProcessedRuleVector's size will be 20, though each of these elements could have a vector of pairs of much larger size*/
	
	int bimolcount = 0;
	for (int i=0;i<Rtlist.size();i++)
	{
		if (Rtlist[i].get_molecularity()==2)bimolcount+=2;
	}
	ProcessedRuleVector.resize(bimolcount);
	
	int number_molecules_processed = 0;
	int rxncount = 0;
	while (unprocessedmol.size()>0)
	{
		if (AllReactions.size()>rxncount*1000)
		{
			cout<<"reactions generated: "<<AllReactions.size()<<endl;
			cout<<"number of molecules to be processed yet: "<<unprocessedmol.size()<<endl;
			rxncount++;
		}
		
		/*while the unprocessed mol is not empty
		1. update the processedmol with this first molecule in the unprocessedmol.
		2. remove this molecule from unprocessedmol*/

		processedmol.push_back(unprocessedmol.begin()->second);
		unprocessedmol.erase(unprocessedmol.begin());
				
		int patt_count=0;	
		/* run through each reaction rule to generate possible reactions */
		for (int i=0;i<Rtlist.size();i++)
		{
			vector<Molecule> reactants;
			reactants.push_back(*processedmol.back());

			if(DoSimultaneousRxns) SimultRxnsGenerated.clear();//clearing up the set so that the new set of simultaneous reactions can be populated
			
			try
			{
				vector<Patternmatch> matches;
				matches.push_back(checkmatch(reactants,i,0));//first check if the molecule can be the first reactant. 
				
				/*Consider different types of rules. Case 1- rule is unimolecular. Proceed with creating a reaction object since the reactant already has the requisite pattern.*/
				
				if (Rtlist[i].get_molecularity()==1)
					GenerateMonoMolecularRxns(i, reactants, matches);
				else //for reactions that are bimolecular
					GenerateBimolecularRxns(i, reactants, matches,ProcessedRuleVector,patt_count);
			}
			catch (pair<string,string> er)
			{
				cout<<"wrong reaction rule "<<i+1<<"( "<<Rtlist[i].getRuleName()<<")"<<endl;
				cout<<"with reactants"<<endl;
				for (int si=0; si<reactants.size();si++)
				{
					cout<<reactants[si].moleculestring()<<"   ";
				}
				cout<<endl;
				cout<<"because "<<er.second<<" has a valency-bond count mismatch in atom "<<er.first<<endl;
				throw 1;
			}
		}
	
	}

	DoPONALumping();
	DoMoreLumping();
	if (LumpStrat.shoudLump()) LumpReactions();
	PrintLumps();
	printOutputInfo();
	if (shouldCalcThermo) calculateThermoValues();

	cout<<"Network generation completed"<<endl;
	
	

}

void rxn_net_gen::AddInitialReactants(std::list<string> & input)
{
	inputStrings =&input;
}

void rxn_net_gen::AddReactionRules(std::vector<Reactiontype> & rules)
{
	Rtlist = rules;
}

void rxn_net_gen::AddGlobalConstraints(ConstrPtr global)
{
	GlobalConstraints = global;
}

void rxn_net_gen::AddLumpingStrategy(LumpingStrategy& Strat)
{
	LumpStrat = Strat;
}

void rxn_net_gen::AddCompositeAtoms(vector<string>& compAtoms)
{
	for (int i=0;i<compAtoms.size();i++)
		CompositeAtomsRegistry::InsertIntoList(compAtoms[i]);
}

void rxn_net_gen::AddCompositeSites(vector<pair<string, SiteType> >& compSites)
{
	CompositeSites = compSites;
}

void rxn_net_gen::SetCalcThermo(bool YesOrNo)
{
	shouldCalcThermo = YesOrNo;
}






bool rxn_net_gen::setInitialReactants()
{
	list<string>::iterator it;
	for (it=inputStrings->begin();it!=inputStrings->end();it++)
	{
		/*for all inputs, check if it satisfies global constraints and if yes, put the molecule in hte string registry and create a new molecule object on the heap and keep the moelcule pointer.
		Update unprocessed mol with the info. The rank is set to be zero. Also start lumping if required*/
		Molecule mol ((*it),moleculesize((*it)));
		mol.unique_smiles();
		//cout<<"unique smiles is "<<mol.moleculestring()<<endl;
		pair<string, int> error = mol.checkValencyError();
		if (error.second==0)
		{
			if (globalconstraintcheck(mol))
			{
			
				string * molptr = StringRegistry::getStringPointer(mol.moleculestring());
				Molecule * molactptr = new Molecule(mol.moleculestring(), mol.getsize());
				(*molactptr).unique_smiles();
				unprocessedmol.insert(pair<int,Molecule*>(0,molactptr));
				//unprocessedmol.push_back(molactptr);
				AllMolecules.insert(pair<string*, int>(molptr, 0));
				cout<<*molptr<<endl;
				
				InitialReactants.insert(molptr);
				if (LumpStrat.shoudLump())
				{
					LumpMolecule(mol,molptr);
					
				}

			}
		}
		else
		{	
			return MolValencyErrorStatements(error);
			//if (error.second==1)cout<<"Nonaromatic atom "<<error.first<<" has a valency-bond count mismatch"<<endl;
			//if (error.second==2)cout<<"Aromatic atom "<<error.first<<" has a Hydrogen count mismatch"<<endl;
			//if (error.second==3)cout<<"Aromatic atom "<<error.first<<" has a valency-bond count mismatch"<<endl;
			//if (error.second==4)cout<<"atom N+ has a Hydrogen count mismatch"<<endl;
			//if (error.second==5)cout<<"atom N+ has a valency-bond count mismatch"<<endl;
			//if (error.second==6)cout<<"carbon atom has a Hydrogen count mismatch"<<endl;
			//return false;

		}	
	}

	return true;
}

bool rxn_net_gen::MolValencyErrorStatements(pair<string, int>& error)
{
	
	switch(error.second)
	{
		case 1:
			cout<<"Nonaromatic atom "<<error.first<<" has a valency-bond count mismatch"<<endl;
			return false;
		case 2:
			cout<<"Aromatic atom "<<error.first<<" has a Hydrogen count mismatch"<<endl;
			return false;
		case 3:
			cout<<"Aromatic atom "<<error.first<<" has a valency-bond count mismatch"<<endl;
			return false;
		case 4:
			cout<<"atom N+ has a Hydrogen count mismatch"<<endl;
			return false;
		case 5:
			cout<<"atom N+ has a valency-bond count mismatch"<<endl;
			return false;
		case 6:
			cout<<"carbon atom has a Hydrogen count mismatch"<<endl;
			return false;
		default:
			return true;
	}
}


void rxn_net_gen::DoPONALumping()
{
	//if PONA lumping is required, do it! 
	if (LumpStrat.getParaffinParameter()>=0)LumpParaffins();
	if (LumpStrat.getOlefinParameter()>=0) LumpOlefins();
	if (LumpStrat.getAromaticsParameter()>=0) LumpHydrocarbonAromatics();
	if (LumpStrat.getNaphthenicsParameter()>=0) LumpNaphthenics();
	SetPONALumpsMap();
}


		

void rxn_net_gen::printOutputInfo()
{
	cout<<"the number of species is "<<processedmol.size()<<endl;

	ofstream ruleRxnsFile;
	ruleRxnsFile.open("ReactionsOfEachRule.txt");

	for (int s=0;s<Rtlist.size();s++)
	{
		cout<<"Number of reactions with rule "<<Rtlist[s].getRuleName()<<" is "<<ReactionsMap.count(s)<<endl;
		ruleRxnsFile<<Rtlist[s].getRuleName()<<"---> ";
		pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ruleRxnsPairIter;
		ruleRxnsPairIter = ReactionsMap.equal_range(s);
		multimap<int,int>::iterator ruleRxnsIter;
		for (ruleRxnsIter = ruleRxnsPairIter.first;ruleRxnsIter!=ruleRxnsPairIter.second;ruleRxnsIter++)
			ruleRxnsFile<<ruleRxnsIter->second<<" ";
		ruleRxnsFile<<endl;
	}
	ruleRxnsFile.close();

	int number_rxns= 0;// total number - including multiplicity
	for (int s=0;s<AllReactions.size();s++)
		number_rxns+=AllReactions[s].get_frequency();

	cout<<"the number of reactions: "<<AllReactions.size()<<", with multiplicity included: "<<number_rxns<<endl;
	cout<<"the number of lumps is "<<MolLumps.size()<<endl;
	cout<<"the additional lump map size is "<<AdditionalLumpMap.size()<<endl;
	cout<<"total number of intermediates is "<<Intermediates.size()<<endl;
	map<int,int>::iterator LumpMap_it;
	set<int>UniqueAdditionalLumpSet;
	for (LumpMap_it=AdditionalLumpMap.begin();LumpMap_it!=AdditionalLumpMap.end();LumpMap_it++)
	{
		UniqueAdditionalLumpSet.insert((*LumpMap_it).second); //just adding the unique ones of the additional lumps
	}
	cout<<"the total number of lumps is "<<MolLumps.size()-AdditionalLumpMap.size()+UniqueAdditionalLumpSet.size()<<endl;
	cout<<"the number of molecules lumped is "<<MolLumpMap.size()<<endl;

	vector<pair<int,int> > LumpedRxnRuleCount;
	LumpedRxnRuleCount.resize(Rtlist.size(),pair<int,int>(0,0));
	cout<<"the number of reactions in the LumpedReactions is "<<LumpedReactionMap.size()<<endl;
	ofstream lumpedfile;
	lumpedfile.open("lumpedNetwork.txt");
	map<LumpedReaction,int,LumpedReactionCompare>::iterator lump_net_it;
	for (lump_net_it=LumpedReactionMap.begin();lump_net_it!=LumpedReactionMap.end();lump_net_it++)
	{
		LumpedRxnRuleCount[(lump_net_it->first).getRule()].first++;
		LumpedRxnRuleCount[(lump_net_it->first).getRule()].second+=lump_net_it->second;
		lumpedfile<<getRxnStringFromLumpedRxn(lump_net_it->first)<<"  "<<lump_net_it->second<<endl;
	}
	for (int s=0;s<LumpedRxnRuleCount.size();s++)
	{
		lumpedfile<<"number of occurences of rule "<<Rtlist[s].getRuleName()<<" is "<<LumpedRxnRuleCount[s].first<<" ("<<LumpedRxnRuleCount[s].second<<")"<<endl;
	}
	lumpedfile.close();


	

}

void rxn_net_gen::calculateThermoValues()
{
	vector<Molecule*>::iterator it2;
	for (it2=processedmol.begin();it2!=processedmol.end();it2++)
	{
		string * molptr = StringRegistry::getStringPointer((*it2)->moleculestring());
		if (shouldCalcThermo)
		{
			
			double dH = 0.0;
			double dS = 0.0;
			double Cp = 0.0;
			double Temp = Temperature;
			bool isIntermediateMol = false;
			if (containsSiteAtom(*molptr)) isIntermediateMol = true; 
			//cout<<"calculating thermo for "<<(*it2)->moleculestring()<<endl;
			try{
				ThermoValues T;

				if (ThermoGA::calculateDeltaH(*(*it2), isIntermediateMol, Temp,dH))
				{
					//cout<<"Values: "<<dH<<" "<<dS<<" "<<Cp<<endl;
					//system("pause");
					T.setThermo(EnthalpyType, Temp, dH);
				}
				else
					throw 1;
				if (ThermoGA::calculateDeltaS(*(*it2), isIntermediateMol, Temp,dS))
				{
					cout<<"Values: "<<dH<<" "<<dS<<" "<<Cp<<endl;
					T.setThermo(EntropyType, Temp, dS);
				}
				else
					throw 2;
				if (ThermoGA::calculateCp(*(*it2), isIntermediateMol, Temp,Cp))
					T.setThermo(CpType, Temp, Cp);
				else
					throw 3;

				T.setThermo(FreeEnergyType,Temp, dH-Temp*(dS/1000));
				
				AllMolThermo[molptr]=T;
				
			}
			catch (int er)
			{
				if (er == 1)
					cout<<"Enthalpy group additivity calculation failed for "<<(*it2)->moleculestring()<<endl;
				else if (er==2) 
					cout<<"Entropy group additivity calculation failed for "<<(*it2)->moleculestring()<<endl;				
				else cout<<"Cp group additivity calculation failed for "<<(*it2)->moleculestring()<<endl;				
				cout<<"adjacency list for the molecule is "<<endl;
				(*it2)->print_adjacency_list();
				throw 1;
			}
		}
	}
}

void rxn_net_gen::GenerateMonoMolecularRxns(int rule, std::vector<Molecule>& reactants, std::vector<Patternmatch>& matches)
{
	int rank = AllMolecules[StringRegistry::getStringPointer(reactants.front().moleculestring())];
	if (rank <= Rtlist[rule].getSpeciesRank()-1)
	{
		if (matches.front().number_of_matches()>0)
		{
			Reaction React(reactants,Rtlist[rule], matches);
			add_unique_molecules_reactions(React,rule, false);
		}
	}
}


void rxn_net_gen::GenerateBimolecularRxns(int i, std::vector<Molecule> & reactants, std::vector<Patternmatch> & matches, std::vector<vector<pair<Molecule*,Patternmatch> > > & ProcessedRuleVector, int &patt_count)
{
	
	try{	
		/*molecule already has the first pattern, now check if it also has the second and keep a note.
		This is important for A+A reactions as well as checking for rxns where this molecule is the second reactant*/
		bool MatchesFirstPattern = false;
		int rank = AllMolecules[StringRegistry::getStringPointer(reactants.front().moleculestring())];
		Patternmatch SecondPatternMatch = checkmatch(reactants,i,1);
		bool CanMatchSecondToo = false;
		if (SecondPatternMatch.number_of_matches()>0)CanMatchSecondToo = true;
			RxnsWithInterchangeableReactants.clear();
		if (matches.back().number_of_matches()>0 && !Rtlist[i].isSelfRxnOnly)
		{				
			MatchesFirstPattern = true;
			if (!Rtlist[i].IntraMolecularRxnOnly)//if the rule is not just intramolecular
			{
				ProcessedRuleVector[patt_count].push_back(pair<Molecule*,Patternmatch>(processedmol.back(),matches.back()));
				//above, the processedRuleVector is updated with info that the molecule can participate as the first reactant in the rule.
				//patt_count keeps track of the correct index
					
				if (ProcessedRuleVector[patt_count+1].size()>0)//checking if there exists potential second reactants yet
				{
					for (int it =0; it< ProcessedRuleVector[patt_count+1].size();it++)
					{
						/*if there are potential co-reactants, then for each one of those, a reaction object is created to generate all possible reactions*/
						reactants.push_back(*(ProcessedRuleVector[patt_count+1].at(it).first));
						int rank2 = AllMolecules[StringRegistry::getStringPointer(reactants.back().moleculestring())];
						if ((rank <= Rtlist[i].getSpeciesRank()-1) && (rank2 <=Rtlist[i].getSpeciesRank()-1))
						{
							if (check_combined_match(reactants,i))//if there are combined constraitns involving both reactants, they are checked here. indiv ones are checked in checkmatch
							{
								
								matches.push_back(ProcessedRuleVector[patt_count+1].at(it).second);
								Reaction React(reactants,Rtlist[i], matches);
								add_unique_molecules_reactions(React,i,CanMatchSecondToo);
								matches.pop_back();
							}
						}
						reactants.pop_back();
					}
				}
			}
		}
												
		patt_count++;

		/*having checked for all reactions wherein the molecule is first reactant, we look for those where the molecule is the second reactant. This swapping ensures that possible reactions
		with a second molecule, say M2, not yet generated, are not missed; when M2 is eventually generated, it's reaction with the current molecule will be taken care of while checking for its reactions*/
		matches.push_back(checkmatch(reactants,i,1));//find out all the matches of the second pattern in the molecule			

		if (matches.back().number_of_matches()>0 && !Rtlist[i].isSelfRxnOnly)
		{
			reactants.pop_back();//we remove the reactant so that we can swap the order of reactants and put hte molecule as second reactant
				
			if (!Rtlist[i].IntraMolecularRxnOnly)
			{
				ProcessedRuleVector[patt_count].push_back(pair<Molecule*, Patternmatch>(processedmol.back(), matches.back()));
						
				if (ProcessedRuleVector[patt_count-1].size()>0)
				{
					for (int it=0;it<ProcessedRuleVector[patt_count-1].size();it++)
					{
								
						reactants.push_back(*(ProcessedRuleVector[patt_count-1].at(it).first));//push the co-reactant as the first reactant
										
						vector<Patternmatch> matches2;
						matches2.push_back(ProcessedRuleVector[patt_count-1].at(it).second);							
						reactants.push_back(*processedmol.back());//the current molecule is the second reactant
						matches2.push_back(matches.back());
							
						int rank2 = AllMolecules[StringRegistry::getStringPointer(reactants.front().moleculestring())];//in reverse order
										
						if (rank <Rtlist[i].getSpeciesRank() && rank2 <Rtlist[i].getSpeciesRank())
						{
							if (check_combined_match(reactants,i))
							{
								//cout<<"Here3.75!!"<<endl;
								//cout<<"i: "<<i<<endl;
								//cout<<"reactants: "<<reactants.size()<<endl;
								//cout<<"matches: "<<matches2.size()<<endl;
								Reaction React(reactants,Rtlist[i], matches2);
								//cout<<"Here3.8!!"<<endl;
								add_unique_molecules_reactions(React,i,false);
							}
							/*NOTE: ProcessedRuleVector has been updated with info regarding the current molecule can participate as either of the pattern. So the case of A+A -->products is taken care of here. 
							Note that if intramolecular rxns are also allowed, then the reaction object checks for possible reactions*/
						}
						reactants.clear();
					}
				}
			}
			else
			{
				/*if the reaction rule is meant ONLY as intramolecular*/
				
				if (MatchesFirstPattern)//checks that the molecule matches both patterns!
				{
					
					reactants.push_back(*processedmol.back());
					reactants.push_back(*processedmol.back());//need to add two copies of the same reactant, only one reactant will be used eventually
									
					if (rank<Rtlist[i].getSpeciesRank())
					{
						if (check_combined_match(reactants,i))//need to check this! - do I need to do a combined check if anyways I am ignoring hte constraints?
						{
							Reaction React(reactants,Rtlist[i],matches);
							add_unique_molecules_reactions(React,i,false);
						}
					}
				}
			}			
		}
		if (Rtlist[i].isSelfRxnOnly && matches.front().number_of_matches()>0 && matches.back().number_of_matches()>0)
		{
			reactants.clear();
			reactants.push_back(*processedmol.back());//the current molecule is the second reactant
			reactants.push_back(*processedmol.back());//the current molecule is the second reactant
			
			if (rank<Rtlist[i].getSpeciesRank() && check_combined_match(reactants,i))
			{
				Reaction React(reactants,Rtlist[i],matches);
				add_unique_molecules_reactions(React,i,false);
			}		
		}

		patt_count++;
	}
	catch(pair<string,string> err)
	{
		throw;
	}
}
rxn_net_gen::rxn_net_gen(rxn_net_gen * Network, LumpingStrategy & Strat)
{
	//TODO - this is still incomplete -> the idea is to allow for pathways and mechanism be done on lumped network. 
	//This will be revisited soon.
	LumpStrat = Strat;
	Rtlist = Network->Rtlist;
	InitialReactants = Network->InitialReactants;

	if (LumpStrat.shoudLump())
	{
		
		map<string*, int, classcomp>::iterator it;
		//go over all molecules
	
		for (it=Network->AllMolecules.begin();it!=Network->AllMolecules.end();it++)
		{
			Molecule mol(*((*it).first), moleculesize(*((*it).first)));
			//lump them!
			string * molptr = StringRegistry::getStringPointer(mol.moleculestring());

			LumpMolecule(mol, molptr);

		}
		
		if (LumpStrat.getParaffinParameter()>=0)LumpParaffins();
		if (LumpStrat.getOlefinParameter()>=0) LumpOlefins();
		if (LumpStrat.getAromaticsParameter()>=0) LumpHydrocarbonAromatics();
		if (LumpStrat.getNaphthenicsParameter()>=0) LumpNaphthenics();
		for (int i=0;i<MolLumps.size();i++)
		{
			AllMolecules.insert(pair<string*,int>(MolLumps[i].getMolStringPtr(), MolLumps[i].getRank()));
			
		}

		//lump reactions and all that! 

		

		cout<<"the number of lumps is "<<MolLumps.size()<<endl;
		cout<<"the additional lump map size is "<<AdditionalLumpMap.size()<<endl;
		map<int,int>::iterator LumpMap_it;
		set<int>UniqueAdditionalLumpSet;
		for (LumpMap_it=AdditionalLumpMap.begin();LumpMap_it!=AdditionalLumpMap.end();LumpMap_it++)
		{
			UniqueAdditionalLumpSet.insert((*LumpMap_it).second);
		}
		cout<<"the total number of lumps is "<<MolLumps.size()-AdditionalLumpMap.size()+UniqueAdditionalLumpSet.size()<<endl;
		cout<<"the number of molecules lumped is "<<MolLumpMap.size()<<endl;
		PrintLumps();
	
		if (LumpStrat.shoudLump()) LumpReactions();
		//write down the lumped reactions! 
		map<LumpedReaction, int, LumpedReactionCompare>::iterator rxn_it;

		for (rxn_it=LumpedReactionMap.begin();rxn_it!=LumpedReactionMap.end();rxn_it++)
		{
			AllReactions.push_back(GetRxnFromLumpedReaction(rxn_it->first));
			for (int j=0;j<AllReactions.back().number_pdcts();j++)
			{
				MolProductMap.insert(pair<string*,int>(AllReactions.back().get_products(j),AllReactions.size()));
				
			}
			for (int j=0;j<AllReactions.back().number_reactants();j++)
			{
				MolReactantMap.insert(pair<string*,int>(AllReactions.back().get_reactants(j),AllReactions.size()));
				
			}
		}

		vector<pair<int,int> > LumpedRxnRuleCount;
	
		LumpedRxnRuleCount.resize(Rtlist.size(),pair<int,int>(0,0));
		cout<<" the number of reactions in the LumpedReactions is "<<LumpedReactionMap.size()<<endl;
		ofstream lumpedfile;
		lumpedfile.open("lumpedNetwork.txt");
		map<LumpedReaction,int,LumpedReactionCompare>::iterator lump_net_it;
		for (lump_net_it=LumpedReactionMap.begin();lump_net_it!=LumpedReactionMap.end();lump_net_it++)
		{
			LumpedRxnRuleCount[(lump_net_it->first).getRule()].first++;
			LumpedRxnRuleCount[(lump_net_it->first).getRule()].second+=lump_net_it->second;
			lumpedfile<<getRxnStringFromLumpedRxn(lump_net_it->first)<<"  "<<lump_net_it->second<<endl;
		}

		for (int s=0;s<LumpedRxnRuleCount.size();s++)
		{
			lumpedfile<<"number of occurences of rule "<<Rtlist[s].getRuleName()<<" is "<<LumpedRxnRuleCount[s].first<<" ("<<LumpedRxnRuleCount[s].second<<")"<<endl;
		}
		
	}
}


Patternmatch rxn_net_gen::checkmatch(vector<Molecule> & M, int & i, int j)
{
	//checks if reaction constraints are satisfied and then finds all patterns! 
	bool result = true;
	//Patternmatch P;
	if (j==0)
		result = (Rtlist[i].RxnConstraints.front())(M[0]);
	else 
		result = (Rtlist[i].RxnConstraints.back())(M[0]);
	if (result)
	{
		return(check_reactant_pattern(M[0],Rtlist[i].reactant_pattern[j]));
	}
	else return Patternmatch();
}
bool rxn_net_gen::check_combined_match(std::vector<Molecule> & M, int &i)
{
	//checks if combined constraints are satisfied.
	//returns true if the rule allows only/ also intramolecular reactions.. for cases of A+A reacting, I take care of it in Reactions' methods.
	if (Rtlist[i].CombinedConstraint!=NULL)
	{
		if ((Rtlist[i].isIntraMolecularAlso() || Rtlist[i].isIntraMolecularOnly()) && (M[0].moleculestring()==M[1].moleculestring()))
			return true;
		else
			return (Rtlist[i].CombinedConstraint)(M[0],M[1]);
	}	
	else return true;
	
}

Patternmatch rxn_net_gen::check_reactant_pattern(Molecule & Molec, Substructure & Sub)
{
	//cout<<Sub.fragmentstring<<endl;
	Patternmatch testmatch(Molec,Sub,0);
	testmatch.unique_matches();
	return testmatch;
}

bool rxn_net_gen::globalconstraintcheck(Molecule &Molec)
{
	return (GlobalConstraints)(Molec);
}


bool rxn_net_gen::check_product_constraints(Molecule& M, int i)
{	
	return (Rtlist[i].ProductConstraints)(M);
}




void rxn_net_gen::add_unique_molecules_reactions(Reaction & React, int RtypeIndex, bool shouldCheckMoreRxns)
{
	/*For each generated rxn
		for each product check if the global constraitns are satisfied
		if not, discard the reactions
		if yes
			check if each product molecule is new and create its lumpHash value
			if one or more of them are new, the reaction has to be new! Add the new reaction and products into the respective data structures
			if none of them are new,  generate a reaction string and check if that string for that rule already occurs
				if it occurs already, increase the frequency
				if they do not occur already, add the new reactions
	*/	
	
	vector<string> BadMolecules;//stores strings of those molecules that fail the Global constraints check! 

	map<int,int> ReactToAllReactions;//tells the index of the ith generated rxn in React in AllReactions (note some reactions in React may not be added into AllReactions!)
	
	
	for (int i=0;i<React.number_rxns_generated();i++)
	{
		/*for each reaction generated, we try to see if products satisfy global constraints, product constraints, and if they dont have valency errors (specific type pertaining to N+) and also check if they satisfy constraints on kinetics*/
		generated_rxn current_rxn(React.get_generated_rxns(i));
		current_rxn.set_rule(RtypeIndex);
		bool IsConstraintSatisfied = true;
		vector<Molecule> products;
		//cout<<"checking reactions "<<i<<endl;


		for (int j=0;j<current_rxn.number_pdcts();j++)
		{
			string* c=current_rxn.get_products(j);
			string molstring = (*c);
			Molecule mol (molstring,moleculesize(molstring));
			products.push_back(mol);
			
			pair<string,int> P = mol.checkValencyError();
			
			
		
			if (!globalconstraintcheck(products.back()) )
			{
				IsConstraintSatisfied = false;
				BadMolecules.push_back(molstring);	
				//cout<<"does not satisfy global constr"<<endl;
				break;
			}
			if (IsConstraintSatisfied)
			{
				if (!check_product_constraints(products.back(),RtypeIndex) || P.second==4 || P.second==5)
				{
					IsConstraintSatisfied = false;
					break;
				}
			}
		}

		//check for reaction k constraints! 

		if (IsConstraintSatisfied)
		{
			if (Rtlist[RtypeIndex].shouldCheckRateConstant())
			{
				//calculate kinetics and check! 

				double PreExp, ActE, TempIndex, kinValue;		
				PreExp = 0.0; ActE=0.0; TempIndex = 0.0; kinValue = 0.0;
				bool calcK = false; bool usesBEP = false; bool usesLFER = false; double alpha = 0.0; double beta = 0.0; bool stick = false;

				double dH = calculateThermoOfRxn(EnthalpyType,current_rxn,Temperature)*1000.0;
				double dS = calculateThermoOfRxn(EntropyType,current_rxn,Temperature);
				double delN = getDeltaNGasPhase(current_rxn);

				KineticsInfo kinInfo(*kinetics,current_rxn,0.0821,Temperature, dH, dS, delN, 0.0,firstGasSpecies(i,0), firstGasSpecies(i,1)); 


				kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK,usesBEP, usesLFER, alpha, beta, stick);
					
				if (kinValue < Rtlist[RtypeIndex].getMinRateConst())
					IsConstraintSatisfied = false;
			}
		}

		
		if (IsConstraintSatisfied)//if all constraints are satisfied, then we can proceed to update the rxns/ new species
		{
			bool HasNewMolecule = false;
			map<int,int> cMolecules; //to refer to closerMolecules, or parents of each product.
			
			//we assign ranks to products as one more than the largest rank of the reactants.

			
			int HighestReactantRank=-1;
			int HighestReactantRankIndex;
			bool DistinctReactantsWithEqualRanks = false;

			GetHighestRankInfoForRxn(current_rxn,HighestReactantRank, HighestReactantRankIndex, DistinctReactantsWithEqualRanks);


			/*
			string* ReactantWithHighestRank;
			
			for (int j=0;j<current_rxn.number_reactants();j++)
			{
				string* CurrentMol = current_rxn.get_reactants(j);
				int CurrentRank = AllMolecules[CurrentMol];
				if (HighestReactantRank < CurrentRank)
				{
					HighestReactantRank = CurrentRank;
					ReactantWithHighestRank = CurrentMol;
					HighestReactantRankIndex = j;
				}
				if (HighestReactantRank == CurrentRank)
				{
					if ((*ReactantWithHighestRank).compare((*CurrentMol))!=0)
					{
						DistinctReactantsWithEqualRanks = true;
					}
				}

			}*/



			for (int j=0;j<current_rxn.number_pdcts();j++)
			{

				string* c=current_rxn.get_products(j);
				GetProductParentsInfoForRxn(current_rxn,cMolecules,HighestReactantRankIndex,DistinctReactantsWithEqualRanks,j);

					
				if (AllMolecules.count(c)==0)//if the molecule is not already present in AllMolecules (the set that stores pointers to all the molecules thus encountered in the network)
				{
					
					HasNewMolecule = true;
					AllMolecules.insert(pair<string*,int>(c, HighestReactantRank+1));//update all molecules.
					
					Molecule* molptr = new Molecule(products[j]);//create a new molecule on the heap and store the pointer in unprocessedmol.
					molptr->unique_smiles();
					if (molptr->isIntermediate() || isSiteIntermediate(molptr->moleculestring())) Intermediates.insert(c);
					unprocessedmol.insert(pair<int, Molecule*>(HighestReactantRank+1, molptr));
					
					if (*c!=molptr->moleculestring())
					{
						cout<<"Oh no the strings dont match "<<*c<<"  "<<molptr->moleculestring()<<endl;
						//molptr->unique_smiles(1);
						cout<<current_rxn.reactionstring()<<"  "<<Rtlist[current_rxn.get_rule()].getRuleName()<<endl;

					}
					if (LumpStrat.shoudLump())//lump, if required.
					{
						LumpMolecule(products[j],c);//NOTE that AllMolecules has been updated --if a new lump is created, then it's rank can be specified to be equal to that of this molecule's rank
					}
					
					
				}	
				else 
				{
					//the molecule has already been encountered
					if ((HighestReactantRank+1) < AllMolecules[c])
					{
						//if this molecule is in unprocessedmol, then update it's rank;

						pair<multimap<int,Molecule*>::iterator,multimap<int,Molecule*>::iterator> iterpair;

						iterpair = unprocessedmol.equal_range(AllMolecules[c]);
						multimap<int,Molecule*>::iterator it;
						Molecule* unpmolptr;
						bool wasPresent = false;
						for (it=iterpair.first;it!=iterpair.second;it++)
						{
							if (it->second->moleculestring().compare(*c)==0)
							{
								unpmolptr = it->second;
								unprocessedmol.erase(it);
								wasPresent = true;
								break;
							}
						}
						
						if (wasPresent)unprocessedmol.insert(pair<int, Molecule*>(HighestReactantRank+1,unpmolptr));
						AllMolecules[c] = HighestReactantRank+1;//set the new, lower, rank. 
					}
				}
			}
			
			
			
			current_rxn.setParentMolecule(cMolecules);
		
			if (HasNewMolecule)
			{
				//Adding new reactions--new molecules obviously means the reaction is new!
				
				AllReactions.push_back(current_rxn);
				if (DoSimultaneousRxns)
				{
					if (Rtlist[RtypeIndex].get_molecularity()==1)SimultRxnsGenerated.insert(AllReactions.size()); // populating a set of reactions that are supposedly to happen simultaneously -- only if unimolecular
					else
					{
						set<string> reactantpairs;
						for (int r_i=0;r_i<current_rxn.number_reactants();r_i++)
							reactantpairs.insert(*current_rxn.get_reactants(r_i));
						reactantpairsForBiMol.insert(reactantpairs);
					}
				}

				AllReactions[AllReactions.size()-1].set_rule(RtypeIndex);
				AllReactions[AllReactions.size()-1].setParentMolecule(cMolecules);
				ReactionsMap.insert(pair<int, int>(RtypeIndex, AllReactions.size()));
				ReactToAllReactions.insert(pair<int,int>(i, AllReactions.size()));//updates this map 
				if (shouldCheckMoreRxns)RxnsWithInterchangeableReactants.push_back(AllReactions.size());
				
				UpdateMolAsReactantsProductsMap(current_rxn);
			}
			else
			{
				//checking if the rxn is new -> the idea is to check only with the previously generated rxn that have already been added to the AllReactions!
				//TODO - could do this without the map too
				bool IsNewRxn = true;
				string rxnstr = current_rxn.reactionstring();
				for (int k=0;k<i;k++)
				{
					if (ReactToAllReactions.count(k)>0)
					{
						//if (rxnstr == AllReactions[ReactToAllReactions[k]-1].reactionstring())
						if (current_rxn.isSameNetReaction(AllReactions[ReactToAllReactions[k]-1]))
						{
							AllReactions[ReactToAllReactions[k]-1].set_frequency(AllReactions[ReactToAllReactions[k]-1].get_frequency() /*+ current_rxn.get_frequency()*/);
							IsNewRxn = false;
							ReactToAllReactions[i] = ReactToAllReactions[k];
							break;
						}
					}
				}
				for (int k=0;k<RxnsWithInterchangeableReactants.size();k++)//checks for cases like A+ + A --> A + A+
				{
					if (current_rxn.isSameNetReaction(AllReactions[RxnsWithInterchangeableReactants[k]-1]))
					{
						AllReactions[RxnsWithInterchangeableReactants[k]-1].set_frequency(AllReactions[RxnsWithInterchangeableReactants[k]-1].get_frequency()/*+current_rxn.get_frequency()*/);
						IsNewRxn = false;
						break;
					}
				}

				if (IsNewRxn)
				{
					//add the reaction if indeed new! 
					//this piece is repetitive. could make a new function! 
					//if (*(current_rxn.get_reactants(0))=="")cout<<"empty reactant!"<<endl;
					AllReactions.push_back(current_rxn);
					if (DoSimultaneousRxns)
					{
						if (Rtlist[RtypeIndex].get_molecularity()==1)SimultRxnsGenerated.insert(AllReactions.size()); // populating a set of reactions that are supposedly to happen simultaneously	--only doing it for unimolecular ones
						else
						{
							set<string> reactantpairs;
							for (int r_i=0;r_i<current_rxn.number_reactants();r_i++)
								reactantpairs.insert(*current_rxn.get_reactants(r_i));
							reactantpairsForBiMol.insert(reactantpairs);
						}
					}

					AllReactions[AllReactions.size()-1].set_rule(RtypeIndex);
					AllReactions[AllReactions.size()-1].setParentMolecule(cMolecules);
					ReactionsMap.insert(pair<int, int>(RtypeIndex, AllReactions.size()));
					ReactToAllReactions.insert(pair<int,int>(i, AllReactions.size()));
					if (shouldCheckMoreRxns)RxnsWithInterchangeableReactants.push_back(AllReactions.size());
					
					UpdateMolAsReactantsProductsMap(current_rxn);

				}
			}
				
			
		}
		
	}
	for (int z=0;z<BadMolecules.size();z++)
	{
		StringRegistry::RemoveFromRegistry(BadMolecules[z]);
	}
	if (DoSimultaneousRxns && SimultRxnsGenerated.size()>1)//I need at least two to say they are simultaneous!
		SimultaneousRxns.push_back(SimultRxnsGenerated);
}


void rxn_net_gen::UpdateMolAsReactantsProductsMap(generated_rxn& current_rxn)
{
	set<string*,classcomp> AlreadyAdded;
	for (int k=0;k<current_rxn.number_reactants();k++)
	{
		if (AlreadyAdded.count(current_rxn.get_reactants(k))==0)
		{
			AlreadyAdded.insert(current_rxn.get_reactants(k));
			MolReactantMap.insert(pair<string*,int>(current_rxn.get_reactants(k), AllReactions.size()));
		}
	}
	AlreadyAdded.clear();

	for (int k=0;k<current_rxn.number_pdcts();k++)
	{
		if (AlreadyAdded.count(current_rxn.get_products(k))==0)
		{
			AlreadyAdded.insert(current_rxn.get_products(k));
			MolProductMap.insert(pair<string*,int>(current_rxn.get_products(k),AllReactions.size()));
		}
	}

}

void rxn_net_gen::GetHighestRankInfoForRxn(generated_rxn & current_rxn, int &HighestReactantRank, int & HighestReactantRankIndex, bool & DistinctReactantsWithEqualRanks)
{
	string* ReactantWithHighestRank;
	HighestReactantRank = -1;
			
	for (int j=0;j<current_rxn.number_reactants();j++)
	{
		string* CurrentMol = current_rxn.get_reactants(j);
		int CurrentRank = AllMolecules[CurrentMol];
		if (HighestReactantRank < CurrentRank)
		{
			HighestReactantRank = CurrentRank;
			ReactantWithHighestRank = CurrentMol;
			HighestReactantRankIndex = j;
		}
		if (HighestReactantRank == CurrentRank)
		{
			if ((*ReactantWithHighestRank).compare((*CurrentMol))!=0)
			{
				DistinctReactantsWithEqualRanks = true;
			}
		}

	}
}

void rxn_net_gen::GetProductParentsInfoForRxn(generated_rxn& current_rxn, map<int,int>& cMolecules, int HighestReactantRankIndex, bool DistinctReactantsWithEqualRanks, int j)
{
	//gets cMolecules for the reaction
	
	string* c=current_rxn.get_products(j);

	if (!DistinctReactantsWithEqualRanks)
	{
		cMolecules[j]=HighestReactantRankIndex;
	}
	else //if two distinct reactants have the highest ranks, we then find out the atom contributions of each reactant to every product and use that as a means of finding parent
	{
		vector<int> AtomsFromOtherReactants = current_rxn.getAtomsFromOthers(j);
		
		if (AtomsFromOtherReactants.empty())
		{
			//cout<<"empty AtomsFromOtherReactants!! "<<*c<<endl;
			cMolecules[j]=0;
		}
		else
		{						
			/*We need to ensure these: (a) Molecule with the highest atom contribution is a potential parent candidate
			(b) In case of equal contribution by two reactants, the first of the two reactants (in the order of specification) is preferred*/
			int MaxValue, MaxReactantIndex;
			MaxValue = 0;
			MaxReactantIndex = 0;
			int sumAtomContributions = 0;
			for (int z = 0; z < AtomsFromOtherReactants.size();z++)
			{
				sumAtomContributions+=AtomsFromOtherReactants[z];
				if (MaxValue < AtomsFromOtherReactants[z])
				{
					MaxValue = AtomsFromOtherReactants[z];
					MaxReactantIndex = z+1;//because we only include from the second reactant in AtomsFromOtherReactants
				}
			}
			int FirstReactantContribution = (moleculesize(*c) - sumAtomContributions);
			if (MaxValue <= FirstReactantContribution)
			{
				MaxValue = FirstReactantContribution;
				MaxReactantIndex = 0;
			}
			cMolecules[j] = MaxReactantIndex;

		}
	}


}

void rxn_net_gen::print_rxnlist()
{
	multimap<int, int>::iterator it;
	ofstream myfile, myfile2, reconstructfile, myfile3;
	myfile.open("reactions_SMILES.txt");
	myfile2.open("reactions_MF.txt");
	myfile3.open("reactions_list.txt");
	reconstructfile.open("reconstructed_network.txt");
	myfile<<"reactions generated are categorized into respective rules"<<endl;
	int prevRule = -1;

	cout<<"writing reactions to file"<<endl;
	for (it=ReactionsMap.begin();it!=ReactionsMap.end();it++)
	{
		
		if (it->first!=prevRule)
		{
			myfile<<"Reactions of rule ";
			myfile<<Rtlist[it->first].getRuleName()<<":"<<endl;
			prevRule = it->first;
		}
		myfile<<AllReactions[(*it).second-1].reactionstring()<<"  "<<AllReactions[(*it).second-1].get_frequency();
		myfile2<<AllReactions[(*it).second-1].get_reactionMF();
		
		if (shouldCalcThermo) myfile<<"   "<<calculateThermoOfRxn(EnthalpyType,AllReactions[(*it).second-1], Temperature)<<"   "<<calculateThermoOfRxn(EntropyType,AllReactions[(*it).second-1],Temperature)<<"   "<<calculateThermoOfRxn(FreeEnergyType, AllReactions[(*it).second-1],Temperature);
		myfile<<endl;
		myfile2<<endl;
	}
	for (int i=0;i<ReconstructedReactions.size();i++)
	{
		reconstructfile<<ReconstructedReactions[i].reactionstring()<<endl;
	}
	reconstructfile<<"the size of the reconstructed network is "<<ReconstructedReactions.size()<<endl;
	reconstructfile.close();

	myfile3<<"Reactions as they were generated are listed below"<<endl;
	for (int i=0;i<AllReactions.size();i++)
		myfile3<<AllReactions[i].reactionstring()<<endl;


	myfile.close();
	myfile2.close();
	myfile3.close();

	vector<Molecule*>::iterator it2;
	ofstream speciesFile;
	speciesFile.open("species.txt");
	speciesFile<<"This file lists all species in the file as SMILES strings, followed by size, rank, and Thermo values if specied "<<endl;
	if (shouldCalcThermo) speciesFile<<"The order for thermochemistry data is: Enthalpy (kJ/mol), Entropy (J/mol/K), Free Energy(kJ/mol), and cp (J/mol/K) at system temperature (and within parenthesis at different temperatures: 300K, 400K, 500K, 600K, 700K, 800K, 1000K, 1500K)"<<endl; 
	for (it2=processedmol.begin();it2!=processedmol.end();it2++)
	{
		string * molptr = StringRegistry::getStringPointer((*it2)->moleculestring());
		speciesFile<<(*it2)->moleculestring()<<"  "<<(*it2)->getsize()<<"  "<<AllMolecules[molptr]<<"  ";
		if (shouldCalcThermo)
		{
			speciesFile<<AllMolThermo[molptr].getThermo(EnthalpyType, Temperature)<<"  "<<AllMolThermo[molptr].getThermo(EntropyType, Temperature)<<"  "<<AllMolThermo[molptr].getThermo(FreeEnergyType, Temperature)<<"  "<<AllMolThermo[molptr].getThermo(CpType,Temperature)<<" (";
			vector<int> cpTemps;
			cpTemps.push_back(300);cpTemps.push_back(400);cpTemps.push_back(500);cpTemps.push_back(600);cpTemps.push_back(700);cpTemps.push_back(800);cpTemps.push_back(1000);cpTemps.push_back(1500);

			bool isIntermediateMol = false;
			if (containsSiteAtom(*molptr)) isIntermediateMol = true; 
			for (int i = 0;i<cpTemps.size();i++)
			{

				double cp = 0.0;
				ThermoGA::calculateCp(*(*it2), isIntermediateMol, cpTemps[i],cp);
				speciesFile<<"  "<<cp;
			}
			speciesFile<<" )"<<endl;
		}
		else
			speciesFile<<endl;
				


	}
	speciesFile.close();

	cout<<"Completed writing reactions and species lists into file"<<endl;
}


int rxn_net_gen::GetAtomValue(string S)
{
	if (S.compare(0,1,"C")==0)return 7;
	else if (S.compare(0,1,"H")==0)return 5;
	else if (S.compare(0,1,"O")==0)return 11;
	else if (S.compare(0,1,"N")==0)return 13;
	else if (S.compare(0,1,"S")==0)return 17;
	else if (S.compare(0,1,"P")==0)return 19;
	else 
	{
		int temp = CompositeAtomsRegistry::getIndexOfAtom("{"+ S + "}");
		if (temp!=-1)
			return (prime(6+temp));
		else return 0;
	}
}


pair<unsigned int, unsigned int> rxn_net_gen::CreateHashValue(const Molecule& M)
{
	unsigned int Hash = 1;
	unsigned int ElectronicHash = 1;
	if (!M.HasNBinteractions())
	{
		for (int i=0;i<M.getsize();i++)
		{
			int AtomtypeIndexCalc=1;
			if (M.getatom(i)->get_nature()==0)//if single atom
				AtomtypeIndexCalc=GetAtomValue(M.getatomtype(i));
			else AtomtypeIndexCalc = 1;
			
			if (AtomtypeIndexCalc==0)
			{
				Hash=0;
				break;
			}
			else
			{

				//if aromatic or ring atom
				if (M.isaromatic(i))AtomtypeIndexCalc=AtomtypeIndexCalc*3;
				if (M.isringatom(i))AtomtypeIndexCalc=AtomtypeIndexCalc*9;//multiplying by 9 is still fine, the product will be a multiple of 9 only if it is ring atom! 

				//neighbors
				for (int j=0;j<M.getNN(i);j++)
				{
					AtomtypeIndexCalc=AtomtypeIndexCalc*2;
					int atomvalue=GetAtomValue(M.getatomtype(M.get_adjacency(i,j)));
					int bovalue = M.get_BO(i,j);
					for (int k=0;k<bovalue;k++)
					{
						AtomtypeIndexCalc=AtomtypeIndexCalc*atomvalue;
					}
					if (M.isringatom(M.get_adjacency(i,j)))AtomtypeIndexCalc=AtomtypeIndexCalc*9;//multiplying by 9 here -- not doing for aromatic atoms - because they can be distinguished in other ways also! 


				
				}

			/*	for (int j=0;j<M.getHydrogens(i);j++)
				{
					AtomtypeIndexCalc=AtomtypeIndexCalc*5;
				}*/
				
				if (AtomtypeIndexCalc!=0)
				{
					map<int,int>::iterator it;
					it =AtomtypeIndex.find(AtomtypeIndexCalc);
					if (it!=AtomtypeIndex.end())
					{
						Hash=Hash*prime((*it).second);
						int eh = M.getElectronicHashValue(i);
						if (eh!=1)
							ElectronicHash=ElectronicHash*prime((*it).second)*eh;//note that the prime number as well as 'eh' could both be 2,3,5, or 7. But, that is still fine because say if it were both 2, then the product 4 is only because of that particular combination 
						//this is so because in each case (the prime number corresponding to the atomtype) and 'eh', the index of a prime number (in their prime factorization) is always one!

					}
					else
					{
						AtomtypeIndex.insert(pair<int,int>(AtomtypeIndexCalc,AtomtypeIndex.size()));
						Hash=Hash*prime(AtomtypeIndex.size()-1);
						int eh = M.getElectronicHashValue(i);
						if (eh!=1)
							ElectronicHash=ElectronicHash*prime(AtomtypeIndex.size()-1)*eh;
					}
				}
				else 
				{
					Hash = 0;
					break;
				}
			}
		}
	}
	else Hash=0;
	//if (M.moleculestring()=="[C+]" || M.moleculestring()=="[H+]")cout<<M.moleculestring()<<" "<<Hash<<endl;
	return pair<unsigned int, unsigned int>(Hash,ElectronicHash);
}

generated_rxn rxn_net_gen::getReaction(int i)
{
	return AllReactions[i];
}

void rxn_net_gen::LumpMolecule(const Molecule& mol, string* molptr)
{
		
	//lumps a molecule according to its hash value, and lumping strategy in general.
	//Initial reactants are not lumped
	UnsignedIntPair hashvalue = CreateHashValue(mol);//creates hash value on which lumping is based
	bool foundLump = false;
	/*if (*molptr=="C1(C=C)(C(C=CC1C)C)C")
		cout<<"C1(C=C)(C(C=CC1C)C)C here: "<<hashvalue.first<<"  "<<hashvalue.second<<endl;
	if (*molptr =="C(=C)C")
		cout<<"it's propene"<<hashvalue.first<<"  "<<hashvalue.second<<endl;

	*/
	if (hashvalue.first>0)
	{
		//Initial reactants are not lumped! 
		//functional lumping constraints must also be satisfied

		bool anyFncConstraints = LumpStrat.isThereFunctionalConstraints();
		bool passesConstraints = (anyFncConstraints && (LumpStrat.getFunctionalLumpingConstraints())(mol)) || !anyFncConstraints;

		
		if (InitialReactants.count(molptr)==0 && LumpHashMap.count(hashvalue)>0 && passesConstraints)//checking if a corresponding lump has already been created
		{
			//if (*molptr=="C1(C=C)(C(C=CC1C)C)C")cout<<"found a lump!"<<endl;

			
			multimap<UnsignedIntPair,int>::iterator iter1;
			pair<multimap<UnsignedIntPair,int>::iterator, multimap<UnsignedIntPair, int>::iterator> pairiter;
			pairiter = LumpHashMap.equal_range(hashvalue);//gets the range of all lumps having the same hash value --> usually only one, rarely may collide
			for (iter1=pairiter.first;iter1!=pairiter.second;iter1++)
			{
				if (MolLumps[iter1->second].getSize()==mol.getsize() && InitialReactants.count(MolLumps[iter1->second].getMolStringPtr())==0)//do a second round of check --molecule size check is good enough at this point! But also include the check to see if this lump does not correspond to the initial reactants
				{
					foundLump = true;
					MolLumpMap.insert(pair<string*, int>(molptr, iter1->second));//updated what lump a mol belongs to
					pair<int,int> lumpvalues = mol.calcBranchandDistanceValue();
					
					//checking if lump representative parameters must be changed, depending on whether molecule is acyclic or cyclic and the lump strategy instructions given.
					if (mol.iscyclicmolecule())
					{
						if ((LumpStrat.getringParameter()==0 && MolLumps[iter1->second].getringDistanceValue()> lumpvalues.second)
						||(LumpStrat.getringParameter()==1 && MolLumps[iter1->second].getringDistanceValue()< lumpvalues.second))
						{
							//change lump properties
							MolLumps[iter1->second].setringDistanceValue(lumpvalues.second);
							MolLumps[iter1->second].setMoleculeString(molptr); //sets molecule string to the new representative 
						}
					}
					
					if ((LumpStrat.getchainParameter()==0 && MolLumps[iter1->second].getBranchValue()> lumpvalues.first && MolLumps[iter1->second].getringDistanceValue() == lumpvalues.second)
					||(LumpStrat.getchainParameter()==1 && MolLumps[iter1->second].getBranchValue()<lumpvalues.first && MolLumps[iter1->second].getringDistanceValue() == lumpvalues.second))
					{
						MolLumps[iter1->second].setMoleculeString(molptr);
						MolLumps[iter1->second].setBranchValue(lumpvalues.first);
					}
					
					//adding molecule rank - rank of the lump is the smallest of all the constituents
					
					if (MolLumps[iter1->second].getRank()>AllMolecules[molptr])MolLumps[iter1->second].setRank(AllMolecules[molptr]);
					
					break;
				}
			}
		}
		if (!foundLump)//if the lump is not found, then we need to create the new lump and set all parameters! 
		{
			pair<int,int> valuepair = mol.calcBranchandDistanceValue();
			int hydrogens=0;

			//if (*molptr=="C1(C=C)(C(C=CC1C)C)C")cout<<"not found lump, so creating a new one!"<<endl;

			for (int z=0;z<mol.getsize();z++)
			{
				hydrogens+=mol.getHydrogens(z);
			}
			LumpInfo L (mol.getsize(), hydrogens,mol.iscyclicmolecule(), valuepair.first, valuepair.second, molptr);
			L.setleaves(mol.NumberOfLeaves());
			L.setRank(AllMolecules[molptr]); //molptr is already defined within AllMolecules!
			L.setPONAcharacteristic(-1);
			L.setBondCounts(mol.totalDoubleBonds(), mol.totalTripleBonds(), mol.totalAromaticBonds());
			L.setMaxRingSize(mol.MaxRingsize());
			L.setRingCount(mol.NumberOfRings());

			if (LumpStrat.getParaffinParameter()>=0)
			{
				if (mol.isparaffinic()) //this ensures that this is acyclic - so naphthenics are not included here
				{
					L.setPONAcharacteristic(0);
					L.setPONAElectronicValue(mol.getTotalElectronicValue());
				}
			}
			if (LumpStrat.getOlefinParameter()>=0)
			{
				if (mol.isolefinic())//this ensures that this is acyclic also so napthenics are not included here
				{
					L.setPONAcharacteristic(1);
					L.setPONAElectronicValue(mol.getTotalElectronicValue());
				}
			}
			if (LumpStrat.getNaphthenicsParameter()>=0)
			{
				if (mol.isNaphthenic())//ensures this is not aromatic. so aromatics are not included here
				{
					L.setPONAcharacteristic(2);
					L.setPONAElectronicValue(mol.getTotalElectronicValue());
				}
			}
			if (LumpStrat.getAromaticsParameter()>=0)
			{
				if (mol.ishydrocarbon() && mol.isaromaticmolecule())
				{
					if (mol.totalDoubleBonds()==0)L.setPONAcharacteristic(3);
					else L.setPONAcharacteristic(4);
					L.setPONAElectronicValue(mol.getTotalElectronicValue());
				}
			}
			string elecformula ="";
			MolLumps.push_back(L);
			LumpHashMap.insert(pair<pair<unsigned int,unsigned int>,int>(hashvalue,MolLumps.size()-1));
			MolLumpMap.insert(pair<string*,int>(molptr, MolLumps.size()-1));
			
		}
	}
}

void rxn_net_gen::DoMoreLumping()
{
	vector<ConstrPtr> MoreAdditionalLumpingConstraints = LumpStrat.getMoreLumpingConstraints();
	vector<Molecule> LumpMolecules;
			
	if (LumpStrat.isThereMoreLumpingConstraints())
	{
		for (int i=0;i<MolLumps.size();i++)
		{
			string molstring = MolLumps[i].getMoleculeString();
			Molecule mol(molstring, moleculesize(molstring));
			mol.unique_smiles();

			LumpMolecules.push_back(mol);

			vector<int> MoreAdditionalLumpingParameters  = LumpStrat.getMoreLumpingParameter();

			for (int j = 0;j<MoreAdditionalLumpingConstraints.size();j++)
			{
				if ((MoreAdditionalLumpingConstraints[j])(mol))//if molecule satisfies the prescribed constraints for more lumping
				{		
					//go ahead and lump this molecule only if it has not already been covered by PONA! 
					bool goAheadAndLump = true;
					if (LumpStrat.getParaffinParameter()>=0 && mol.isparaffinic())goAheadAndLump=false;
					if (LumpStrat.getOlefinParameter()>=0 && mol.isolefinic()) goAheadAndLump =false;
					if (LumpStrat.getAromaticsParameter()>=0 && mol.isaromaticmolecule()) goAheadAndLump = false;
					if (LumpStrat.getNaphthenicsParameter()>=0 && mol.isNaphthenic()) goAheadAndLump = false;
					

					if (goAheadAndLump) 
					{
						//define P that has the PONALumpInfo 	
						MoreAdditionalLumpingInfo P;

						P.Size = MolLumps[i].getSize(); P.HydCount = MolLumps[i].getHydrogens();P.DoubleBonds = MolLumps[i].getDoubleBondCount();
						P.TripleBonds = MolLumps[i].getTripleBondCount(); P.NumberOfRings = MolLumps[i].getRingCount(); P.MaxRingSize = MolLumps[i].getMaxRingSize();
						P.AromaticBonds = MolLumps[i].getAromaticBondCount(); P.WhichMoreLumpingDcl = j;

						//MoreLumpingInfo P(MolLumps[i].getSize(),pair<int,pair<int,int> >(MolLumps[i].getHydrogens(), pair<int,int> (MolLumps[i].getScore(), MolLumps[i].getPONAElectronicValue())));

						
						if (MoreAdditionalLumpingSizeLumpMap.count(P)>0) //if already present in MoreLumpingSizeLumpMap
						{
							//checking the leaves in the lump to do additional lumping to the most or least branched
				
							if ((MoreAdditionalLumpingParameters[j]==0 && MolLumps[i].getleaves() < MolLumps[MoreAdditionalLumpingSizeLumpMap[P]].getleaves())
								|| (MoreAdditionalLumpingParameters[j]==1 && MolLumps[i].getleaves() > MolLumps[MoreAdditionalLumpingSizeLumpMap[P]].getleaves()))	
								MoreAdditionalLumpingSizeLumpMap[P] = i;//setting the representative paraffin to this lump if it satisfies the branching requirements
						}
						else MoreAdditionalLumpingSizeLumpMap[P] = i;//else set MoreLumpingSizeLumpMap
						break;//if a lump has been identified or new lump created, then work's done, break out of the inner for loop
					}	
				}
			}
					
		}
	}

	//updating additional lump  map -- need to do this at the end because the actual representative lump can change! 
	for (int i=0;i< MolLumps.size();i++)
	{
		//MoreLumpingInfo P(MolLumps[i].getSize(),pair<int,pair<int,int> >(MolLumps[i].getHydrogens(), pair<int,int> (MolLumps[i].getScore(), MolLumps[i].getPONAElectronicValue())));

		MoreAdditionalLumpingInfo P;

		P.Size = MolLumps[i].getSize(); P.HydCount = MolLumps[i].getHydrogens();P.DoubleBonds = MolLumps[i].getDoubleBondCount();
		P.TripleBonds = MolLumps[i].getTripleBondCount(); P.NumberOfRings = MolLumps[i].getRingCount(); P.MaxRingSize = MolLumps[i].getMaxRingSize();
		P.AromaticBonds = MolLumps[i].getAromaticBondCount();
		P.WhichMoreLumpingDcl=-1;

		for (int j = 0;j<MoreAdditionalLumpingConstraints.size();j++)
		{
			if ((MoreAdditionalLumpingConstraints[j])(LumpMolecules[i]))//if molecule satisfies the prescribed constraints for more lumping
			{
				P.WhichMoreLumpingDcl=j;
				break;
			}
		}

		if (MoreAdditionalLumpingSizeLumpMap.count(P)>0) //this is similar to what I do in SetPONALumpsMap! 
			AdditionalLumpMap[i] = MoreAdditionalLumpingSizeLumpMap[P];
	}
}

void rxn_net_gen::LumpParaffins()
{
	//Initial reactants are NOT lumped.
	for (int i=0;i< MolLumps.size();i++)
	{
		if (MolLumps[i].getPONAcharacteristic()==0 && InitialReactants.count(MolLumps[i].getMolStringPtr())==0)
		{
			string molstring = MolLumps[i].getMoleculeString();

			Molecule mol(molstring, moleculesize(molstring));
			mol.unique_smiles();

			
			

			bool AnyParaffinConstraints = LumpStrat.isThereParaffinConstraints();
			bool goAheadAndLump = (AnyParaffinConstraints && (LumpStrat.getParaffinConstraints())(mol)) || !AnyParaffinConstraints; //does molecule satisfy paraffin lumping constraints if there's one? 


			if (goAheadAndLump) 
			{
				ToUndergoAdditionalLumping.insert(i);
				//define P that has the PONALumpInfo 			
				PONALumpingInfo P(pair<int,int>(MolLumps[i].getSize(),MolLumps[i].getHydrogens()),MolLumps[i].getPONAElectronicValue());
			
				if (mol.getsize()==4)
					cout<<"go ahead and lump "<<mol.moleculestring()<<endl;
				if (ParaffinSizeLumpMap.count(P)>0) //if already present in ParaffinSizeLumpMap
				{
					//checking the leaves in the lump to do additional lumping to the most or least branched
			
					if ((LumpStrat.getParaffinParameter()==0 && MolLumps[i].getleaves() < MolLumps[ParaffinSizeLumpMap[P]].getleaves())
						|| (LumpStrat.getParaffinParameter()==1 && MolLumps[i].getleaves() > MolLumps[ParaffinSizeLumpMap[P]].getleaves()))
						ParaffinSizeLumpMap[P] = i;//setting the representative paraffin to this lump if it satisfies the branching requirements
				}
				else ParaffinSizeLumpMap[P] = i;//else set ParaffinSizeLumpMap. 
				
				//NOTE: I haven't yet noted down anywhere that this lump is to undergo additional lump. I need to say which representative lump it needs to be lumped to. 
				//This is done in setPONAlumps because I cannot do it here as I don't know what my representative is yet. 
				
				//AdditionalLumpMap[i] = ParaffinSizeLumpMap[P];
			}
			
		}
	}
	
}
void rxn_net_gen::LumpOlefins()
{
	
	//Initial reactants are NOT lumped
	for (int i=0;i< MolLumps.size();i++)
	{
		if (MolLumps[i].getPONAcharacteristic()==1 && InitialReactants.count(MolLumps[i].getMolStringPtr())==0)
		{
			string molstring = MolLumps[i].getMoleculeString();

			Molecule mol(molstring, moleculesize(molstring));
			mol.unique_smiles();

			bool AnyOlefinConstraints = LumpStrat.isThereOlefinConstraints();
			bool goAheadAndLump = (AnyOlefinConstraints && (LumpStrat.getOlefinConstraints())(mol)) || !AnyOlefinConstraints; //does molecule satisfy olefin lumping constraints if there's one? 

			//if (mol.moleculestring()=="C1(C=C)(C(C=CC1C)C)C")
			//	cout<<"Oh no! it's classified as olefins"<<endl;
			if (goAheadAndLump) 
			{
			
				ToUndergoAdditionalLumping.insert(i);
				PONALumpingInfo P(pair<int,int>(MolLumps[i].getSize(), MolLumps[i].getHydrogens()),MolLumps[i].getPONAElectronicValue());
				

				if (OlefinSizeLumpMap.count(P)>0)
				{
					if ((LumpStrat.getOlefinParameter()==0 && MolLumps[i].getleaves() < MolLumps[OlefinSizeLumpMap[P]].getleaves())
						|| (LumpStrat.getOlefinParameter()==1 && MolLumps[i].getleaves() > MolLumps[OlefinSizeLumpMap[P]].getleaves()))
						OlefinSizeLumpMap[P] = i;//checking the leaves like in LumpOlefins
				}
				else OlefinSizeLumpMap[P] = i;
				//NOTE: I haven't yet noted down anywhere that this lump is to undergo additional lump. I need to say which representative lump it needs to be lumped to. 
				//This is done in setPONAlumps because I cannot do it here as I don't know what my representative is yet. 
				
				//AdditionalLumpMap[i] = OlefinSizeLumpMap[P];
			}
		}
	}
}
void rxn_net_gen::LumpNaphthenics()
{
	//Initial reactants are not lumped
	ofstream naphthenicsfile("naphthenicsfile.txt", ios::app);

	for (int i=0;i< MolLumps.size();i++)
	{
		if (MolLumps[i].getPONAcharacteristic()==2 && InitialReactants.count(MolLumps[i].getMolStringPtr())==0)
		{
			string molstring = MolLumps[i].getMoleculeString();

			Molecule mol(molstring, moleculesize(molstring));
			mol.unique_smiles();

			//if (mol.moleculestring()=="C1(=CC(CCC1)C)C")
				//cout<<"phew! it is under naphthenics "<<endl;
			//if (mol.moleculestring()=="C(=C)C")
				//cout<<"propene is here! Oh No!"<<endl;

			


			bool AnyNaphthenicsConstraints = LumpStrat.isThereNaphthenicsConstraints();
			bool goAheadAndLump = (AnyNaphthenicsConstraints && (LumpStrat.getNaphthenicsConstraints())(mol)) || !AnyNaphthenicsConstraints; //does molecule satisfy naphthenics lumping constraints if there's one? 


			if (goAheadAndLump) 
			{

				ToUndergoAdditionalLumping.insert(i);
				PONALumpingInfo P(pair<int,int>(MolLumps[i].getSize(), MolLumps[i].getHydrogens()),MolLumps[i].getPONAElectronicValue()*(int)(pow(11.0,mol.MaxRingsize()-2)));
				/*NOTE and TODO: the last value that typically only stores the electronic value has been modified to include a score for max number of rings as 2^maxRings -- this is to ensure that, for eg., five and six memebered rings are differentiated. */
			/*	if (mol.moleculestring()=="C1(=CC(CCC1)C)C")
				{
					cout<<"molecule is C1(=CC(CCC1)C)C "<<MolLumps[i].getSize()<<" "<<MolLumps[i].getHydrogens()<<"  "<<MolLumps[i].getPONAElectronicValue()<<"  "<<(int)(pow(11.0,mol.MaxRingsize()-2))<<endl;
					naphthenicsfile<<"molecule is C1(=CC(CCC1)C)C "<<MolLumps[i].getSize()<<" "<<MolLumps[i].getHydrogens()<<"  "<<MolLumps[i].getPONAElectronicValue()<<"  "<<(int)(pow(11.0,mol.MaxRingsize()-2))<<endl;
				}
				if (mol.moleculestring()=="C1(=C)CCC(CC1)C")
				{
					cout<<"molecule is C1(=C)CCC(CC1)C "<<MolLumps[i].getSize()<<" "<<MolLumps[i].getHydrogens()<<"  "<<MolLumps[i].getPONAElectronicValue()<<"  "<<(int)(pow(11.0,mol.MaxRingsize()-2))<<endl;
					naphthenicsfile<<"molecule is C1(=C)CCC(CC1)C "<<MolLumps[i].getSize()<<" "<<MolLumps[i].getHydrogens()<<"  "<<MolLumps[i].getPONAElectronicValue()<<"  "<<(int)(pow(11.0,mol.MaxRingsize()-2))<<endl;
				}*/
				
				
				if (NapthenicsSizeLumpMap.count(P)>0)
				{
					if ((LumpStrat.getNaphthenicsParameter()==0 && MolLumps[i].getleaves() < MolLumps[NapthenicsSizeLumpMap[P]].getleaves())
						|| (LumpStrat.getNaphthenicsParameter()==1 && MolLumps[i].getleaves() > MolLumps[NapthenicsSizeLumpMap[P]].getleaves()))
						NapthenicsSizeLumpMap[P] = i; //like in LumpParaffins/Olefins
				}
				else NapthenicsSizeLumpMap[P] = i;
				//NOTE: I haven't yet noted down anywhere that this lump is to undergo additional lump. I need to say which representative lump it needs to be lumped to. 
				//This is done in setPONAlumps because I cannot do it here as I don't know what my representative is yet. 
				
				//AdditionalLumpMap[i] = NapthenicsSizeLumpMap[P];
			}
		}
	}
}
void rxn_net_gen::LumpHydrocarbonAromatics()
{
	//Initial reactants are NOT lumped
	for (int i=0;i< MolLumps.size();i++)
	{
		string molstring = MolLumps[i].getMoleculeString();

		Molecule mol(molstring, moleculesize(molstring));
		mol.unique_smiles();

		bool AnyAromaticsConstraints = LumpStrat.isThereAromaticsConstraints();
		bool goAheadAndLump = (AnyAromaticsConstraints && (LumpStrat.getAromaticsConstraints())(mol)) || !AnyAromaticsConstraints; //does molecule satisfy aromatics lumping constraints if there's one? 

		//if (mol.moleculestring()=="C1(C=C)(C(C=CC1C)C)C")
				//cout<<"Oh no! it's classified as aromatics"<<endl;

		if (goAheadAndLump) 
		{
			ToUndergoAdditionalLumping.insert(i);
			//TODO: I don't remember why I call the two maps AromaticsBranchLumpMap and not AromaticsSizeLumpMap --perhaps because the leaves directly correspond to alkyl groups?!

			PONALumpingInfo P(pair<int,int>(MolLumps[i].getSize(), MolLumps[i].getHydrogens()),MolLumps[i].getPONAElectronicValue());
			//pair<int,int> P(MolLumps[i].getSize(),MolLumps[i].getPONAElectronicValue());
			if (MolLumps[i].getPONAcharacteristic()==3 && InitialReactants.count(MolLumps[i].getMolStringPtr())==0)
			{
				if (AlkylAromaticsBranchLumpMap.count(P)>0)
				{
					if ((LumpStrat.getAromaticsParameter()==0 && MolLumps[i].getleaves() < MolLumps[AlkylAromaticsBranchLumpMap[P]].getleaves())
						|| (LumpStrat.getAromaticsParameter()==1 && MolLumps[i].getleaves() > MolLumps[AlkylAromaticsBranchLumpMap[P]].getleaves()))
						AlkylAromaticsBranchLumpMap[P] = i;// like in LumpParaffins/Olefins/Naphthenics
				}
				else AlkylAromaticsBranchLumpMap[P] = i;
				//AdditionalLumpMap[i] = AlkylAromaticsBranchLumpMap[P];

				//NOTE: I haven't yet noted down anywhere that this lump is to undergo additional lump. I need to say which representative lump it needs to be lumped to. 
				//This is done in setPONAlumps because I cannot do it here as I don't know what my representative is yet. 
				
			}
			if (MolLumps[i].getPONAcharacteristic()==4 && InitialReactants.count(MolLumps[i].getMolStringPtr())==0)
			{
				if (AlkenylAromaticsBranchLumpMap.count(P)>0)
				{
					if ((LumpStrat.getAromaticsParameter()==0 && MolLumps[i].getleaves() < MolLumps[AlkenylAromaticsBranchLumpMap[P]].getleaves())
						|| (LumpStrat.getAromaticsParameter()==1 && MolLumps[i].getleaves() > MolLumps[AlkenylAromaticsBranchLumpMap[P]].getleaves()))
						AlkenylAromaticsBranchLumpMap[P] = i;// like in LumpParaffins/Olefins/Naphthenics
				}
				else AlkenylAromaticsBranchLumpMap[P] = i;
				//AdditionalLumpMap[i] = AlkenylAromaticsBranchLumpMap[P];

				//NOTE: I haven't yet noted down anywhere that this lump is to undergo additional lump. I need to say which representative lump it needs to be lumped to. 
				//This is done in setPONAlumps because I cannot do it here as I don't know what my representative is yet. 
				
			}
		}
	}
}

void rxn_net_gen::SetPONALumpsMap()
{
	for (int i=0;i< MolLumps.size();i++)
	{
		
		if (ToUndergoAdditionalLumping.count(i)>0)
		{
			string molstring = MolLumps[i].getMoleculeString();

			Molecule mol(molstring, moleculesize(molstring));
			mol.unique_smiles();
			PONALumpingInfo P(pair<int,int>(MolLumps[i].getSize(), MolLumps[i].getHydrogens()),MolLumps[i].getPONAElectronicValue());
			
			if (MolLumps[i].getPONAcharacteristic()==0 && ParaffinSizeLumpMap.count(P)>0)
				AdditionalLumpMap[i] = ParaffinSizeLumpMap[P];
			else if (MolLumps[i].getPONAcharacteristic()==1 && OlefinSizeLumpMap.count(P)>0)
				AdditionalLumpMap[i] = OlefinSizeLumpMap[P];
			else if (MolLumps[i].getPONAcharacteristic()==2)
			{
				P.second*=(int)(pow(11.0,mol.MaxRingsize()-2));
				if (NapthenicsSizeLumpMap.count(P)>0)
					AdditionalLumpMap[i] = NapthenicsSizeLumpMap[P];
			}
			else if (MolLumps[i].getPONAcharacteristic()==3 && AlkylAromaticsBranchLumpMap.count(P)>0)
				AdditionalLumpMap[i] = AlkylAromaticsBranchLumpMap[P];
			else if (MolLumps[i].getPONAcharacteristic()==4 && AlkenylAromaticsBranchLumpMap.count(P)>0)
				AdditionalLumpMap[i] = AlkenylAromaticsBranchLumpMap[P];

		}
	}
}





LumpedReaction rxn_net_gen::GenerateLumpedReaction(generated_rxn g)
{
	vector<int> reactants;
	vector<int> products;
	for (int i=0;i<g.number_reactants();i++)
	{
		int lumpIndex = MolLumpMap[g.get_reactants(i)];
		if (AdditionalLumpMap.count(lumpIndex)>0)reactants.push_back(AdditionalLumpMap[lumpIndex]);
		else reactants.push_back(lumpIndex);
	}
	for (int i=0;i<g.number_pdcts();i++)
	{
		int lumpIndex = MolLumpMap[g.get_products(i)];
		if (AdditionalLumpMap.count(lumpIndex)>0)products.push_back(AdditionalLumpMap[lumpIndex]);
		else products.push_back(lumpIndex);
	}
	
	return LumpedReaction(reactants, products, g.get_rule());
}

generated_rxn rxn_net_gen::GetRxnFromLumpedReaction(const LumpedReaction& L)
{
	//NOTE: this generates a generated_rxn object from a lump. Note that the order of reactant may not be preserved. 
	// Also 
	
	generated_rxn G;

	vector<int> reactants = L.getReactantLumps();
	vector<int> products = L.getProductLumps();

	
	
	for (int i =0;i<reactants.size();i++)
	{
		G.add_reactants(MolLumps[reactants[i]].getMolStringPtr());
	}
	for (int i =0;i<products.size();i++)
	{
		G.add_products(MolLumps[products[i]].getMolStringPtr());
	}
	G.set_rule(L.getRule());

	//add everything that I would add when I generate a new reaction! */

	return G;
}
	

	

void rxn_net_gen::LumpReactions()
{
	map<LumpedReaction,int,LumpedReactionCompare> HowManyLumped;

	for (int i=0;i<AllReactions.size();i++)
	{
		
		LumpedReaction L = GenerateLumpedReaction(AllReactions[i]);

		if (LumpedReactionMap.count(L)>0)
		{
			LumpedReactionMap[L]+=AllReactions[i].get_frequency();
			HowManyLumped[L]+=1;
		}
		else 
		{
			LumpedReactionMap[L]=AllReactions[i].get_frequency();
			HowManyLumped[L] = 1;
		}
	}

	//now set an average frequency -- this is because adding up creates a problem with calculating kinetics.
	//note that lumping implies that the concentration is already higher - so an average is a much better way to keep track of symmetry numbers
	//else we will be double counting the factor coming from symmetry

	map<LumpedReaction,int,LumpedReactionCompare>::iterator it;

	for (it=LumpedReactionMap.begin();it!=LumpedReactionMap.end();it++)
	{
		it->second=it->second/HowManyLumped[it->first];
	}




}


int rxn_net_gen::findLumpofMolecule(string* S)
{
	int lump = MolLumpMap[S];
	if (AdditionalLumpMap.count(lump)>0)
		lump=AdditionalLumpMap[lump];
	return lump;
}


void rxn_net_gen::PrintLumps()
{
	multimap<int,string*> LumpMoleculesMap;//reverse of MolLumpMap;

	map<string*,int, classcomp>::iterator it;

	for (it=MolLumpMap.begin();it!=MolLumpMap.end();it++)
	{
		if ((*it->first)=="C1(C=C)(C(C=CC1C)C)C")
		{
			if (AdditionalLumpMap.count(it->second)>0)
				cout<<"the molecule underwent additional lump to get into AdditionalLumpMap[it->second]"<<endl;
			else
				cout<<"no additional lumping for this molecule"<<endl;
		}
		if (AdditionalLumpMap.count((*it).second)>0)
			LumpMoleculesMap.insert(pair<int,string*>(AdditionalLumpMap[(*it).second], (*it).first));
		else LumpMoleculesMap.insert(pair<int,string*>((*it).second, (*it).first));
	}

	ofstream myfile;
	myfile.open("Lumps.txt");

	multimap<int,string*>::iterator it2;
	int prevLump =-1;
	for (it2=LumpMoleculesMap.begin();it2!=LumpMoleculesMap.end();it2++)
	{
		
		if ((*it2).first!=prevLump)
		{	
			if(prevLump!=-1)myfile<<endl;
			myfile<<MolLumps[(*it2).first].getMoleculeString()<<" ("<<LumpMoleculesMap.count(it2->first)<<") "<<"--->  ";
		}
		prevLump =(*it2).first;
		if ((*it2).first==prevLump)myfile<<*((*it2).second)<<"   ";
	}
}

int rxn_net_gen::findSpeciesRank(std::string S)
{
	Molecule mol(S,moleculesize(S));
	mol.unique_smiles();
	if (AllMolecules.count(StringRegistry::getStringPointer(mol.moleculestring()))>0)
		return AllMolecules[StringRegistry::getStringPointer(mol.moleculestring())];
	else return -1;
}

bool rxn_net_gen::isSiteIntermediate(string s)
{

	//if the composite site pattern exists and it is not an initial reactant, return true 
	return containsSiteAtom(s) && (InitialReactants.count(StringRegistry::getStringPointer(s))==0);
}
	
bool rxn_net_gen::containsSiteAtom(string s)
{
	//if molecule has any composite atom that has been declared as a site
	if (CompositeSites.empty())return false;
	else
	{
		Molecule m(s, moleculesize(s));

		for (int i =0; i< CompositeSites.size();i++)
		{
			if (Patternmatch(m,Substructure(CompositeSites[i].first,patternsize(CompositeSites[i].first)),1).number_of_matches()>0)
				return true;
		}
		return false;
	}
}


int rxn_net_gen::getDeltaNGasPhase(generated_rxn& rxn)
{
	int deltaN = 0;

	for (int i =0;i<rxn.number_pdcts();i++)
		if (!containsSiteAtom(*rxn.get_products(i)))
			deltaN++;
	for (int i = 0;i<rxn.number_reactants();i++)
		if (!containsSiteAtom(*rxn.get_reactants(i)))
			deltaN--;
	
	return deltaN;
}


string rxn_net_gen::getRxnStringFromLumpedRxn(LumpedReaction L)
{
	
	//generated_rxn G;
	string rxnstring ="";
	vector<int> R = L.getReactantLumps();
	vector<int> P = L.getProductLumps();
	for (int i=0; i<R.size(); i++)
	{
		if (i!=0)rxnstring+=".";
		rxnstring+=MolLumps[R[i]].getMoleculeString();
	}
	rxnstring+=">>";
	for (int i=0; i<P.size(); i++)
	{
		if (i!=0)rxnstring+=".";
		rxnstring+=MolLumps[P[i]].getMoleculeString();
	}

	return rxnstring;
}

void rxn_net_gen::ReconstructReactions(int r1, int r2, ConstrPtr Constr1, ConstrPtr Constr2, CombinedConstrPtr combConstr, string rulename)
{
	RulesRemovedByReconstruction.insert(r1);
	RulesRemovedByReconstruction.insert(r2);
	if (LumpStrat.shoudLump())
		ReconstructForLumpedNetwork(r1,r2,Constr1,Constr2,combConstr, rulename);
	else
		ReconstructForOriginalNetwork(r1,r2,Constr1,Constr2, combConstr, rulename);
}

void rxn_net_gen::SetAllLumpedRxns()
{

	//We assume (and force) that this is called AFTER network generation and reconstruction (if any). 

	map<LumpedReaction, int, LumpedReactionCompare>::iterator it;

	for (it = LumpedReactionMap.begin();it!=LumpedReactionMap.end();it++)
	{
		if (RulesRemovedByReconstruction.count(it->first.getRule())==0)
		{
			generated_rxn r1 = GetRxnFromLumpedReaction(it->first);
			r1.set_frequency(LumpedReactionMap[it->first]);
			LumpedReconstructedNetwork.push_back(r1);
			SetAllLumpedMoleculesInfo(r1);

		}
	}

	//just making sure that this AllLumpedMolecules does not miss any of the initial reactants! 

	set<string*, classcomp>::iterator setIt;

	for (setIt = InitialReactants.begin();setIt!=InitialReactants.end();setIt++)
	{
		cout<<*(*setIt)<<" ";
		cout<<findLumpofMolecule(*setIt)<<endl;
		cout<<MolLumps[findLumpofMolecule(*setIt)].getRank()<<endl;
		AllLumpedMolecules.insert(pair<string*,int>(*setIt, MolLumps[findLumpofMolecule(*setIt)].getRank()));
	}

	
	for (int i =0;i<LumpedReconstructedNetwork.size();i++)
	{
		LumpedReconstructedReactionsMap.insert(pair<int,int>(LumpedReconstructedNetwork[i].get_rule(),i+1));
	}

	//Now, we can print THE final list of lumped and reconstructed network as well as final list of species lumps
	printFinalLumpSpeciesAndRxns();
	cout<<"The final set of lumps (after lumping and reconstruction, if any) is stored in LumpSpecies.txt"<<endl;
	cout<<"There are "<<AllLumpedMolecules.size()<<" molecules in this set."<<endl;
	cout<<"The final list of reactions (classified according to rules and after lumping and reconstruction, if any) is stored in LumpedReconstructedNetwork.txt"<<endl;
	cout<<"There are "<<LumpedReconstructedNetwork.size()<<" reactions is this list"<<endl;

}
								
void rxn_net_gen::SetAllLumpedMoleculesInfo(generated_rxn& r1)
{

	for (int i =0;i<r1.number_pdcts();i++)
	{
		MolLumpProductMap.insert(pair<string*,int> (r1.get_products(i),LumpedReconstructedNetwork.size()));
		AllLumpedMolecules.insert(pair<string*,int>(r1.get_products(i), MolLumps[findLumpofMolecule(r1.get_products(i))].getRank()));
	}
	for (int i = 0;i<r1.number_reactants();i++)
	{
		MolLumpReactantMap.insert(pair<string*,int> (r1.get_reactants(i),LumpedReconstructedNetwork.size()));
		AllLumpedMolecules.insert(pair<string*,int> (r1.get_reactants(i), MolLumps[findLumpofMolecule(r1.get_reactants(i))].getRank()));
	}
}



void rxn_net_gen::ReconstructForLumpedNetwork(int r1, int r2, ConstrPtr Constr1, ConstrPtr Constr2, CombinedConstrPtr combConstr, string rulename)
{
	map<LumpedReaction, int, LumpedReactionCompare>::iterator it;

	Reactiontype Rt;
	Rt.setRuleName(rulename);
	Rtlist.push_back(Rt);

	for (it = LumpedReactionMap.begin();it!=LumpedReactionMap.end();it++)
	{
		if (it->first.getRule()==r1)
		{
			vector<int> R = it->first.getReactantLumps();
			string molstring = MolLumps[R[0]].getMoleculeString();
			Molecule mol1((molstring),moleculesize(molstring));
			if ((Constr1)(mol1))
			{
				map<LumpedReaction, int, LumpedReactionCompare>::iterator it2;
				generated_rxn r1 = GetRxnFromLumpedReaction(it->first);
				r1.set_frequency(LumpedReactionMap[it->first]);
						
				
				for (it2 = LumpedReactionMap.begin();it2!=LumpedReactionMap.end();it2++)
				{

					if (it2->first.getRule()==r2)
					{
						vector<int> R2 = it2->first.getReactantLumps();
						string molstring2 = MolLumps[R2[1]].getMoleculeString();
						Molecule mol2((molstring2),moleculesize(molstring2));
						if ((Constr2)(mol2))
						{
							if ((combConstr)(mol1,mol2))
							{
								generated_rxn r2 = GetRxnFromLumpedReaction(it2->first);
								r2.set_frequency(LumpedReactionMap[it2->first]);
			

								generated_rxn r3 = GetOneRxnFromTwo(r1,r2);

								
								r3.set_rule(Rtlist.size()-1);

								LumpedReconstructedNetwork.push_back(r3);

								//for each reactant update MolLumpReactantMap
								// for each product update MolLumpProductMap
								SetAllLumpedMoleculesInfo(r3);
								
							}
						}
					}
				}
			}
		}
	}
}

generated_rxn rxn_net_gen::GetOneRxnFromTwo(generated_rxn& r1, generated_rxn& r2)
{
	map<string*,int> stoich;
	vector<pair<generated_rxn,int> > rxnPair;
	rxnPair.push_back(pair<generated_rxn,int>(r1,1));
	rxnPair.push_back(pair<generated_rxn,int>(r2,1));
					
	partialMechanism pM(&rxnPair);
					
	stoich = pM.getStoichiometry();
	map<string*,int>::iterator map_it;

	generated_rxn G;

	//Need to preserve hte order of reactants. First reactant of first reaction should be the first reactant in the final
	//reaction, and second reactant of the second reactant should be the second reactant of the reaction. No such constraints
	//on products

	for (int i =0;i<r1.number_reactants();i++)
	{
		if (stoich.count(r1.get_reactants(i))>0)
		{
			if (i==0)
				G.add_reactants(r1.get_reactants(i));
			else
			{
				if (stoich[r1.get_reactants(i)]>0)
				{
					for (int i=0;i<abs(stoich[r1.get_reactants(i)]);i++)
						G.add_reactants(r1.get_reactants(i));
				}
			}
		}
	}

	for (int i =0;i<r2.number_reactants();i++)
	{
		if (stoich.count(r2.get_reactants(i))>0)
		{
			if (i==1)
				G.add_reactants(r2.get_reactants(i));
			else
			{
				if (stoich[r2.get_reactants(i)]>0)
				{
					for (int i=0;i<abs(stoich[r2.get_reactants(i)]);i++)
						G.add_reactants(r2.get_reactants(i));
				}
			}
		}
	}


	for (map_it=stoich.begin();map_it!=stoich.end();map_it++)
	{
		if (map_it->second>0)
		{
			for (int i=0;i<abs(map_it->second);i++)
				G.add_products(map_it->first);							
		}
		
	
		if (map_it->second==0 && (*(map_it->first)== *rxnPair[0].first.get_reactants(0) || *(map_it->first)== *rxnPair[1].first.get_reactants(1) ))
			G.add_products(map_it->first);
	}
	G.set_frequency(rxnPair[0].first.get_frequency()*rxnPair[1].first.get_frequency());
	
	return G;
}
							





void rxn_net_gen::ReconstructForOriginalNetwork(int r1, int r2, ConstrPtr Constr1, ConstrPtr Constr2, CombinedConstrPtr combConstr, string rulename)
{
	cout<<"molecularity of r1 and r2 are "<<Rtlist[r1].get_molecularity()<<"  "<<Rtlist[r2].get_molecularity()<<endl;
	
	Reactiontype Rt;
	Rt.setRuleName(rulename);
	Rtlist.push_back(Rt);

	if (ReactionsMap.count(r1)>0 && ReactionsMap.count(r2)>0)
	{
		multimap<int, int>::iterator it1;
		multimap<int, int>::iterator it2;
		pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret1;
		pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret2;


		ret1 = ReactionsMap.equal_range(r1);
		ret2 = ReactionsMap.equal_range(r2);

		for (it1=ret1.first; it1!=ret1.second;it1++)
		{
			string molstring = *(AllReactions[it1->second-1].get_reactants(0));
			Molecule mol1((molstring),moleculesize(molstring));
			if ((Constr1)(mol1))
			{
				for (it2=ret2.first; it2!=ret2.second;it2++)
				{
					string molstring2 = *(AllReactions[it2->second-1].get_reactants(1));
					Molecule mol2((molstring2),moleculesize(molstring2));
					if ((Constr2)(mol2))
					{
						if ((combConstr)(mol1,mol2))
						{
							generated_rxn G = GetOneRxnFromTwo(AllReactions[it1->second-1],AllReactions[it2->second-1]);
							G.set_rule(Rtlist.size()-1);
							ReconstructedReactions.push_back(G);
						}
					}
				}
			}
		}
	}
}



void rxn_net_gen::generateStoichMatrix()
{
    ofstream file1, file2,file3, file4;
    file1.open("StoichRxns.txt");
    file2.open("StoichCoeff.txt");
	file3.open("Visualization.txt");
	file4.open("StoichMatrix.txt");

	

	vector<map<int,int> > StoichMatrixInfo;

	map<int,int> StoichMapForASpecie;

		
	StoichMatrixInfo.resize(AllMolecules.size(),StoichMapForASpecie);
	
	 
    pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator> pairIter1;
	file3<<"digraph G {"<<endl;
	
	map<string*,int,classcomp> StringIndexMap;
	

	for (int i=0; i<AllReactions.size();i++)
	{
		file3<<"  r"<<i<<"[shape=circle,color=black,style=filled,label=\"\",width=0.3];"<<endl;
	}
	

    string previousMol = "";

	map<string*,int,classcomp>::iterator it;
    int i = 0;
	for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		
		StringIndexMap.insert(pair<string*,int>(it->first,StringIndexMap.size()));
		file1<<*it->first;
		file2<<*it->first;
		file4<<*it->first;
		file3<<"  m"<<StringIndexMap[it->first];

		if (InitialReactants.count(it->first)>0)
			file3<<"[shape=hexagon,color=green,label=\"\",style=filled];"<<endl;
		else if (isSiteIntermediate(*it->first))
			file3<<"[shape=oval,color=red,label=\"\",style=filled];"<<endl;
		else file3<<"[shape=square,color=blue,label=\"\",style=filled];"<<endl;

		pairIter1=MolReactantMap.equal_range(it->first);

		multimap<string*,int>::iterator it2;

		for (it2 = pairIter1.first;it2!=pairIter1.second;it2++)
		{
			file1<<" "<<it2->second;
			int netStoich = AllReactions[it2->second-1].NetOccurence(it->first);
			if (netStoich < 0)
			{
				file2<<" "<<netStoich;
				StoichMatrixInfo[StringIndexMap.size()-1][it2->second-1]=netStoich;
			}
			file3<<"  m"<<StringIndexMap[it->first]<<"->r"<<(it2->second-1)<<" [label=\""<<AllReactions[it2->second-1].getReactantStoich(it->first)<<"\"];"<<endl;
		}
		
		pairIter1=MolProductMap.equal_range(it->first);

		for (it2 = pairIter1.first;it2!=pairIter1.second;it2++)
		{
			file1<<" "<<it2->second;
			int netStoich = AllReactions[it2->second-1].NetOccurence(it->first);
			if (netStoich > 0)
			{
				file2<<" "<<netStoich;
				StoichMatrixInfo[StringIndexMap.size()-1][it2->second-1]=netStoich;
			}
			file3<<"  r"<<(it2->second-1)<<"->m"<<StringIndexMap[it->first]<<" [label=\""<<AllReactions[it2->second-1].getProductStoich(it->first)<<"\"];"<<endl;
		}
		file1<<endl;
		file2<<endl;

		int rxnCounter = 0;

		map<int,int>::iterator map_it;
		

		for (map_it=StoichMatrixInfo[StringIndexMap.size()-1].begin();map_it!=StoichMatrixInfo[StringIndexMap.size()-1].end();map_it++)
		{
			while (rxnCounter<map_it->first)
			{
				file4<<" 0";
				rxnCounter++;
			}
			file4<<" "<<map_it->second;
			rxnCounter++;
			
		}
		while (rxnCounter<AllReactions.size())
		{
			file4<<" 0";
			
			rxnCounter++;
		}
		file4<<endl;
	}
	file3<<"}"<<endl;

		 
    file1.close();
    file2.close();
	file3.close();
	file4.close();

	cout<<"finished generating stoichiometric matrix"<<endl;

}

void rxn_net_gen::GenerateInfoForAthena()
{
	ofstream file,file2;

	file.open("StoichForAthena.txt");
	file2.open("RateExpressionForAthena.txt");

	
	for (int i=0;i<AllReactions.size();i++)
	{
		file<<"Stoich(:nc,"<<i+1<<")=(/";
		file2<<"r("<<i+1<<")=k("<<i+1<<")";
		
		

		map<string*,int,classcomp>::iterator it;
		
		for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
		{
			int IndexOfMol = AllReactions[i].NetOccurence(it->first);
			if (it!=AllMolecules.begin())
				file<<","<<IndexOfMol;
			else
				file<<IndexOfMol;
			for (int j=0;j<IndexOfMol;j++)
				file2<<"*c("<<distance(AllMolecules.begin(),it)+1<<")";

		}
		file<<"/)"<<endl;
		file2<<endl;

	}

	file.close();
	file2.close();

	cout<<"finished generating files for Athena"<<endl;

}


void rxn_net_gen::CalculateParametersFileForGAMS()
{
	ofstream parametersFile;
	
	parametersFile.open("ParametersForGAMS.txt");

	parametersFile<<"*------------------"<<endl;
	parametersFile<<"* Model sets"<<endl;
	parametersFile<<"*------------------"<<endl;
	parametersFile<<endl;
	parametersFile<<"Sets"<<endl;
	parametersFile<<endl;
	parametersFile<<"i	species in the reaction network /"<<endl;

	map<string*,int,classcomp>::iterator it;
		
	for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		parametersFile<<"\'"<<SMILESGAMSSpeciesMap[it->first]<<"\'"<<endl;
	}
	parametersFile<<"/"<<endl;
	parametersFile<<endl;
	parametersFile<<"i_r(i)	species that are reactants /"<<endl;

	set<string*,classcomp>::iterator reactants_it;

	for (reactants_it = InitialReactants.begin();reactants_it!=InitialReactants.end();reactants_it++)
		parametersFile<<"\'"<<SMILESGAMSSpeciesMap[*reactants_it]<<"\'"<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;
	parametersFile<<endl;

	parametersFile<<"j        reactions in the reaction network  /"<<endl;
	parametersFile<<"rxn1*rxn"<<AllReactions.size()<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;

	/*parametersFile<<"j_r      reactant inputs /"<<endl;
	parametersFile<<"rxnr1*rxnr"<<InitialReactants.size()<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;*/

	/*parametersFile<<"j_waste  waste and byproduct outputs /"<<endl;
	parametersFile<<"rxnw1*rxnw"<<AllMolecules.size()<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;*/

	parametersFile<<"j_type   reaction types /"<<endl;

	for (int i=0;i<Rtlist.size();i++)
		parametersFile<<"\'"<<Rtlist[i].getRuleName()<<"\'"<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;

	parametersFile<<"j_corr(j,j_type) corresponding reaction type for each reaction /"<<endl;

	for (int i=0;i<AllReactions.size();i++)
		parametersFile<<"rxn"<<i+1<<" . \'"<<Rtlist[AllReactions[i].get_rule()].getRuleName()<<"\'"<<endl;
	parametersFile<<"/"<<endl;

	parametersFile<<endl;
	parametersFile<<endl;
	parametersFile<<endl;

	parametersFile<<"*------------------"<<endl;
	parametersFile<<"* Model parameters"<<endl;
	parametersFile<<"*------------------"<<endl;

	parametersFile<<endl;

	parametersFile<<"Parameters"<<endl;

	parametersFile<<endl;

	/*parametersFile<<"A_1(i,j_r)       reactant inputs stoichiometric coefficient matrix /"<<endl;

	for (reactants_it = InitialReactants.begin();reactants_it!=InitialReactants.end();reactants_it++)
		parametersFile<<"\'"<<*(*reactants_it)<<"\' . rxnr"<<distance(InitialReactants.begin(),reactants_it)+1<<" 1"<<endl;
	parametersFile<<"/"<<endl;
	
	parametersFile<<endl;*/

	parametersFile<<"A_2(i,j)         stoichiometric coefficient matrix /"<<endl;

	pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator> pairIter1;
	
	for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		pairIter1=MolReactantMap.equal_range(it->first);

		multimap<string*,int>::iterator it2;
		int netStoich;
		for (it2 = pairIter1.first;it2!=pairIter1.second;it2++)
		{
			if (DenticityInfoMap.count(*it->first)>0)
			{
				netStoich = - CalculateNetStoichDifferenceWithDenticity(AllReactions[it2->second-1],*it->first);
				//if (it2->second==4)
				//	cout<<AllReactions[it2->second-1].reactionstring()<<"  "<<netStoich<<endl;
			}
			else
				netStoich = AllReactions[it2->second-1].NetOccurence(it->first);

			if (netStoich < 0)
				parametersFile<<"\'"<<SMILESGAMSSpeciesMap[it->first]<<"\' . rxn"<<it2->second<<" "<<netStoich<<endl;	
		}
	
		pairIter1=MolProductMap.equal_range(it->first);

		for (it2 = pairIter1.first;it2!=pairIter1.second;it2++)
		{
			if (DenticityInfoMap.count(*it->first)>0)
				netStoich = - CalculateNetStoichDifferenceWithDenticity(AllReactions[it2->second-1],*it->first);
			else
				netStoich = AllReactions[it2->second-1].NetOccurence(it->first);

			if (netStoich > 0)
				parametersFile<<"\'"<<SMILESGAMSSpeciesMap[it->first]<<"\' . rxn"<<it2->second<<" "<<netStoich<<endl;	
		}
	}
	parametersFile<<"/"<<endl;

	parametersFile<<endl;

/*	parametersFile<<"A_3(i,j_waste)   waste and byproduct stoichiometric coefficient matrix /"<<endl;

	for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		parametersFile<<"\'"<<(*it->first)<<"\' . rxnw"<<distance(AllMolecules.begin(),it)+1<<" 1"<<endl;
	}
	parametersFile<<"/"<<endl;*/

	parametersFile.close();
}


void rxn_net_gen::CheckDensityInfo(multimap<string,string>& denticity)
{
	DenticityInfoMap = denticity;
}

int rxn_net_gen::CalculateNetStoichDifferenceWithDenticity(generated_rxn& rxn, string freeSite)
{
	//string denticity = "$[|~" + pattern + "|>=1]";
	//string denticity2 = pattern + "[|-H|>=1]";
	//string denticity3 ="C[|-" + pattern + "|>=1]C[|-H|==1]=O";

	int netStoich = 0;

	pair<multimap<string,string>::iterator, multimap<string,string>::iterator > pairIter;

	pairIter = DenticityInfoMap.equal_range(freeSite);

	for (multimap<string,string>::iterator it = pairIter.first;it!=pairIter.second;it++)
	{
		netStoich += rxn.NetPatternDiff(it->second);
	}

	return netStoich;
}


double rxn_net_gen::calculateThermoOfRxn(ThermoType type, generated_rxn& rxn, double Temp)
{
	double ThermoRxn = 0.0;
	for (int j=0;j<rxn.number_pdcts();j++)
	{
		if (AllMolThermo.count(rxn.get_products(j))>0 && AllMolThermo[rxn.get_products(j)].isTempAvailable(type,Temp)>=0)
			ThermoRxn+=AllMolThermo[rxn.get_products(j)].getThermo(type, Temp);
		else
		{
			double Thermo = 0.0;
			calcMolThermo(rxn.get_products(j),type,Temp,Thermo);
			ThermoRxn+=Thermo;
			AllMolThermo[rxn.get_products(j)].setThermo(type,Temp,Thermo);
			
		}
	}
	for (int j=0;j<rxn.number_reactants();j++)
	{

		if (AllMolThermo.count(rxn.get_reactants(j))>0 && AllMolThermo[rxn.get_reactants(j)].isTempAvailable(type,Temp)>=0)
			ThermoRxn-=AllMolThermo[rxn.get_reactants(j)].getThermo(type, Temp);
		else
		{
			double Thermo = 0.0;
			calcMolThermo(rxn.get_reactants(j),type,Temp,Thermo);
			ThermoRxn-=Thermo;
			AllMolThermo[rxn.get_reactants(j)].setThermo(type,Temp,Thermo);
			//cout<<"calculating thermo"<<endl;
		}
	}

	return ThermoRxn;
}

bool rxn_net_gen::calcMolThermo(string* molptr,ThermoType type, double Temp, double& value)
{
	Molecule mol(*molptr, moleculesize(*molptr));
	mol.unique_smiles();

	if (type==EnthalpyType)
		ThermoGA::calculateDeltaH(mol,containsSiteAtom(*molptr), Temp, value);
	else if (type ==EntropyType)
		ThermoGA::calculateDeltaS(mol,containsSiteAtom(*molptr), Temp, value);
	else if (type ==CpType)
		ThermoGA::calculateCp(mol,containsSiteAtom(*molptr), Temp, value);
	else if (type ==FreeEnergyType)
		ThermoGA::calculateDeltaG(mol,containsSiteAtom(*molptr), Temp, value);
	else {
		cout<<"wrong type"<<endl;
		return false;
	}
	return true;
}

int rxn_net_gen::FindRxnWithReactantsAndRule(const set<string>& reactants, int rule)
{
	multimap<int,int>::iterator mmap_it;
	pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;
	ret = ReactionsMap.equal_range(rule);

	for (mmap_it= ret.first;mmap_it!=ret.second;mmap_it++)
	{
		set<string> reactantsset;

		for (int j=0;j<AllReactions[mmap_it->second-1].number_reactants();j++)
			reactantsset.insert(*AllReactions[mmap_it->second-1].get_reactants(j));
		if (reactants == reactantsset)
			return mmap_it->second;
	}
	return -1;
}

set<int> rxn_net_gen::FindReactionsWithReactantsAndRule(const set<string>& reactants, int rule)
{
	set<int> setOfRxns;
	multimap<int,int>::iterator mmap_it;
	pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;
	ret = ReactionsMap.equal_range(rule);

	for (mmap_it= ret.first;mmap_it!=ret.second;mmap_it++)
	{
		set<string> reactantsset;

		for (int j=0;j<AllReactions[mmap_it->second-1].number_reactants();j++)
			reactantsset.insert(*AllReactions[mmap_it->second-1].get_reactants(j));
		if (reactants == reactantsset)
			setOfRxns.insert(mmap_it->second);
	}
	return setOfRxns;
}

void rxn_net_gen::findRxnsOfSimultRules(vector<pair<int,int> > & simultaneousrules)
{
	//for each set of reactants
		//for each pair, take teh first element - which indicates a rule
		  // find the first reaction consuming the reactants in the set from amongst reactions with this rule
		  // take the second element 
			 // if bimolecular
			 //again find the first reaction of this type that consumes the reactants
		    // if unimolecular
		     // find the reactions that consume of the reactants.
	//if both rules are unimolecular, find reactions of the first rule, and check if there are reactions of the second type that
	//have the same set of reactions. 

	set<set<string> >::iterator it;

	for (it = reactantpairsForBiMol.begin(); it!=reactantpairsForBiMol.end();it++)
	{
		//go through the simultaneous rules
		for (int i=0;i<simultaneousrules.size();i++)
		{
			int firstRxn = FindRxnWithReactantsAndRule(*it, simultaneousrules[i].first);
			if (firstRxn>=0)
			{
				if (Rtlist[simultaneousrules[i].second].get_molecularity()>1)
				{		
					int secondRxn = FindRxnWithReactantsAndRule(*it,simultaneousrules[i].second);
					if (secondRxn>=0)
					{
						set<int> simultRxns;
						simultRxns.insert(firstRxn);
						simultRxns.insert(secondRxn);
						SimultaneousRxns.push_back(simultRxns);
					}
				}
				else
				{
					set<string>::iterator mset_it;
					for (mset_it = it->begin();mset_it!=it->end();mset_it++)
					{
						set<string> OneRxn;
						OneRxn.insert(*mset_it);
						int secondRxn  = FindRxnWithReactantsAndRule(OneRxn,simultaneousrules[i].second);
						if (secondRxn>=0)
						{
							set<int> simultRxns;
							simultRxns.insert(firstRxn);
							simultRxns.insert(secondRxn);
							SimultaneousRxns.push_back(simultRxns);
							
						}
					}
				}
			}
		}
	}

	//this is for the case where both rules are unimolecular
	for (int i=0;i<simultaneousrules.size();i++)
	{
		if (Rtlist[simultaneousrules[i].first].get_molecularity()==1 && Rtlist[simultaneousrules[i].second].get_molecularity()==1)
		{
			multimap<int,int>::iterator rule_it;

			pair<multimap<int,int>::iterator, multimap<int,int>::iterator > pair_it;

			pair_it = ReactionsMap.equal_range(simultaneousrules[i].first);

			for (rule_it=pair_it.first;rule_it!=pair_it.second;rule_it++)
			{
				set<string> reactantsSet;
				reactantsSet.insert(*AllReactions[rule_it->second-1].get_reactants(0));//only one reactant - unimolecular!
				
				//find a reaction belonging to the second rule
				int secondRxn = FindRxnWithReactantsAndRule(reactantsSet,simultaneousrules[i].second);

				//these two reactions are simultaneous

				if (secondRxn>=0)
				{
					set<int> simultRxns;
					simultRxns.insert(rule_it->second);//note, it's not rule_it->second-1
					simultRxns.insert(secondRxn);
					SimultaneousRxns.push_back(simultRxns);
				}
			}
		}
	}




}

void rxn_net_gen::GatherSimultaneousAndJointReactionsForGAMS(vector<pair<int,int> >& simultaneousrules, vector<pair<int,int> >& ccformationsteps)
{
	IdentifySimultaneousBimolRxns();//needed only for bimolecular rxns because simultaneous unimolecular have already been found during network generation
	findRxnsOfSimultRules(simultaneousrules);
	findTransformRxns(ccformationsteps);
	
	ofstream jointsfile("j_joint.txt");
	ofstream jtransformfile("j_transform.txt");

	
	jointsfile<<"sets    joints /1*"<<SimultaneousRxns.size()<<"/"<<endl;
	jointsfile<<endl;
	jointsfile<<"        j_joint(joints,j)   reactions that occur simultaneously /"<<endl;

	for (int i=0;i<SimultaneousRxns.size();i++)
	{
		jointsfile<<i+1<<". (";

		set<int>::iterator it;

		for (it=SimultaneousRxns[i].begin();it!=SimultaneousRxns[i].end();it++)
		{
			if (it!=SimultaneousRxns[i].begin())jointsfile<<", ";
			jointsfile<<"rxn"<<*it;

		}
		jointsfile<<")"<<endl;
	}
	jointsfile<<"/"<<endl;
	jointsfile.close();

	jtransformfile<<"set      j_transforms(j,jj)    reactions that cannot occur without related reactions /"<<endl;
	map<int, vector<int> >::iterator it;

	for (it = RxnsToBeSummedUpMap.begin();it!=RxnsToBeSummedUpMap.end();it++)
	{
		jtransformfile<<"rxn"<<it->first<<" . (rxn"<<it->first;
		for (int i=0;i<it->second.size();i++)
		{
			jtransformfile<<",rxn"<<it->second.at(i);
		}
		jtransformfile<<")"<<endl;
	}
	jtransformfile<<"/"<<endl;
	jtransformfile.close();
}

void rxn_net_gen::findTransformRxns(std::vector<pair<int,int> > &TransformSteps)
{
	set<set<string> >::iterator it;

	for (it = reactantpairsForBiMol.begin(); it!=reactantpairsForBiMol.end();it++)
	{
		//go through the TransformSteps rules
		for (int i=0;i<TransformSteps.size();i++)
		{
			int firstRxn = FindRxnWithReactantsAndRule(*it, TransformSteps[i].first);
			if (firstRxn>=0)
			{
				set<string>::iterator mset_it;
				for (mset_it = it->begin();mset_it!=it->end();mset_it++)
				{
					set<string> OneRxn;
					OneRxn.insert(*mset_it);
					int secondRxn  = FindRxnWithReactantsAndRule(OneRxn,TransformSteps[i].second);
					if (secondRxn>=0)
					{
						RxnsToBeSummedUpMap[firstRxn];//this will insert a new value if firstRxn is not alreay there
						RxnsToBeSummedUpMap[firstRxn].push_back(secondRxn);
					}
				}
			}
		}
	}
}

void rxn_net_gen::IdentifySimultaneousBimolRxns()
{

	set<set<string> >::iterator it;

	for (it = reactantpairsForBiMol.begin(); it!=reactantpairsForBiMol.end();it++) //reactantpairsForBiMol is populated during network generation
	{
		for (int i = 0;i<Rtlist.size();i++)
		{
			if (Rtlist[i].get_molecularity()>1)// if bimolecular
			{
				set<int> rxns = FindReactionsWithReactantsAndRule(*it,i);
				if (rxns.size()>1)SimultaneousRxns.push_back(rxns);
			}
		}
	}
}


struct IsLowerChar{
	bool operator()(int value)
	{
		return ::islower((unsigned char)value);
	}
};

void rxn_net_gen::PrepareSpeciesInfoForGAMS()
{
	map<string*,int,classcomp>::iterator it;

	for (it=AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		string GAMSSpeciesString = *it->first;
		string::iterator mol_it = find_if(it->first->begin(), it->first->end(), IsLowerChar());//there is lower case char
		if (mol_it!=it->first->end())	
			GAMSSpeciesString+="A";
		SMILESGAMSSpeciesMap[it->first]=GAMSSpeciesString;
	}
}

void rxn_net_gen::SetTemp(double temp)
{
	Temperature = temp;

}

void rxn_net_gen::printFinalLumpSpeciesAndRxns()
{
	ofstream LumpSpeciesFile("LumpSpecies.txt");
	ofstream LumpedReconstructedNetworkFile("LumpedReconstructedNetwork.txt");

	LumpSpeciesFile<<"This file contains the final list of lumps"<<endl;
	LumpedReconstructedNetworkFile<<"This file consists of the final list of lumped (and reconstructed) network"<<endl;

	
	multimap<int,int>::iterator it;

	int prevRule = -1;
	for (it=LumpedReconstructedReactionsMap.begin();it!=LumpedReconstructedReactionsMap.end();it++)
	{
		
		if (it->first!=prevRule)
		{
			LumpedReconstructedNetworkFile<<"Reactions of rule ";
			LumpedReconstructedNetworkFile<<Rtlist[it->first].getRuleName()<<":"<<endl;
			prevRule = it->first;
		}
		LumpedReconstructedNetworkFile<<LumpedReconstructedNetwork[(*it).second-1].reactionstring()<<"  "<<LumpedReconstructedNetwork[(*it).second-1].get_frequency()<<endl;
		

	}

	map<string*,int, classcomp>::iterator lumpedIt;
	
	for (lumpedIt = AllLumpedMolecules.begin();lumpedIt!=AllLumpedMolecules.end();lumpedIt++)
		LumpSpeciesFile<<*(lumpedIt->first)<<endl;

	LumpSpeciesFile.close();
	LumpedReconstructedNetworkFile.close();

}

void rxn_net_gen::GetNetworkFromFile(const char * reactionsfile, const char * speciesfile)
{

	//First, set initial reactants in the normal way. 
	//Get all species from the species file. 
	//For each species, determine if it is a proper SMILES-like string
	//Construct the molecule, check for valency errors and throw any errors
	// Check global constraints
	// Add into AllMolecules. Initialize all ranks to 0?
	
	//Move into reactions file 
	//Check for the keyword "Rule" and pick the rule name - figure out the rule index.
	//Start reading reaction strings, figure out hte reaction SMILES and product SMILES. 
	//Create the "generated_rxn" object.
	//Populate AllReactions, ReactionsMap, MolReactantsMap, and MolProductsMap -- Anything else?

	//Set Molecule ranks by traversing the network
	//Lump molecule if required. 
	

	if (!setInitialReactants()) throw 1;
	unprocessedmol.clear(); //setInitialReactants populates unprocessedmol -- need to clear it (to remove Molecule instance from the heap)
	
	

	
	map<int, string*> SpeciesByIndex; 
	cout<<"Proceeding to read files"<<endl;
	
	string SpeciesListFile(speciesfile);
	cout<<"reading species from "<<SpeciesListFile<<endl;
	if (!ReadSpeciesFile(speciesfile, SpeciesByIndex))
	{
		cout<<"error reading species file"<<endl;
		throw 1;
	}

	string RxnNetwork(reactionsfile);
	cout<<"reading reactions from "<<RxnNetwork<<endl;
	if (!ReadReactionFile(reactionsfile, SpeciesByIndex))
	{
		cout<<"error reading reactions file"<<endl;
		throw 1;
	}

	TraverseNetworkToGetRanks();
	if (LumpStrat.shoudLump())LumpAllMolecules();

	DoPONALumping();
	DoMoreLumping();
	if (LumpStrat.shoudLump()) LumpReactions();
	PrintLumps();
	printOutputInfo();
	if (shouldCalcThermo) calculateThermoValues();

	cout<<"Network generation from files completed"<<endl;
}

bool rxn_net_gen::ReadSpeciesFile(const char* speciesfile, map<int,string*>& SpeciesByIndex)
{
	ifstream speciesFile(speciesfile);
	

	string line = "";
	int linenumber = 0; 
	bool isFileGood = true;
	
	if (speciesFile.is_open())
	{
		while (speciesFile.good())
		{
			getline(speciesFile,line);
			vector<string> tokens;
			linenumber++;

			if (TokenizeIntoWords(line, 0, tokens))
			{
				if (tokens.size()==1)
				{
					if (tokens[0] !="END")
					{
						throwsmileserror(tokens[0]);

						Molecule* mol = new Molecule(tokens[0], moleculesize(tokens[0]));
						mol->unique_smiles();
						pair<string, int> error = mol->checkValencyError();

						if (error.second ==0)
						{
							if (globalconstraintcheck(*mol))
							{
								string * molptr = StringRegistry::getStringPointer(mol->moleculestring());
								AllMolecules.insert(pair<string*, int>(molptr, -1));
								SpeciesByIndex.insert(pair<int,string*>(SpeciesByIndex.size(), molptr));
								if (mol->isIntermediate() || isSiteIntermediate(*molptr)) Intermediates.insert(molptr);
								processedmol.push_back(mol);
						
							}
							else
							{
								cout<<"error in line: "<<linenumber<<": Molecule does not satisfy global constraints!"<<endl;
								isFileGood = false;
								break;
							}
								
						}
						else
						{
							cout<<"error in line: "<<linenumber<<": ";
							MolValencyErrorStatements(error);
							isFileGood = false;
							break;
						}
					}
					else
					{
						cout<<"End of file reached!"<<endl;
					}
				}
				else{
					cout<<"error in line: "<<linenumber<<": One (and only one) word per line"<<endl;
					isFileGood = false;
					break;
				}
			}
		}
	}
	else
	{
		cout<<"file not open!"<<endl;
		isFileGood = false;
	}

	return isFileGood;

	
}

bool rxn_net_gen::ReadReactionFile(const char* reactionsfile, map<int, string*>& SpeciesByIndex)
{

	ifstream rxnfile(reactionsfile);
	

	map<string, int> RuleAndIndexMap; //stores the rule index of each rule known by its rule name;

	for (int i =0;i<Rtlist.size();i++)
	{
		RuleAndIndexMap[Rtlist[i].getRuleName()] = i;
	}

	string line = "";

	int RuleIndex = -1;
	int linenumber = 0;

	bool isGoodFile = true;

	if (rxnfile.is_open())
	{
		while (rxnfile.good())
		{
			getline(rxnfile,line);
			linenumber++;
			vector<string> tokens;

			if(TokenizeIntoWords(line,0,tokens))
			{
				if (tokens.size()>0) 
				{
					if (tokens[0]=="Rule")
					{
						if (tokens.size()==2)
						{
							if (RuleAndIndexMap.count(tokens[1])>0)
								RuleIndex = RuleAndIndexMap[tokens[1]];
							else
							{
								cout<<"error in line "<<linenumber<<": Rule not specified"<<endl;
								isGoodFile = false;
								break;
							}
						}
						else
						{
							cout<<"error in line "<<linenumber<<":There can only be two words in this line"<<endl;
							isGoodFile = false;
							break;
						}
					}
					else if (tokens[0]=="END")
					{
						cout<<"end of file reached!"<<endl;
					}
					else
					{
						//checking if string has a reaction delimiter -- a quick way to check if we do indeed have a reaction

						if (!InterpretRxnsAndUpdate(tokens,linenumber,SpeciesByIndex,RuleIndex))
						{
							isGoodFile = false;
							break;
						}
					}
				}
				else
				{
					cout<<"error in line "<<linenumber<<": This line is empty!"<<endl;
					isGoodFile =false;
					break;
				}
			}
		}
	}
	else
	{
		cout<<"file not open!"<<endl;
		isGoodFile = false;
	}

	return isGoodFile;
	
}

bool rxn_net_gen::InterpretRxnsAndUpdate(vector<string>& tokens, int linenumber, map<int,string*>& SpeciesByIndex, int RuleIndex)
{
	bool isRxnDelimiterPresent = false;
	bool isRxnFrequencyDelimiterPresent = false;
	for (int i =0;i<tokens.size();i++)
	{
		if (tokens[i]==">>"){
			isRxnDelimiterPresent = true;
		}
		if (tokens[i]=="|") { isRxnFrequencyDelimiterPresent =true;
			break;
		}

	}

	bool isReactant = true;
	bool isValidRxn = true;
	bool isFrequency = false;
	generated_rxn G;

	if (!isRxnFrequencyDelimiterPresent)
	{
		cout<<"no frequency given in line "<<linenumber<<"!"<<endl;
		return false;
	}

	if (isRxnDelimiterPresent)
	{
		
		for (int i =0;i<tokens.size();i++)
		{
			if (tokens[i]==">>")
				isReactant = false;
			else if (tokens[i]=="|")
				isFrequency = true;
			else
			{
				if (isFrequency)
				{
					int freq = -1;
					try{freq = StringToInt(tokens[i]);}
					catch (int error)
					{ if (error ==1) cout<<"error in line "<<linenumber<<": cannot convert "<<tokens[i]<<"into integer!"<<endl;}

					if (freq <=0)
					{
						cout<<"Not a valid frequency, hence not a valid reaction!"<<endl;
						isValidRxn = false;
						break;
					}
					else G.set_frequency(freq);


				}
				else
				{
					int species = -1;
					bool CorrectSpecies = true;
					try{species = StringToInt(tokens[i]);}
					catch (int error)
					{ if (error ==1) cout<<"error in line "<<linenumber<<": cannot convert "<<tokens[i]<<"into integer!"<<endl; CorrectSpecies = false;}

					if (SpeciesByIndex.count(species)==0)
					{
						cout<<"error in line "<<linenumber<<": cannot locate species with index "<<species<<endl;
						CorrectSpecies = false;
					}
					
					if (CorrectSpecies)
					{
						if (isReactant)
							G.add_reactants(SpeciesByIndex[species]);
						else
							G.add_products(SpeciesByIndex[species]);
					}
					else
					{
						cout<<"error in line "<<linenumber<<": Not a valid species, hence not a valid reaction!"<<endl;
						isValidRxn = false;
						break;
					}
				}
			}
		}
	}
	else
	{
		cout<<"error in line "<<linenumber<<": There is no reaction delimiter >> in the line"<<endl;
		isValidRxn = false;
	}

		
	if (isValidRxn)
	{

		//first check if it is net zero - if it is, remove!
		vector<pair<generated_rxn, int> > Rxns;

		Rxns.push_back(pair<generated_rxn, int> (G, 1));

		partialMechanism P(&Rxns);

		ofstream wrongrxn;

		if(linenumber==0)
			wrongrxn.open("wrongRxns.txt");
		else
			wrongrxn.open("wrongRxns.txt",ios::app);


		if (P.isNetZero())
		{
			cout<<"Reaction on line "<<linenumber<<" (rule: "<<Rtlist[RuleIndex].getRuleName()<<") is a net zero reaction! Omitting this reaction!"<<endl;
		}
		else if (G.netMassDiff()!=0)
		{
			cout<<"Reaction on line "<<linenumber<<" (rule: "<<Rtlist[RuleIndex].getRuleName()<<") does not conserve mass! Omitting this reaction!"<<endl;
			wrongrxn<<G.reactionstring()<<"  "<<linenumber<<endl;
		}
		else
		{

			G.set_rule(RuleIndex);
			AllReactions.push_back(G);
			ReactionsMap.insert(pair<int,int> (RuleIndex, AllReactions.size()));
			UpdateMolAsReactantsProductsMap(G);
		}
	}



	return isValidRxn;
}


bool rxn_net_gen::TraverseNetworkToGetRanks()
{
	//Do a breadth-first search of the network.

	//Start from initial reactants and successively proceed to cover all species;

	multimap<int, string*> RankSpeciesMap; //keeps a sorted list of all species in terms of its ranks.

	int SpeciesAssignedRanks = 0; 

	set<string*, classcomp>::iterator it;

	
	for (it = InitialReactants.begin();it!=InitialReactants.end();it++)
	{
		RankSpeciesMap.insert(pair<int, string*> (0,*it));
		SpeciesAssignedRanks++;
		//cout<<"Adding initial reactants"<<endl;
	}

	multimap<int,string*>::iterator mIt;
	mIt = RankSpeciesMap.begin();
	int NumberSpeciesChecked = 0;

	bool continueTraversing = true;

	while (continueTraversing)
	{
		//since this is sorted according to the ranks, iterating through this map is as good as doing a BFS!

		//go over all reactions wherein this species (mIt) is a reactant

		multimap<string*,int>::iterator it;

		pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator > pairIt;

		//cout<<*mIt->second<<" "<<MolReactantMap.count(mIt->second)<<endl;
		//cout<<NumberSpeciesChecked<<" "<<AllMolecules.size()<<" "<<SpeciesAssignedRanks<<" "<<RankSpeciesMap.size()<<endl;
		pairIt = MolReactantMap.equal_range(mIt->second);
		//cout<<"found equal range"<<endl;

		for (it=pairIt.first;it!=pairIt.second;it++)
		{
			//get the reaction and iterate through its products and assign ranks if not already assigned

			generated_rxn G = AllReactions[it->second-1];

			for (int i =0;i<G.number_pdcts();i++)
			{
				string* product = G.get_products(i);
				int pdctRank = AllMolecules[G.get_products(i)];
				if (pdctRank==-1 || pdctRank > mIt->first+1)
				{
					AllMolecules[product]=mIt->first+1;
					SpeciesAssignedRanks++;
					RankSpeciesMap.insert(pair<int, string*> (mIt->first+1, product));
					
					//TODO: Unfortunately, I cannot update information on parents for reactions because I dont know
					//how many atoms came from each reactant. Need to get this from input file, perhaps! 
				
				}
			}
				
		}
		NumberSpeciesChecked++;

		if (NumberSpeciesChecked>=RankSpeciesMap.size())
			continueTraversing =false;
		else
			mIt++;
	}

	map<string*,int, classcomp>::iterator molIt;

	for (molIt = AllMolecules.begin();molIt!=AllMolecules.end();molIt++)
	{
		int molRank = molIt->second;
		if (RankSpeciesMap.count(molRank)==0)
			cout<<"species "<<*(molIt->first)<<" is not reachable in the network traversal!"<<endl;
		int rxnsAsReactant, rxnsAsProduct;
		rxnsAsReactant = MolReactantMap.count(molIt->first);
		rxnsAsProduct = MolProductMap.count(molIt->first);

		if (rxnsAsReactant > 0 && rxnsAsProduct ==0)
			cout<<"species "<<*(molIt->first)<<" is never a product!"<<endl;
		if (rxnsAsReactant == 0 && rxnsAsProduct >0)
			cout<<"species "<<*(molIt->first)<<" is never a reactant!"<<endl;
	}

	


	return true;
}


void rxn_net_gen::LumpAllMolecules()
{
	map<string*,int,classcomp>::iterator it;

	for (it = AllMolecules.begin();it!=AllMolecules.end();it++)
	{
		if (InitialReactants.count(it->first)==0)
		{
			Molecule mol (*it->first, moleculesize(*it->first));
			mol.unique_smiles();
			LumpMolecule(mol,it->first);
		}
	}
}			
				

void rxn_net_gen::StoreRxnsAndSpecies(const char* speciesfileName, const char* rxnfileName, bool ForLumpedNetwork)
{
	ofstream speciesfile(speciesfileName);
	ofstream rxnfile(rxnfileName);

	map<string*, int> SpeciesIndex; 

	map<string*,int,classcomp>::iterator Mapit;

	map<string*,int,classcomp>* Molecules;
	vector<generated_rxn>* Reactions;
	multimap<int,int>* ReactionsMapPtr;

	if (ForLumpedNetwork)
	{
		Molecules = &AllLumpedMolecules;
		Reactions = &LumpedReconstructedNetwork;
		ReactionsMapPtr = &LumpedReconstructedReactionsMap;
	}
	else
	{
		Molecules = &AllMolecules;
		Reactions = &AllReactions;
		ReactionsMapPtr = &ReactionsMap;
	}


	for (Mapit = Molecules->begin();Mapit!=Molecules->end();Mapit++)
	{
		speciesfile<<*(Mapit->first)<<endl;

		SpeciesIndex[Mapit->first] = SpeciesIndex.size();
	}
	speciesfile<<"END";

	int prevRule = -1;
	multimap<int,int>::iterator it;
	for (it=ReactionsMapPtr->begin();it!=ReactionsMapPtr->end();it++)
	{
		
		if (it->first!=prevRule)
		{
			rxnfile<<"Rule ";
			rxnfile<<Rtlist[it->first].getRuleName()<<endl;
			prevRule = it->first;
		}

		if (RulesRemovedByReconstruction.count(it->first)==0)
		{

			for (int i =0; i<Reactions->at((*it).second-1).number_reactants();i++)
				rxnfile<<SpeciesIndex[Reactions->at((*it).second-1).get_reactants(i)]<<" ";
		
			rxnfile<<">>";

			for (int i =0; i<Reactions->at((*it).second-1).number_pdcts();i++)
				rxnfile<<" "<<SpeciesIndex[Reactions->at((*it).second-1).get_products(i)];

			rxnfile<<" | "<<Reactions->at((*it).second-1).get_frequency()<<endl;
		}

	}

	if (!ForLumpedNetwork && RulesRemovedByReconstruction.size()>0)
	{
		for (int i = 0; i<ReconstructedReactions.size(); i++)
		{
			if (ReconstructedReactions[i].get_rule()!=prevRule)
			{
				rxnfile<<"Rule ";
				rxnfile<<Rtlist[ReconstructedReactions[i].get_rule()].getRuleName()<<endl;
				prevRule = ReconstructedReactions[i].get_rule();
			}
			for (int j =0; j<ReconstructedReactions[i].number_reactants();j++)
				rxnfile<<SpeciesIndex[ReconstructedReactions[i].get_reactants(j)]<<" ";
		
			rxnfile<<">>";

			for (int j =0; j<ReconstructedReactions[i].number_pdcts();j++)
				rxnfile<<" "<<SpeciesIndex[ReconstructedReactions[i].get_products(j)];

			rxnfile<<" | "<<ReconstructedReactions[i].get_frequency()<<endl;
		}
	}


	rxnfile<<"END";

	speciesfile.close();
	rxnfile.close();

}

void rxn_net_gen::setKinetics(vector<KineticParamPtr>& kinFns)
{
	kinetics = &kinFns;
}

int rxn_net_gen::firstGasSpecies(int rxn, int rOrp)
{
	if (rOrp<0 || rOrp >1)
		return -1;
	else
	{
		if (rOrp ==0)
		{
			for (int i =0;i<AllReactions[rxn].number_reactants();i++)
			{
				if (!containsSiteAtom(*AllReactions[rxn].get_reactants(i)))
					return i;
			}
		}
		else
		{
			for (int i =0;i<AllReactions[rxn].number_pdcts();i++)
			{
				if (!containsSiteAtom(*AllReactions[rxn].get_products(i)))
					return i;
			}
		}

	}
	return -1;
}

				

	

	








			

       
