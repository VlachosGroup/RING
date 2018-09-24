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
#include "additionalfunc.h"
#include "stringreg.h"


Mechanisms::Mechanisms(rxn_net_gen* inputNetwork, const char* fname)
{
	Network = inputNetwork;
	filename = fname;
	DistinctDirectMechs = false;
	canCalculateActE = false;

}

bool Mechanisms::isComplete(partialMechanism* pMech)
{
	bool value = true;
	map<string*,int> MechStoich = pMech->getStoichiometry();

	map<string*,int>::iterator it;
	if (MechStoich.size()>0)
	{
		for (it=MechStoich.begin();it!=MechStoich.end();it++)
		{
			if (Network->Intermediates.count(it->first)>0)
			{
				if (it->second!=0)
				{
					value = false;
					break;
				}
			}
		}
	}
	else value = false;
	return value;
}

map<string*,int> Mechanisms::getIntermediates(partialMechanism* pMech)
{
	map<string*,int> Intermediates;
	map<string*,int>::iterator it;
	for (it=pMech->stoichiometry.begin();it!=pMech->stoichiometry.end();it++)
	{
		if (Network->Intermediates.count(it->first)>0 && it->second!=0)
			Intermediates[it->first] = it->second;
	}
	return Intermediates;
	
}

bool Mechanisms::ContainsInitialReactantsAsProducts(partialMechanism* pMech)
{
	bool contains = false;
	set<string*, classcomp>::iterator it;

	for (it=Network->InitialReactants.begin();it!=Network->InitialReactants.end();it++)
	{
		if (pMech->getFreq(*it)>0)
		{
			contains = true;
			break;
		}
	}
	return contains;
}

map<int,int> Mechanisms::getMechanismMap(vector <pair<int, int> >& ReactionIndices)
{
	map<int,int> mechanismMap;
	int MinValue = 10000000000;
			
	for (int i=0;i<ReactionIndices.size();i++)
	{
		if (mechanismMap.count(ReactionIndices[i].first)>0)
			mechanismMap[ReactionIndices[i].first]+=ReactionIndices[i].second;
		else 
			mechanismMap[ReactionIndices[i].first]=ReactionIndices[i].second;
	}

	for (map<int,int>::iterator it = mechanismMap.begin();it!=mechanismMap.end();it++)
		if (it->second<MinValue)MinValue = it->second;
	
	//cout<<"map size is "<<mechanismMap.size()<<endl;
	
	//find GCD of all multipliers
	if (MinValue >1)
	{
		for (int j=MinValue;j>0;j--)
		{
			bool allAreMultiplesofJ = true;
			for (map<int,int>::iterator it = mechanismMap.begin();it!=mechanismMap.end();it++)
			{
				if (it->second%j!=0)allAreMultiplesofJ=false;
				break;
			}
		
			if (allAreMultiplesofJ)
			{
				//This is the gcd. Now all by it.
				for (map<int,int>::iterator it = mechanismMap.begin();it!=mechanismMap.end();it++)
					it->second=it->second/j;
				break;
			}
		}

	}
				

	return mechanismMap;
}


vector< map<int,int> > Mechanisms::FindDirectMechanisms(string* S, string* S2, bool DistinctOverallRxns)
{
	vector< map<int,int> > Mechs; 
	vector <pair<int, int> > ReactionIndices;//first integer is reaction index, the value is frequency of occurence -- IMPORTANT NOTE: to pull this reaction from AllReactions in the network, refer to this reaction by Index -1
	vector< pair<generated_rxn, int> > Reactions;//the integer gives the frequency of occurence
	vector<string*> ProdAndIntermediates;
	map<int, pair<string*, int> > RxnCounterMap;//says for each reaction which molecule is being eliminated and how many of the reactions forming or consuming it have been considered so far
	map<int, int > StoichiometryMap;// says what the stoichiometry of each of those molecules above is
	RxnCounterMap[0] = pair<string*,int> (S,0);//this says that for the first reaction, S is th molecule being eliminated and so far none of the reactions forming it or consuming it has been considered
	StoichiometryMap[0] = -1;
	ProdAndIntermediates.push_back(S);
	int Mechcount = 0;

	//ofstream logfile("mechanismLogs.txt");
	
	while (!ProdAndIntermediates.empty() )
	{
		partialMechanism M(&Reactions);
		
		if (isComplete(&M))
		{
			if (!M.isNetZero() && Reactions.size()>=DirectConstr.getMinLength() && getTotalCost(&Reactions)<= DirectConstr.getMaxCost() && getTotalCost(&Reactions) >= DirectConstr.getMinCost() &&(M.getFreq(S)>0) /*&& !ContainsInitialReactantsAsProducts(&M)*/ && !(M.getFreq(S2)<0))
			{
				//add into Mechs if it is new. 
				//logfile<<"M is complete and not net zero!"<<endl;
				map<int,int> OneMech;
				

				OneMech = getMechanismMap(ReactionIndices);


							
				bool IsPresent = false;
				bool IsUniqueOverall = true;
				//logfile<<"total number of mechs now is "<<Mechs.size()<<endl;
				for (int k=0;k<Mechs.size();k++)
				{
					map<int,int>::iterator map_it;
					int count = 0;
					vector<pair<generated_rxn, int> > MechRxns;
					for (map_it=Mechs[k].begin();map_it!=Mechs[k].end();map_it++)
					{
						
						MechRxns.push_back(pair<generated_rxn,int>(Network->AllReactions[map_it->first-1],map_it->second));
						if (OneMech.count(map_it->first)>0)
						{
							if (OneMech[map_it->first]>=map_it->second)
								count++;
						}
																		
					}
					if (partialMechanism(&MechRxns).isEqual(&M))
						IsUniqueOverall = false;
					if (count==Mechs[k].size())
					{
							IsPresent = true;
							break;
					}
				}
				//if it is not already present and is unique (if required to be distinct) & satisfies all constraints! 
				if (!IsPresent && ((DistinctOverallRxns &&IsUniqueOverall) || !DistinctOverallRxns)) 
				{
					vector<double> ActEnergy;
					
					
					if (canCalculateActE)
					{
						calculateActEnergy(ReactionIndices, ActEnergy);
						M.setRxnActE(&ActEnergy);
						for (int k = 0;k<ReactionIndices.size();k++)
							RxnActEnergyMap[ReactionIndices[k].first]=ActEnergy[k];//storing the activation energies - useful for printing them out later
					}
					//cout<<"checking to see if the mechanism satisfies constraints"<<endl;
					if (DirectConstr.isSatisfied(M))
					{
						
						Mechcount++;
						//cout<<"found a mechanism, total so far: "<<Mechcount<<endl;
						Mechs.push_back(OneMech);

						for (int k=0;k<ActEnergy.size();k++)
							cout<<ActEnergy[k]<<endl;
					}
					//else cout<<"not adding..."<<endl;
				}
				
			}
			
			if (Reactions.size()>0)
			{
				
				int freq = Reactions.back().first.occurence(ProdAndIntermediates.back());
				for (int i=0;i<Reactions.size()-1;i++)
				{
					
					Reactions[i].second=Reactions[i].second/freq;
					ReactionIndices[i].second = Reactions[i].second;
				}
			}

			
			if (!Reactions.empty())
				Reactions.pop_back();
			if (!ReactionIndices.empty())
				ReactionIndices.pop_back();
		}
		else 
		{
			//determine if consumption or production reaction has to be generated
			//if the counter is less than the total number of reactions available for that species, 
				//find a reaction from the set that produces or consumes it
				//add it and check if the reaction is zero. If zero remove that reaction, else add it and update. 
				//determine which one is to be removed next and update accordingly
			//else
				//remove that species from that vector; 
		//	cout<<"i am here"<<endl;
			pair<multimap<string*,int>::iterator,multimap<string*,int>::iterator> ret;
			int totalCount;
			int howMany;

			//cout<<Reactions.size()<<endl;
			
			if (StoichiometryMap[Reactions.size()]<0)
			{
				ret =Network->MolProductMap.equal_range(ProdAndIntermediates.back());
				totalCount = Network->MolProductMap.count(ProdAndIntermediates.back());
				howMany = -StoichiometryMap[Reactions.size()];
			}
			else
			{
				ret =Network->MolReactantMap.equal_range(ProdAndIntermediates.back());
				totalCount = Network->MolReactantMap.count(ProdAndIntermediates.back());
				howMany = StoichiometryMap[Reactions.size()];
			}

			bool isSpeciesAlreadyinMech =false;

			
			
			
			//we need to get a reaction not already considered yet (hence the totalCount condition below), check if the unique reactions are fewer than or equal to the max length (without taking into consideration the freq of each step - equal to because I can now keep adding more instances of reactions already in the partial mechanism), and check if total cost is less than max cost (note that total cost takes into consideration the reaction frequency) 
			int numberUniqueRxns = getMechanismMap(ReactionIndices).size();
			
			if (totalCount > RxnCounterMap[Reactions.size()].second && numberUniqueRxns<= DirectConstr.getMaxLength() && getTotalCost(&Reactions)<= DirectConstr.getMaxCost())/*the condition with getMechanismMap checks pretty much if the unique reactions number are in total less than max allowed --this specifically is helpful in the case the next intermediate is not new and there's already some reaction consuming it*/
			{

				for (int i =0;i<ProdAndIntermediates.size()-1;i++)
				{
					//checking if the last one is equal to anyone (other than the last of course)
					if (*ProdAndIntermediates[i]==*ProdAndIntermediates[ProdAndIntermediates.size()-1])
					{
						isSpeciesAlreadyinMech = true;
						advance(ret.first, RxnCounterMap[i].second-1);
						RxnCounterMap[Reactions.size()].second=totalCount-1;//this is a hack to prevent other reactions of this species from further consideration
						break;
					}
				}
				if (!isSpeciesAlreadyinMech)//MAYBE I can move this into the for loop as an else statement 
					advance (ret.first,RxnCounterMap[Reactions.size()].second);
						
				//cout<<"in here "<<totalCount<<"  "<<Reactions.size()<<"  "<<RxnCounterMap[Reactions.size()].second<<endl;
				generated_rxn R = Network->AllReactions[(*ret.first).second-1];
				int occurence = R.occurence(ProdAndIntermediates.back());
				RxnCounterMap[Reactions.size()].second++;//Updating the counter - reactions has been increased by one
				int incrementCost = 0;
				incrementCost = Network->Rtlist[R.get_rule()].getCost();//NOTE: not multiplied by howMany even though getTotalCost in general is - See note2 just below
				
				//cout<<occurence<<"  "<<RxnCounterMap[Reactions.size()].second<<endl;
				
				/*bool canAddOneMoreRxn = true; //determining if this reaction can be safely added
				
				if (numberUniqueRxns==DirectConstr.getMaxLength()) //TODO: can be converted into a different function?
				{
					//check if R is already present in Reactions -- if present and the number of reactions already added is equal to max, then we can still add this reaction
					canAddOneMoreRxn = false;
					for (int i =0;i<ReactionIndices.size();i++)
					{
						if (ReactionIndices[i].first == (*ret.first).second)
						{
							canAddOneMoreRxn = true;
							break;
						}
					}
				}*/

				//cout<<"trying to add the reaction "<<R.reactionstring()<<endl;

				if (occurence!=0 && !isReverseAlreadyPresent(&Reactions,&R) && !willFormInternalCycles(R,Reactions) && (getTotalCost(&Reactions) + incrementCost<= DirectConstr.getMaxCost()))
				{
					//cout<<"adding reactions"<<endl;
					
					for (int i=0;i<Reactions.size();i++)
					{
					
						Reactions[i].second=Reactions[i].second*occurence;
						ReactionIndices[i].second = Reactions[i].second;
					
					}
					
					Reactions.push_back(pair<generated_rxn,int>(R, howMany));
					ReactionIndices.push_back(pair<int,int>((*ret.first).second, howMany));
					
					//Note2: Instead of adding this reaction "howMany" times, we add only once. Also, we don't update the multiplicity coefficient of Reactions as well. Anyways, since an intermediate can show up multiple times, there may be more of the same reaction added later on
					//This is important because what if the intermediate has a stoichiometry of +2 in Reactions so far while R has only -1 and hence howMany = 2. 
					//By adding R howMany times, we exclude the possibility of this R and some other R2 that shows up later on together making up the -2 that is required to negate the intermediate
					//Reactions.push_back(pair<generated_rxn,int>(R, 1));
					//ReactionIndices.push_back(pair<int,int>((*ret.first).second, 1));

					
					
	
					
					partialMechanism M2(&Reactions);
					
				
					map<string*,int> M2Intermediates;
					M2Intermediates = getIntermediates(&M2);
					if (M2Intermediates.size()>0)
					{
						map<string*,int>::iterator it;
						string* NextIntermediate = M2Intermediates.begin()->first;
						int Value = 100000000;
						for (it=M2Intermediates.begin();it!=M2Intermediates.end();it++)
						{
							int TCount;
							if (it->second<0)
							{
								TCount = Network->MolProductMap.count(it->first);
							}
							else
								TCount = Network->MolReactantMap.count(it->first);
							
							if (TCount <Value)
							{
								Value = TCount;
								NextIntermediate = it->first;
							}
						}
						
						ProdAndIntermediates.push_back(NextIntermediate);
						RxnCounterMap[Reactions.size()] = pair<string*,int>(NextIntermediate,0);
						StoichiometryMap[Reactions.size()] = M2Intermediates[NextIntermediate];
						//cout<<"next intermediate "<<*NextIntermediate<<"  "<<M2Intermediates[NextIntermediate]<<endl;
						
					}
					
				
				}
				//else logfile<<"skipping reaction"<<endl;
			}
			else
			{
				//cout<<"erasing"<<endl;
				//remove teh rxn counter entry
				RxnCounterMap.erase (Reactions.size());
				
				//remove the stoichiometry entry
				StoichiometryMap.erase(Reactions.size());
				//remove the reactions entry
				ProdAndIntermediates.pop_back();

				if (Reactions.size()>0)
				{
					int freq = Reactions.back().first.occurence(ProdAndIntermediates.back());
					for (int i=0;i<Reactions.size()-1;i++)
					{
					
						Reactions[i].second=Reactions[i].second/freq;
						ReactionIndices[i].second = Reactions[i].second;
					}
					Reactions.pop_back();
					ReactionIndices.pop_back();
				}
				
				
				
				
			
			}		
		}
	
	}

	return Mechs;
}


bool Mechanisms::willFormInternalCycles(generated_rxn R, vector< pair<generated_rxn, int> >& Reactions)
{

	//ASSUMPTIONS: For each intermediate there is exactly one unique reaction that forms it and one that consumes it
	//This means that I have only one choice of reaction for each intermediate I want to eliminate

	//cout<<"checking for cycles"<<endl;
	if (Reactions.size()==0)return false;

	bool intermediateExists = true;
	string* nextIntermediate;

	vector<pair<generated_rxn, int> > newRxnVector;

	newRxnVector.push_back(pair<generated_rxn, int> (R,1));

	map<int,int> RxnsAdded;

	bool isReactionAlreadyPresent = false;//is R already in Reactions

	for (int i =0;i<Reactions.size();i++)
	{
		/*if (R.reactionstring()==Reactions[i].first.reactionstring() && R.get_rule() == Reactions[i].first.get_rule())
		{
			isReactionAlreadyPresent = true;
		}*/
		RxnsAdded[i] = 0;
	}

	while (intermediateExists)
	{
		
		partialMechanism M(&newRxnVector);

		map<string*,int> MIntermediates;
		MIntermediates = getIntermediates(&M);
		if (MIntermediates.size()>0)
		{
			intermediateExists = true;
			nextIntermediate = MIntermediates.begin()->first;

			//if it is a net reactant, find a reaction in Reaction that generates it

			bool AddedRxn = false;
			if (MIntermediates.begin()->second<0)
			{
				for (int i =0;i<Reactions.size();i++)
				{
					if (Reactions[i].first.IsProduct(nextIntermediate))
					{
						
						if (RxnsAdded[i]<Reactions[i].second)
						{
							newRxnVector.push_back(pair<generated_rxn,int> (Reactions[i].first, 1));
							RxnsAdded[i]+=1;
							AddedRxn = true;
							break;
						}
						else return false;//if you have consumed as many of this reactions as max allowed, there's no other way you are going to be able to negate this intermediate 
					}
				}
			}
			
			//if it is a product, find a reaction in Reaction that consumes it
			if (MIntermediates.begin()->second>0)
			{
				for (int i =0;i<Reactions.size();i++)
				{
					if (Reactions[i].first.IsReactant(nextIntermediate))
					{
						if (RxnsAdded[i]<Reactions[i].second)
						{
							newRxnVector.push_back(pair<generated_rxn,int> (Reactions[i].first, 1));
							RxnsAdded[i]+=1;
							AddedRxn = true;
							break;
						}
						else return false;
					}
				}
			}

			if (!AddedRxn) return false;//if no reaction was added, no way this intermediate can be removed! 

		}
		else
		{
			intermediateExists = false;
			if (M.isNetZero())
			{
				//cout<<"will form cycles"<<endl;
				return true;
			}
		}
	}

	return false; 
}
		







void Mechanisms::GenerateDirectMechanismsOlny(ConstrPtr C)
{
	multimap<string*, int>::iterator it;
	
	string previous ="";
	ofstream mechfile;
	
	string molecule_strfile(filename);
	cout<<"direct mechanisms to be stored in "<<molecule_strfile<<endl;
	
	mechfile.open(filename);
	for (it=Network->MolProductMap.begin();it!=Network->MolProductMap.end();it++)
	{
		
		vector<Molecule> M;
		
		bool flag = true;
		
		if (previous==*((*it).first))flag = false;
		else
		{
			previous = *((*it).first);
			Molecule mol(*((*it).first), moleculesize(*((*it).first)));
			
			
			flag = (C)(mol);
			
		}
		if (flag)
		{
			MoleculeMechanismsMap[it->first] = FindDirectMechanisms(it->first,it->first, DistinctDirectMechs);

			for (int i=0;i<MoleculeMechanismsMap[it->first].size();i++)
			{
				map<int,int>::iterator mapit;
				vector<pair<generated_rxn,int> > pMech;
				for (mapit=MoleculeMechanismsMap[it->first][i].begin();mapit!=MoleculeMechanismsMap[it->first][i].end();mapit++)
				{
					
					pMech.push_back(pair<generated_rxn,int>(Network->AllReactions[mapit->first-1], mapit->second));
					mechfile<<Network->AllReactions[mapit->first-1].reactionstring();
					if (mapit->second>1) mechfile<<" x "<<mapit->second;
					double dHrxn = 0.0;
					double dSrxn = 0.0;
					double dGrxn =0.0;
					if (Network->shouldCalcThermo)
					{
						dHrxn = Network->calculateThermoOfRxn(EnthalpyType, Network->AllReactions[mapit->first-1], Network->Temperature);
						dSrxn = Network->calculateThermoOfRxn(EntropyType, Network->AllReactions[mapit->first-1], Network->Temperature);
						dGrxn = dHrxn - Network->Temperature*dSrxn/1000;
						mechfile<<"  "<<dHrxn<<"  "<<dSrxn<<"  "<<dGrxn;
					}
					if (canCalculateActE)
					{
						//NOTE map_it->first starts from 1 too - because that's how it is stored in ReactionIndices. RxnActEnergyMap stores it the same way

						mechfile<<"  "<<RxnActEnergyMap[mapit->first];
					}

					double rule = Network->AllReactions[mapit->first-1].get_rule();
					mechfile<<"  "<<Network->Rtlist.at(rule).getRuleName();
					
					
					mechfile<<endl;
				}
				mechfile<<"------------------------------------------------"<<endl;
				mechfile<<endl;
				mechfile<<partialMechanism(&pMech).printformula()<<endl;
				mechfile<<"------------------------------------------------"<<endl;
			}
		
		}
	}
	mechfile.close();
}
void Mechanisms::GenerateMechanisms(ConstrPtr S)
{
	multimap<string*, int>::iterator it;
	int counter=0;
	string previous ="";
	ofstream mechfile, AllProductsFile;
	
	string molecule_strfile(filename);
	cout<<"mechanisms to be stored in "<<molecule_strfile<<endl;
	
	mechfile.open(filename);
	AllProductsFile.open("AllProductsFromMech.txt");
	for (it=Network->MolProductMap.begin();it!=Network->MolProductMap.end();it++)
	{
		
		vector<Molecule> M;
		
		bool flag = true;
		
		if (previous==*((*it).first))flag = false;
		else
		{
			previous = *((*it).first);
			Molecule mol(*((*it).first), moleculesize(*((*it).first)));
			
			
			flag = (S)(mol);
			
		}
		if (flag && Network->AllMolecules[(*it).first] <= CompleteConstr.getMaxLength())
		{
			int num_mechs = 0;
			cout<<"checking mechanisms for "<<*(it->first)<<endl;
			num_mechs =	GenerateOverallMechanisms(it->first);
			if (num_mechs>0)
			{
				counter++;
				mechfile<<num_mechs<<" number of mechanisms found!"<<endl;
			}
		}
	}

	cout<<"mechanism generated for "<<counter<<" number of molecules of molecules"<<endl;
	set<string> AllProducts;
	for (int i=0;i<AllMechanisms.size();i++)
	{
		vector<pair<generated_rxn,int> > fullMech;
		for (int j=AllMechanisms[i].NumberOfDirectMechs()-1;j>=0;j--)
		{
			pair<string*,int> index = AllMechanisms[i].getDirectMechanism(j);
			int freq = AllMechanisms[i].getFrequency(j);

			map<int,int>::iterator it;
			vector<pair<generated_rxn,int> > pMech;
			for (it = MoleculeMechanismsMap[index.first][index.second].begin();it!=MoleculeMechanismsMap[index.first][index.second].end();it++)
			{
				pMech.push_back(pair<generated_rxn,int>(Network->AllReactions[it->first-1], 1));
				fullMech.push_back(pair<generated_rxn,int>(Network->AllReactions[it->first-1], it->second*freq));
			}

			//mechfile<<partialMechanism(&pMech).printformula()<<"  "<<endl;
			mechfile<<partialMechanism(&pMech).printSMILES()<<"  ";
			if (freq!=1) mechfile<<"x"<<freq<<endl;
			else mechfile<<endl;
			
		}
		mechfile<<"-----------------------------------------------------------------------"<<endl;
		partialMechanism fM(&fullMech);
		mechfile<<fM.printformula()<<endl;
		string SMILESstringForMech;
		SMILESstringForMech= fM.printSMILES();
		mechfile<<SMILESstringForMech<<endl;
		ofstream AthenaFile;
		if (i==0)AthenaFile.open("MechInfoForAthena.txt");
		else AthenaFile.open("MechInfoForAthena.txt", ios::app);

		AthenaFile<<"Details for mechanism "<<i+1<<" is as follows"<<endl;
		
		AthenaFile<<"The overall mechanism is "<<endl;
		AthenaFile<<SMILESstringForMech<<endl;
		AthenaFile<<"------------------------------------------------------------"<<endl;

		fM.MechinfoForAthena();
		AthenaFile.close();

		mechfile<<endl;

		set<string>::iterator it;
		set<string> FullMechProducts = fM.getProducts();
		pair<set<string>::iterator, bool> NewOrNot;
		for (it=FullMechProducts.begin();it!=FullMechProducts.end();it++)
		{
			NewOrNot = AllProducts.insert(*it);
			if (NewOrNot.second)AllProductsFile<<(*it)<<endl;
		}
	}
	
	mechfile.close();
	AllProductsFile.close();
}

int Mechanisms::GenerateOverallMechanisms(string* S)
{
	vector< pair<generated_rxn, int> > Reactions;//the integer gives the frequency of occurence
	vector< pair<OverallRxn,int> > OverallReactions;
	OverallMechanism OM;
	vector<string*> ProdMolecules;
	map<int, int> RxnCounterMap;
	map<int, int > StoichioMap;
	RxnCounterMap[0] = 0;
	int OverallmechanismCount = 0;
	
	StoichioMap[0] = 1;
	ProdMolecules.push_back(S);
	

	
	while (!ProdMolecules.empty())
	{
		partialMechanism M(&Reactions);

		if (FoundOverallMechanism(&M) && M.getFreq(S)>0 && Reactions.size() >= CompleteConstr.getMinLength() && OM.NumberOfDirectMechs()>= CompleteConstr.getMinCycleCount() && getTotalCost(&Reactions)>=CompleteConstr.getMinCost() && getTotalCost(&Reactions) <= CompleteConstr.getMaxCost())
		{
			if (!M.isNetZero() && CompleteConstr.isSatisfied(M))
			{
				//cout<<Reactions.size()<<endl;
				map<string*,pair<int,int> > m = OM.generateStringPairMap();
				if (UniqueMechs.count(m)==0 )
				{
					AllMechanisms.push_back(OM);
					OverallmechanismCount++;
					//cout<<"mechanism found "<<OverallmechanismCount<<endl;
					//OverallRxn OvR(&M);
					//cout<<"mech is "<<M.printformula()<<endl;
					//OM.printMechDetails();
					UniqueMechs.insert(m);
				}
			}
			
			if (Reactions.size()>0)
			{
				pair<string*, int> index = OM.getDirectMechanism(OM.NumberOfDirectMechs()-1);
				int freq = (OverallReactions.back().first.getFreq(ProdMolecules.back()))/OverallReactions.back().second;
				for (int i=0;i<MoleculeMechanismsMap[index.first].at(index.second).size();i++)
					Reactions.pop_back();
				for (int i=0;i<Reactions.size();i++)
				{
					Reactions[i].second=Reactions[i].second/freq;
				}
				
				OM.removeLastDirectMech();
				
				OverallReactions.pop_back();
				for (int i=0;i<OverallReactions.size();i++)
				{
					OverallReactions[i].second=OverallReactions[i].second/freq;
				}

				
			}
			
		}
		else
		{
			
			
			if (MoleculeMechanismsMap.count(ProdMolecules.back())==0)
				MoleculeMechanismsMap[ProdMolecules.back()] = FindDirectMechanisms(ProdMolecules.back(),S,DistinctDirectMechs);

			int totalcount = MoleculeMechanismsMap[ProdMolecules.back()].size();
			int Indivcount = RxnCounterMap[OM.NumberOfDirectMechs()];
			

			if (totalcount > Indivcount && OM.NumberOfDirectMechs()<CompleteConstr.getMaxCycleCount())
			{
				RxnCounterMap[OM.NumberOfDirectMechs()]++;
				map<int,int> directMech = MoleculeMechanismsMap[ProdMolecules.back()].at(Indivcount);
				vector<pair <generated_rxn,int> > rxn;
				map<int,int>::iterator it;
				bool IsReversePresent = false;
				for (it=directMech.begin();it!=directMech.end();it++)
				{
					rxn.push_back(pair<generated_rxn,int>(Network->AllReactions[it->first-1],it->second*StoichioMap[OM.NumberOfDirectMechs()]));
					if (it->second<=0)cout<<"it->second is "<<it->second<<endl;
					
				}
				int addedOrNot = 0;

				if (Reactions.size() + rxn.size() <=CompleteConstr.getMaxLength())
				{
					partialMechanism pM(&rxn);
					if (pM.getFreq(S)>=0)
					{
						OverallRxn oR(&pM);
						int occurence = oR.getFreq(ProdMolecules.back())/StoichioMap[OM.NumberOfDirectMechs()];
						for (int i=0;i<Reactions.size();i++)
						{
							Reactions[i].second=Reactions[i].second*occurence;
						}
						for (int i=0;i<rxn.size();i++)
						{
							Reactions.push_back(rxn[i]);
						}
						for (int i=0;i<OverallReactions.size();i++)
						{
							OverallReactions[i].second=OverallReactions[i].second*occurence;
						}
						OM.insertDirectMech(ProdMolecules.back(),RxnCounterMap[OM.NumberOfDirectMechs()]-1,StoichioMap[OM.NumberOfDirectMechs()]);
						if (StoichioMap[OM.NumberOfDirectMechs()-1]==0)cout<<"StoichiometryMap is zero"<<endl;
						
						OverallReactions.push_back(pair<OverallRxn,int>(oR,StoichioMap[OM.NumberOfDirectMechs()-1]));
						if (oR.getFreq(ProdMolecules.back())==0)
						{
							cout<<"oR is zero!"<<endl;
							cout<<"rxn size is "<<rxn.size()<<endl;
							for (int i=0;i<rxn.size();i++)
							{
								cout<<rxn[i].first.reactionstring()<<"  "<<rxn[i].second<<endl;
							}
							
						}
						partialMechanism M2(&Reactions);

									
						map<string*,int> RemainingSpecies = M2.getStoichiometry();

						map<string*,int>::iterator it2;
						for (it2=RemainingSpecies.begin();it2!=RemainingSpecies.end();it2++)
						{
							
							if (it2->second<0 && Network->InitialReactants.count(it2->first)==0 && RemainingSpecies.size()>0)
							{
															
								ProdMolecules.push_back(it2->first);
								RxnCounterMap[OM.NumberOfDirectMechs()]=0;
								StoichioMap[StoichioMap.size()]=-it2->second;
								addedOrNot = 1;
								break;
							}
						}
						if (M2.isNetZero())
						{
							for (int i=0;i<Reactions.size();i++)
							{
								cout<<Reactions[i].first.reactionstring()<<endl;
							}
						}
					}
				}
			}

				
				
			else
			{
				

				//remove teh rxn counter entry
				RxnCounterMap.erase (RxnCounterMap.size()-1);
				//cout<<"before popping "<<*ProdMolecules.back()<<endl;
				ProdMolecules.pop_back();
				//remove the stoichiometry entry
				StoichioMap.erase(StoichioMap.size()-1);
				//remove the reactions entry
				if (Reactions.size()>0)
				{
					pair<string*, int> index = OM.getDirectMechanism(OM.NumberOfDirectMechs()-1);
					
					int freq = (OverallReactions.back().first.getFreq(ProdMolecules.back()))/OverallReactions.back().second;
				
					for (int i=0;i<MoleculeMechanismsMap[index.first].at(index.second).size();i++)
						Reactions.pop_back();
					for (int i=0;i<Reactions.size();i++)
					{
						Reactions[i].second=Reactions[i].second/freq;
					}
					
					OM.removeLastDirectMech();
					OverallReactions.pop_back();
					for (int i=0;i<OverallReactions.size();i++)
					{
						OverallReactions[i].second=OverallReactions[i].second/freq;
						
						
					}

					
				}
				
				
			}
		}
	}
	cout<<"Mechanisms found: "<<OverallmechanismCount<<endl;
	return OverallmechanismCount;
}
	
bool Mechanisms::FoundOverallMechanism(partialMechanism * p)
{
	bool value = true;
	map<string*,int> m = p->getStoichiometry();
	if (m.size()>0)
	{
		map<string*,int>::iterator it;
		for (it=m.begin();it!=m.end();it++)
		{
			if (it->second<0 && Network->InitialReactants.count(it->first)==0)
			{
				value = false;
				break;
			}
		}
		return value;
	}
	else return false;
}

bool Mechanisms::isReverseAlreadyPresent(std::vector<pair<generated_rxn,int> >* Rxn, generated_rxn* G)
{
	bool returnvalue = false;

	for (int i=0;i<Rxn->size();i++)
	{
		if (Rxn->at(i).first.isReverseReaction(*G))
		{
			returnvalue = true;
			break;
		}
	}
	return returnvalue;
}

int Mechanisms::getTotalCost(vector<pair<generated_rxn,int> >* R)
{
	int TotalCost = 0;

	for (int i=0;i<(*R).size();i++)
	{
		TotalCost+=Network->Rtlist[(*R)[i].first.get_rule()].getCost()*((*R)[i].second);
	}
	
	return (TotalCost);
}

void Mechanisms::SetCompleteMechConstraints(CompleteMechConstraints cM)
{
	CompleteConstr = cM;
}

void Mechanisms::AddKineticFunctions(std::vector<KineticParamPtr> & KineticFunc)
{
	KineticFunctions = &KineticFunc;
	canCalculateActE = true;
}

bool Mechanisms::calculateActEnergy(vector<pair<int,int> >& rxns, vector<double>& actE)
{
	actE.clear();

	for (int i=0;i<rxns.size();i++)
	{		
		double PreExp, ActE, TempIndex, kinValue;
		int rule = Network->AllReactions[rxns[i].first-1].get_rule();
		

		PreExp = 0.0; ActE=0.0; TempIndex = 0.0; kinValue = 0.0;
		bool calcK = false; bool usesBEP = false; bool usesLFER = false; double alpha = 0.0; double beta = 0.0; bool stick = false;

		double dH = Network->calculateThermoOfRxn(EnthalpyType,Network->AllReactions[rxns[i].first-1],Network->Temperature)*1000.0;
		double dS = Network->calculateThermoOfRxn(EntropyType,Network->AllReactions[rxns[i].first-1],Network->Temperature);
		double delN = Network->getDeltaNGasPhase(Network->AllReactions[rxns[i].first-1]);

		KineticsInfo kinInfo(*KineticFunctions,Network->AllReactions[rxns[i].first-1],0.0821,Network->Temperature, dH, dS, delN, 0.0, Network->firstGasSpecies(rxns[i].first-1,0), Network->firstGasSpecies(rxns[i].first-1,1)); 


		if (!kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK,usesBEP,usesLFER, alpha, beta, stick))
			return false;
		actE.push_back(ActE);
	}
	return true;

}

void Mechanisms::SetDistinctDirectMech(bool c)
{
	DistinctDirectMechs = c;
}

void Mechanisms::SetDirectMechConstraints(DirectMechConstraints dM)
{
	DirectConstr = dM;
}


partialMechanism::partialMechanism(vector < pair<generated_rxn,int> >* pMech)
{
	rxn_index = pMech;

	
	for (int i = 0;i < rxn_index->size();i++)
	{
		
		for (int j=0;j<(*rxn_index)[i].first.number_pdcts();j++)
		{
			if (stoichiometry.count((*rxn_index)[i].first.get_products(j))>0)
				stoichiometry[(*rxn_index)[i].first.get_products(j)]+=(*rxn_index)[i].second;
			else
				stoichiometry[(*rxn_index)[i].first.get_products(j)]=(*rxn_index)[i].second;
		}
		for (int j=0;j<(*rxn_index)[i].first.number_reactants();j++)
		{
			if (stoichiometry.count((*rxn_index)[i].first.get_reactants(j))>0)
				stoichiometry[(*rxn_index)[i].first.get_reactants(j)]+=-(*rxn_index)[i].second;
			else
				stoichiometry[(*rxn_index)[i].first.get_reactants(j)]=-(*rxn_index)[i].second;
		}
	}

}

bool partialMechanism::isNetZero()
{
	bool value = true;
	map<string*,int>::iterator it;
	if (stoichiometry.size()>0)
	{
		for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
		{
			if (it->second !=0)
			{
				value = false;
				break;
			}
		}
	}
	else value = false;
	return value;
}

map<string*,int> partialMechanism::getStoichiometry()
{
	return stoichiometry;
}

bool partialMechanism::isEqual(partialMechanism* pMech)
{
	bool isequal = true;

	map<string*,int> pStoich = pMech->getStoichiometry();
	map<string*,int>::iterator it;
	
	for (it=pStoich.begin();it!=pStoich.end();it++)
	{
		if (stoichiometry.count(it->first)==0)
		{
			if (it->second!=0)
			{
				isequal = false;
				break;
			}
		}
		else
		{
			if (stoichiometry[it->first]!=it->second)
			{
				isequal =false;
				break;
			}
		}
	}
	return isequal;
}


string partialMechanism::printformula()
{
	map<string*, int>::iterator it;
	int num_reactants=0;
	int num_products=0;
	string reactants = "";
	string products = "";
	for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
	{
		if (it->second>=1)
		{
			if (num_products>=1)products+=" + ";
			if (it->second>1)
			{
				char out[5];
				
				sprintf(out, "%d", it->second);
				products+=out[0];
				if (it->second>=10)
					products+=out[1];
				products+=" ";
			}
			
			products+=(*(it->first));
			
			num_products++;
		}
		
		if (it->second<0)
		{
			if (num_reactants>=1)reactants+=" + ";
			if (it->second<-1)
			{
				char out[5];
				
				sprintf(out, "%d", -it->second);
				reactants+=out[0];
				if (it->second<-9)
					reactants+=out[1];
				reactants+=" ";
			}
			
			reactants+=(*(it->first));
			num_reactants++;

		}
	}

	return reactants+" --> "+ products;
}


string partialMechanism::printSMILES()
{
	map<string*, int>::iterator it;
	generated_rxn G;
	for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
	{
		for (int i=0;i<=abs(it->second)-1;i++)
		{
			if (it->second>=1)G.add_products(it->first);
			else G.add_reactants(it->first);
		}
	}
	
	return G.reactionstring();
}

int partialMechanism::getFreq(string* S)
{
	if (stoichiometry.count(S)==0)
		return 0;
	else
		return stoichiometry[S];
}

void partialMechanism::MechinfoForAthena()
{
	ofstream file;

	file.open("MechInfoForAthena.txt", ios::app);

	//first list all reactions - TODO unique preferrably

	for (int i=0;i<rxn_index->size();i++)
	{
		file<<i+1<<". "<<rxn_index->at(i).first.reactionstring()<<"	x"<<rxn_index->at(i).second<<endl;
	}
	file<<endl;
	// second list all species

	map<string*, int>::iterator it;
	for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
	{
		file<<distance(stoichiometry.begin(),it)+1<<". "<<*(it->first)<<endl;
	}
	file<<endl;
	
	// list stoich
	for (int i=0; i<rxn_index->size();i++)
	{
		file<<"Stoich(:nc,"<<i+1<<")=(/";

		
		for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
		{
			int IndexOfMol = rxn_index->at(i).first.NetOccurence(it->first);
			if (it!=stoichiometry.begin())
				file<<","<<IndexOfMol;
			else
				file<<IndexOfMol;
		}
		file<<"/)"<<endl;
	}
	file<<endl;

	// list rate expression
	for (int i=0; i<rxn_index->size();i++)
	{
		file<<"r("<<i+1<<")=k("<<i+1<<")";

		for (it=stoichiometry.begin();it!=stoichiometry.end();it++)
		{
			int IndexOfMol = rxn_index->at(i).first.NetOccurence(it->first);
			for (int j=0;j<IndexOfMol;j++)
				file<<"*c("<<distance(stoichiometry.begin(),it)+1<<")";
		}
		file<<endl;
	}
	file<<endl;
	file<<endl;

	file.close();
	

}


int partialMechanism::numMolecule(ConstrPtr C)
{
	int counter =0;
	for (int i=0;i<(*rxn_index).size();i++)
	{
		int number = (*rxn_index)[i].first.number_reactants();
		for (int j=0;j< number ; j++)
		{
			string mol = *(*rxn_index)[i].first.get_reactants(j);
			if ((C)(Molecule(mol,moleculesize(mol))))
			{
				counter+= (*rxn_index)[i].second;
				break;
			}
		}
	}
	return counter;
}

int partialMechanism::numRule(int Rule)
{
	int Counter=0;
	for (int i = 0;i<(*rxn_index).size();i++)
	{
		if ((*rxn_index)[i].first.get_rule()==Rule)Counter+=(*rxn_index)[i].second;
	}
	return Counter;
}

int partialMechanism::numRuleWithParticipants(int Rule , ConstrPtr c, int p)
{
	int counter = 0;
	for (int i=0;i<(*rxn_index).size();i++)
	{
		if ((*rxn_index)[i].first.get_rule()==Rule)
		{
			
			int number;
			if (p==0)number = (*rxn_index)[i].first.number_reactants();
			else number = (*rxn_index)[i].first.number_pdcts();

			for (int j=0;j<number;j++)
			{
				string mol="";
				if (p==0)mol = *(*rxn_index)[i].first.get_reactants(j);
				else mol = *(*rxn_index)[i].first.get_reactants(j);

				if ((c)(Molecule(mol,moleculesize(mol))))
				{
					counter+=(*rxn_index)[i].second;
					break;
				}
			}
		}
	}
	return counter;
}

int partialMechanism::numRuleWithParticipants(int Rule, CombinedConstrPtr c, int p)
{
	int counter = 0;
	for (int i=0;i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].first.get_rule()==Rule)
		{
			if (p==0 && (*rxn_index)[i].first.number_reactants()==2)
			{
				Molecule mol1(*(*rxn_index)[i].first.get_reactants(0), moleculesize(*(*rxn_index)[i].first.get_reactants(0)));
				Molecule mol2(*(*rxn_index)[i].first.get_reactants(1), moleculesize(*(*rxn_index)[i].first.get_reactants(1)));
				if ((c)(mol1,mol2))
				{
					counter+=(*rxn_index)[i].second;
					break;
				}
			}
			if (p!=0 && (*rxn_index)[i].first.number_pdcts()>=2)
			{
				bool satisfied=false;
				for (int j=0;j<(*rxn_index)[i].first.number_pdcts()-1;j++)
				{
					for (int k=j+1;k<(*rxn_index)[i].first.number_pdcts();j++)
					{

						Molecule mol1(*(*rxn_index)[i].first.get_products(j), moleculesize(*(*rxn_index)[i].first.get_products(j)));
						Molecule mol2(*(*rxn_index)[i].first.get_products(k), moleculesize(*(*rxn_index)[i].first.get_products(k)));
						if ((c)(mol1,mol2))
						{
							counter+=(*rxn_index)[i].second;
							satisfied = true;
							break;
						}
					}
				}
				if (satisfied)break;
			}
		}
	}
	return counter;
}

int partialMechanism::numMoleculeOverall(ConstrPtr c)
{
	map<string*, int> m = getStoichiometry();
	int counter = 0;
	if (m.size()>0)
	{
		map<string*,int>::iterator it;
		for (it=m.begin();it!=m.end();it++)
		{
			Molecule mol (*(it->first), moleculesize(*(it->first)));
			if ((c)(mol))
				counter+= it->second;					
		}
		return counter;
		
	}
	else return 0;
}

int partialMechanism::numSatisfying(int r, RxnConstraintPtr rcp)
{
	int counter = 0;
	for (int i =0; i<rxn_index->size();i++)
	{
		if ((*rxn_index)[i].first.get_rule()==r)
		{
			RxnInfo Rinfo(&rxn_index->at(i).first, 0.0, 0.0, 0.0,RxnActE->at(i));
			if ((rcp)(Rinfo))
				counter+=(*rxn_index)[i].second;
		}
	}
	return counter;
}
	
void partialMechanism::setRxnActE(vector<double> * ActE)
{
	RxnActE = ActE;
}


set<string> partialMechanism::getProducts()
{
	map<string*, int>::iterator it;
	set<string> products;
	for (it = stoichiometry.begin();it!=stoichiometry.end();it++)
	{
		if (it->second>0)
			products.insert(*(it->first));
	}
	return products;
}
	 
	


OverallMechanism::OverallMechanism()
{
	MoleculeDirectMechPair.clear();
	Frequency.clear();
}
void OverallMechanism::insertDirectMech(string* S, int index, int freq)
{
	MoleculeDirectMechPair.push_back(pair<string*,int> (S,index));
	Frequency.push_back(freq);
}

int OverallMechanism::NumberOfDirectMechs()
{
	return MoleculeDirectMechPair.size();
}

pair<string*,int> OverallMechanism::getDirectMechanism(int i)
{
	return MoleculeDirectMechPair[i];
}

int OverallMechanism::getFrequency(int i)
{
	return Frequency[i];
}

void OverallMechanism::removeLastDirectMech()
{
	MoleculeDirectMechPair.pop_back();
	Frequency.pop_back();
}

map<string*, pair<int,int> > OverallMechanism::generateStringPairMap()
{
	map<string*,pair<int,int> > m;
	for (int i=0;i<MoleculeDirectMechPair.size();i++)
	{
		m[MoleculeDirectMechPair[i].first]= pair<int,int>(MoleculeDirectMechPair[i].second,Frequency[i]);
	}
	return m;
}

void OverallMechanism::printMechDetails()
{
	for (int i=0;i<MoleculeDirectMechPair.size();i++)
	{
		cout<<*(MoleculeDirectMechPair[i].first)<<"  "<<MoleculeDirectMechPair[i].second<<"  "<<Frequency[i]<<endl;
	}
}


OverallRxn::OverallRxn(partialMechanism* p)
{
	map<string*,int> m = (*p).getStoichiometry();

	map<string*,int>::iterator it;

	for (it=m.begin();it!=m.end();it++)
	{
		if (it->second!=0)
			SpeciesStoichio[it->first]=it->second;
	}
}

int OverallRxn::getFreq(string* S)
{
	if (SpeciesStoichio.count(S)==0)
	{
		cout<<"not found, freq 0"<<endl;
		return 0;
	}
	else
		return SpeciesStoichio[S];
}


DirectMechConstraints::DirectMechConstraints()
{
	MaxLength = 5;
	MaxCost = 10000000;
	MinLength = 0;
	MinCost= 0;
};

void DirectMechConstraints::SetMaxCost(int c)
{
	MaxCost = c;
}

void DirectMechConstraints::SetMaxLength(int c)
{
	MaxLength = c;
}

void DirectMechConstraints::SetMinLengthAndCost(int c, int d)
{
	MinLength = c;
	MinCost = d;
}

bool DirectMechConstraints::isSatisfied(partialMechanism& pM)
{
	return (DMConstraints)(pM);
}

void DirectMechConstraints::AddDMConstraints(DirectMechConstrPtr mC)
{
	DMConstraints = mC;
}


int DirectMechConstraints::getMaxCost()
{
	return MaxCost;
}

int DirectMechConstraints::getMaxLength()
{
	return MaxLength;
}

int DirectMechConstraints::getMinCost()
{
	return MinCost;
}

int DirectMechConstraints::getMinLength()
{
	return MinLength;
}


CompleteMechConstraints::CompleteMechConstraints()
{
	MaxLength= 20;
	MaxCost = 10000000;
	MinLength = 0;
	MinCost = 0;
	MaxCycleCount = 20;
	MinCycleCount = 0;
}

void CompleteMechConstraints::AddMechConstraints(CompleteMechConstrPtr mC)
{
	MechConstraints = mC;
}

void CompleteMechConstraints::SetMaxCycleCount(int c)
{
	MaxCycleCount = c;
}

void CompleteMechConstraints::SetMinCycleCount(int c)
{
	MinCycleCount = c;
}

void CompleteMechConstraints::SetMaxCost(int c)
{
	MaxCost = c;
}


void CompleteMechConstraints::SetMaxLength(int c)
{
	MaxLength = c;
}

void CompleteMechConstraints::SetMinLengthAndCost(int c, int d)
{
	MinLength = c;
	MinCost = d;
}

bool CompleteMechConstraints::isSatisfied(partialMechanism & pM)
{
	return (MechConstraints)(pM);
}

int CompleteMechConstraints::getMaxCycleCount()
{
	return MaxCycleCount;
}

int CompleteMechConstraints::getMaxLength()
{
	return MaxLength;
}

int CompleteMechConstraints::getMaxCost()
{
	return MaxCost;
}

int CompleteMechConstraints::getMinCycleCount()
{
	return MinCycleCount;
}
int CompleteMechConstraints::getMinCost()
{
	return MinCost;
}

int CompleteMechConstraints::getMinLength()
{
	return MinLength;
}
