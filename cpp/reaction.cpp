#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <map>
#include <vector>
#include <utility>
//#include <deque>
//#include <cctype>
//using namespace std;

//#include <stdio.h>
//#include <stdlib.h>

//#include "common.h"
//#include "additionalfunc.h"
#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "singleatom.h"
//#include "compositeatom.h"
//#include "atomcontainer.h"
#include "molecule.h"
//#include "patternmatch.h"
//#include "substructure.h"
//#include "reactiontype.h"
//#include "generated_rxn.h"
#include "reaction.h"
//#include "rng.h"

#include <Python.h>
#include "..\python\pyring\reactiontype.h"

using std::string; 
using std::pair; 
using std::vector; 
using std::map; 
//using std::deque; 
using std::cout; using std::endl; 


Reactiontype::Reactiontype()
{
	reactant_pattern.clear();
	bondchanges.clear();// a vector of bond_operations to store connectivity changes
	mod_atomtype.clear();//a map indicating how an atomtype of an atom with the given label changes
	
	IntraMolecularRxnOnly = false;
	IntraMolecularRxnAlso = false;
	isSelfRxnOnly = false;
	checkRateConstant = false;

	CombinedConstraint = NULL;
	Cost = 0;
	RuleName = "";
	maxSpeciesRank =1000000;
	
		 
}
void Reactiontype::add_reactant_pattern(Substructure S)
{
	
	reactant_pattern.push_back(S);
	
}

void Reactiontype::AllowIntraMolecularRxnOnly()
{
	IntraMolecularRxnOnly = true;
}

void Reactiontype::AllowIntraMolecularRxnAlso()
{
	IntraMolecularRxnAlso = true;
}

void Reactiontype::add_bondchanges(bond_operations B)
{
	bondchanges.push_back(B);
}
void Reactiontype::add_mod_atomtype(int i, std::string S)
{
	string S2="";
	for (int j=0;j<S.length();j++)
	{
		if (S.compare(j,1,"'")!=0)S2+=S[j];
	}


	mod_atomtype[i]=S2;

}

int Reactiontype::get_molecularity()
{
	return reactant_pattern.size();
}
bool Reactiontype::BreaksAromaticity(int j)
{
	int atomnumber = 0;
	int patternindex = -1;
	for (int i=0;i<reactant_pattern.size();i++)
	{
		atomnumber = reactant_pattern[i].findatomwithlabel(j);
		if (atomnumber>0)patternindex =i;
	}
	if (atomnumber<=0)
	{
		for (int z = 0;z<FragmentCopyLabels.size();z++)
		{
			map<int, int>::iterator it;
			for (it=FragmentCopyLabels.at(z).begin();it!=FragmentCopyLabels.at(z).end();it++)
			{
				if ((*it).second==j)
				{
					atomnumber=(*it).first;
					patternindex = FragmentCopyIndex[z];
				}
			}
			if (atomnumber>=0)break;
		}
	}

	if (atomnumber<=0)return false;
	else
	{
		bool returnvalue = false;
		if (!reactant_pattern[patternindex].atoms[atomnumber]->IsSymbolAliphatic()|| (reactant_pattern[patternindex].atom_flag[atomnumber][0]==0))
		{
			map<int,string>::iterator it2;
			it2=mod_atomtype.find(j);
			if (it2!=mod_atomtype.end())
			{
				string c = (*it2).second;
				if (isupper(c[0]))returnvalue = true;
			}
		}
		return returnvalue;
	}
}
void Reactiontype::setCost(int i)
{
	Cost = i;
}

int Reactiontype::getCost()
{
	return Cost;
}

				
void Reactiontype::AddReactantCopy(int Index, map<int,int> Labels)
{
	FragmentCopyIndex.push_back(Index);
	FragmentCopyLabels.push_back(Labels);
}

void Reactiontype::add_reactantconstraint(ConstrPtr CP)
{
	RxnConstraints.push_back(CP);
}

void Reactiontype::add_combined_constraint(CombinedConstrPtr CP)
{
	CombinedConstraint = CP;
}

void Reactiontype::add_productconstraint(ConstrPtr CP)
{
	ProductConstraints = CP;
}

void Reactiontype::disconnect_bond(int i, int j)
{
	bond_operations B;
	B.first_atom_label=i;
	B.second_atom_label=j;
	B.set_bond_order=0;
	B.change_bond_order=0;
	add_bondchanges(B);

}

void Reactiontype::setRuleName(string S)
{
	RuleName = S;
}

string Reactiontype::getRuleName()
{
	return RuleName;
}

void Reactiontype::connect_bond(int i, int j)
{
	bond_operations B;
	B.first_atom_label=i;
	B.second_atom_label=j;
	B.set_bond_order=1;
	B.change_bond_order=0;
	add_bondchanges(B);
}

void Reactiontype::increaseBO(int i, int j, int k)
{
	bond_operations B;
	B.first_atom_label = i;
	B.second_atom_label  = j;
	B.set_bond_order = -1;
	B.change_bond_order = k;
	add_bondchanges(B);
}

void Reactiontype::decreaseBO(int i,int j, int k)
{
	bond_operations B;
	B.first_atom_label = i;
	B.second_atom_label  = j;
	B.set_bond_order = -1;
	B.change_bond_order = -k;
	add_bondchanges(B);
}

int Reactiontype::getFragmentCopyIndex(int i)
{
	return FragmentCopyIndex[i];
}

void Reactiontype::setSpeciesRank(int i)
{
	maxSpeciesRank = i;
}

int Reactiontype::getSpeciesRank()
{
	return maxSpeciesRank;
}

void Reactiontype::AllowSelfRxnOnly()
{
	isSelfRxnOnly = true;
}

void Reactiontype::setMinRateConst(double k)
{
	Min_k_value = k;
	checkRateConstant = true;
}

double Reactiontype::getMinRateConst()
{
	return Min_k_value;
}

bool Reactiontype::shouldCheckRateConstant()
{
	return checkRateConstant;
}

//end of Reactiontype implementation

//start of Reaction
Reaction::Reaction(vector<Molecule>& Mol, Reactiontype & React, vector<Patternmatch>& matches)
{
	cout<<"in here"<<endl;
	M= Mol;
	cout<<"M size: "<<M.size()<<endl;
	
	Rt=React;
	reactpattlist = matches;
	GenerateIntraMolecularRxnsAlso = false;
	GenerateIntraMolecularRxnsOnly = false;
	OnlyOneReactantUsed = false;
	if (Mol.size()==2)
	{
		if (React.IntraMolecularRxnAlso)
		{
			if (Mol[0].moleculestring()==Mol[1].moleculestring())//check for intramolecular rxns only when both reactants are the same. This 
				//prevents unnecessary duplication of intramolecular reactions. We are also assured of atleast one case where we check for A+A kind of reactants. 
				GenerateIntraMolecularRxnsAlso = true;
		}
		if (React.IntraMolecularRxnOnly)
		{
			
			GenerateIntraMolecularRxnsOnly = true;
		}
		
	}
	try{
	generate_reactions();
	}
	catch(pair<string,string> err)
	{
		throw;
	}
	
}

void Reaction::generate_reactions()
{
	/*for (int i=0;i<M.size();i++)
	{
		
		findmatch(i);
		
	}*/

	/*There are many for-loops in the code below with seemingly similar blocks inside. This code needs a little bit of cleaning up. But the basic idea is as follows
	There are many options with respect to a reaction. They are

	1. Only one reactant
	2. Only one reactant and its copy
	3. Two reactants. 
	4. Two reactants and a copy of the first
	5. Two reactants and a copy of the second
	6. Two reactants- but both are same and Intramolecular reaction is allowed! Note that possibility of intramolecular reactions is checked only within Reaction and not in rxn_net_gen for example. 
	So the code below, does 

		For each of unique matches of the first reactant pattern, 
			If the reaction has two reactants 
				If the intramolecular reaction is allowed
					for each of unique matches of the second reactant pattern such that there are no atom overlaps in the molecule with that of hte first reactant
						Find all reactions
				Merge the atom containers of the two reactants
				for each of the unique matches of the second reactant pattern
					If a reactant copy exists
						for each set of unique matches of the appropriate reactant pattern
							find all reactions comprising the three reactants
					Else 
						find all reactions comprising these two reactants
			Else
				If a reactant copy exists
					for each set of unique matches of the first reactant pattern 
						find all reactions comprising the two reactants
				Else 
					find reactions comprising this single reactant

						

					 
	*/
				
	try{
	int Number_reactions=0;
	if (M.size()==2 && (GenerateIntraMolecularRxnsOnly || GenerateIntraMolecularRxnsAlso))
	{
		for (int i=0; i<reactpattlist[0].number_of_matches();i++)
		{
			Atomcontainer A=(Atomcontainer)M[0];//set the atomcontainer first as that of just the first reactant
			

			for (int k=0;k<reactpattlist[0].get_match(i).size();k++)
			{
				int atomnumber=(reactpattlist[0].get_match(i))[k];
				int atomlabel=Rt.reactant_pattern[0].getlabel(k);
				vector<int> atomHlabel = Rt.reactant_pattern[0].getHlabel(k);
				A.setlabel(atomnumber,atomlabel);
				A.setHlabel(atomnumber,atomHlabel);
			
			}

			OnlyOneReactantUsed = true;
			
			for (int j=0;j<reactpattlist[1].number_of_matches();j++)
			{
				Atomcontainer C=A;
				bool IsAtomOverlap = false;

				for (int k=0;k<reactpattlist[1].get_match(j).size();k++)
				{
					int atomnumber = (reactpattlist[1].get_match(j))[k];
					//We need to check if there are overlaps in teh atoms of the first and second pattern and avoid that. Below, we check if the label of the atom is already set or not. A fresh atom will have its label as 0.
					if (C.getlabel(atomnumber)>0)
					{
						IsAtomOverlap = true;
						break;
					}
					else
					{
						
						int atomlabel=Rt.reactant_pattern[1].getlabel(k);
						vector<int> atomHlabel = Rt.reactant_pattern[1].getHlabel(k);
						C.setlabel(atomnumber,atomlabel);
						C.setHlabel(atomnumber,atomHlabel);
					}
				}
				
				if (!IsAtomOverlap)
				{
					vector<pair<string,int> >Hyd_products;
					vector<Molecule> prod_mol;
					bool ChangesWithoutErrors;
					ChangesWithoutErrors = perform_changes(C,Hyd_products);
					Number_reactions = reactpattlist[0].H_factor(i)*reactpattlist[1].H_factor(j);				
					if (ChangesWithoutErrors) addproducts(C,prod_mol,Hyd_products, Number_reactions);
				}
			}
		}
	}
	
	for (int i=0;i<reactpattlist[0].number_of_unique_matches();i++)//find all the unique matches
	{
		Atomcontainer A=(Atomcontainer)M[0];//set the atomcontainer first as that of just the first reactant
		

		for (int k=0;k<reactpattlist[0].get_unique_matches(i).size();k++)
		{
			int atomnumber=(reactpattlist[0].get_unique_matches(i))[k];
			int atomlabel=Rt.reactant_pattern[0].getlabel(k);
			vector<int> atomHlabel = Rt.reactant_pattern[0].getHlabel(k);
			A.setlabel(atomnumber,atomlabel);
			A.setHlabel(atomnumber,atomHlabel);
			
		}	
		
		Number_reactions =0;				
		if (M.size()==2)
		{
			/*if (GenerateIntraMolecularRxnsOnly || GenerateIntraMolecularRxnsAlso)
			{
				OnlyOneReactantUsed = true;
				int counting = 0;
				for (int j=0;j<reactpattlist[1].number_of_unique_matches();j++)
				{
					Atomcontainer C=A;
					bool IsAtomOverlap = false;
					cout<<"second pattern number "<<j<<endl;
					for (int k=0;k<reactpattlist[1].get_unique_matches(j).size();k++)
					{
						int atomnumber = (reactpattlist[1].get_unique_matches(j))[k];
						//We need to check if there are overlaps in teh atoms of the first and second pattern and avoid that. Below, we check if the label of the atom is already set or not. A fresh atom will have its label as 0.

						if (C.getlabel(atomnumber)>0)
						{
							IsAtomOverlap = true;
							break;
						}
						else
						{
							
							int atomlabel=Rt.reactant_pattern[1].getlabel(k);
							vector<int> atomHlabel = Rt.reactant_pattern[1].getHlabel(k);
							C.setlabel(atomnumber,atomlabel);
							C.setHlabel(atomnumber,atomHlabel);
						}
					}
					
					//if (IsAtomOverlap)break;
					if (!IsAtomOverlap)
					{
						
						counting++;
						vector<pair<string,int> >Hyd_products;
						vector<Molecule> prod_mol;
						perform_changes(C,Hyd_products);
						Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor(i);
						Number_reactions = Number_reactions*reactpattlist[1].getMatchFrequency(j)*reactpattlist[1].H_factor(j);						
						addproducts(C,prod_mol,Hyd_products, Number_reactions);
					}
				}
				cout<<"counting is "<<counting<<endl;
				
			}*/
			cout<<"Till here"<<endl;
			if (!GenerateIntraMolecularRxnsOnly && !GenerateIntraMolecularRxnsAlso && check_combined_constraint_for_reactiontype(Rt, M[0], M[1]))
			//if (!GenerateIntraMolecularRxnsOnly && !GenerateIntraMolecularRxnsAlso && (Rt.CombinedConstraint)(M[0], M[1]))//NOTE: this way, effectively, combined constraints cannot be applied to intramolecular reactions!!!!
			{
				Atomcontainer B=(Atomcontainer)M[1];
				A.merge(B);
				OnlyOneReactantUsed = false;

				for (int j=0;j<reactpattlist[1].number_of_unique_matches();j++)
				{
				
					Atomcontainer C=A;
					for (int k=0;k<reactpattlist[1].get_unique_matches(j).size();k++)
					{
						int atomnumber=(reactpattlist[1].get_unique_matches(j))[k];
						int atomlabel=Rt.reactant_pattern[1].getlabel(k);
						vector<int> atomHlabel = Rt.reactant_pattern[1].getHlabel(k);
						C.setlabel(atomnumber+M.at(0).getsize(),atomlabel);
						C.setHlabel(atomnumber+M.at(0).getsize(),atomHlabel);
					}
					if (Rt.FragmentCopyIndex.size()>0)
					{
						/*int totalSize = M.at(0).getsize() + M.at(1).getsize();
						for (int z = 0;z<Rt.FragmentCopyIndex.size();z++)
						{
							int CurrentFragmentCopyIndex = Rt.FragmentCopyIndex.at(i);
							Atomcontainer B2=(Atomcontainer) M[CurrentFragmentCopyIndex];
							C.merge(B2);
						
							for (int j2=0;j2<reactpattlist[CurrentFragmentCopyIndex].number_of_unique_matches();j2++)
							{
								Atomcontainer C2=C;
								for (int k2=0;k2<reactpattlist[CurrentFragmentCopyIndex].get_unique_matches(j2).size();k2++)
								{
									int atomnumber = (reactpattlist[CurrentFragmentCopyIndex].get_unique_matches(j2))[k2];
									int atomlabel = Rt.reactant_pattern[CurrentFragmentCopyIndex].getlabel(k2);
							
									vector<int> atomHlabel = Rt.reactant_pattern[CurrentFragmentCopyIndex].getHlabel(k2);
								
									map<int,int>::iterator it;
								
									it = Rt.FragmentCopyLabels.find(atomlabel);
									atomlabel= (*it).second;
									C2.setlabel(atomnumber+M.at(0).getsize()+M.at(1).getsize(), atomlabel);
									C2.setHlabel(atomnumber+M.at(0).getsize()+M.at(1).getsize(), atomHlabel);
									for (it=Rt.FragmentCopyLabels.begin();it!=Rt.FragmentCopyLabels.end();it++)
									{
										C2.changeHlabel((*it).first,(*it).second);
									}

								}
								vector<pair<string,int> >Hyd_products;
								vector<Molecule> prod_mol;
								perform_changes(C2,Hyd_products);
								Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
								Number_reactions = Number_reactions*reactpattlist[1].getMatchFrequency(j)*reactpattlist[1].H_factor_unique_match(j);
								Number_reactions = Number_reactions*reactpattlist[Rt.FragmentCopyIndex].getMatchFrequency(j2)*reactpattlist[Rt.FragmentCopyIndex].H_factor_unique_match(j2);
								addproducts(C2,prod_mol,Hyd_products, Number_reactions);
							}
						}*/
						
						Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
						Number_reactions = Number_reactions*reactpattlist[1].getMatchFrequency(j)*reactpattlist[1].H_factor_unique_match(j);
						GenerateRxnWithFragmentCopies(C,0,Number_reactions);
					}
					else
					{
				
						vector<pair<string,int> >Hyd_products;
						vector<Molecule>prod_mol;
						perform_changes(C,Hyd_products);
						Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
						Number_reactions = Number_reactions*reactpattlist[1].getMatchFrequency(j)*reactpattlist[1].H_factor_unique_match(j);
						addproducts(C, prod_mol,Hyd_products, Number_reactions);
					}
				}				
			}
		}
		else
		{
			if (Rt.FragmentCopyIndex.size()>0)
			{
				/*Atomcontainer A2=(Atomcontainer)M[0];
				A.merge(A2);

				for (int j=0;j<reactpattlist[0].number_of_unique_matches();j++)
				{
					Atomcontainer C=A;
					for (int k=0;k<reactpattlist[0].get_unique_matches(j).size();k++)
					{
						int atomnumber = (reactpattlist[0].get_unique_matches(j))[k];
						int atomlabel = Rt.reactant_pattern[0].getlabel(k);
						vector<int> atomHlabel = Rt.reactant_pattern[0].getHlabel(k);
						map<int,int>::iterator it;
						it = Rt.FragmentCopyLabels.find(atomlabel);
						atomlabel= (*it).second;
						C.setlabel(atomnumber+M.at(0).getsize(), atomlabel);//sets the label
						C.setHlabel(atomnumber+M.at(0).getsize(), atomHlabel);//sets the H label
						for (it =Rt.FragmentCopyLabels.begin();it!=Rt.FragmentCopyLabels.end();it++)
						{
							C.changeHlabel((*it).first,(*it).second);
						}
					}
					vector<pair<string,int> >Hyd_products;
					vector<Molecule> prod_mol;
					perform_changes(C,Hyd_products);
					Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
					Number_reactions = Number_reactions*reactpattlist[0].getMatchFrequency(j)*reactpattlist[0].H_factor_unique_match(j);
					addproducts(C,prod_mol,Hyd_products, Number_reactions);
				}*/
				Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
				GenerateRxnWithFragmentCopies(A,0,Number_reactions);
			}
			else
			{
				vector<pair<string,int> >Hyd_products;
				vector<Molecule>prod_mol;
				perform_changes(A,Hyd_products);
				Number_reactions = reactpattlist[0].getMatchFrequency(i)*reactpattlist[0].H_factor_unique_match(i);
				
				addproducts(A,prod_mol,Hyd_products, Number_reactions);
				
				A.resetlabel(0);
				
			}
		}		
	}
	}
	catch (pair<string,string> er){
		throw;
	}
}


void Reaction::addproducts(Atomcontainer & A, vector<Molecule>& prod_mol, vector<pair<string, int> >& Hyd_products, int reaction_freq)
{
    
    
    pair< vector<Atomcontainer>, vector< vector <int> > >prod;
    prod=A.connectedcomponents();
    vector<int> ReactantSizes;
    
    ReactantSizes.push_back(M.at(0).getsize());
    if (!OnlyOneReactantUsed)
    {
        if (M.size()==2)
            ReactantSizes.push_back(M.at(1).getsize());
		else ReactantSizes.push_back(0);
        if (Rt.FragmentCopyIndex.size()>0)
        {
			for (int i=0;i<Rt.FragmentCopyIndex.size();i++)
			{
				ReactantSizes.push_back(M.at(Rt.FragmentCopyIndex.at(i)).getsize());
			}
        }
    }
	bool isOneReactantAromatic = false;
	int totalInitialAromaticRings = 0;
	int totalFinalAromaticRings = 0;
	if (M.front().isaromaticmolecule() || M.back().isaromaticmolecule())
	{
		isOneReactantAromatic = true;
		totalInitialAromaticRings+=M.front().NumberOfAromaticRings();
		if (M.size()>1)totalInitialAromaticRings+=M.at(1).NumberOfAromaticRings();
		if (Rt.FragmentCopyIndex.size()>0)
		{
			for (int i=0;i<Rt.FragmentCopyIndex.size();i++)
				totalInitialAromaticRings+=M.at(Rt.FragmentCopyIndex.at(i)).NumberOfAromaticRings();
		}
	}
	bool isOneProductAromatic = false;

	bool HasValencyMismatch = false;
	pair<string,string> ErrorInfo;

    for (int i=0;i<prod.first.size();i++)
    {
       
		Molecule mol((prod.first)[i]);
        mol.unique_smiles();

		if (mol.moleculestring()=="C1(C(C(C(C(C1CC)[{Pt}])([{Pt}])[{Pt}])[{Pt}])[{Pt}])(O_[{Pt}])[{Pt}]")
		{
			mol.print_adjacency_list();
			cout<<"trying this again"<<endl;
			Molecule mol2((prod.first)[i]);
			mol2.unique_smiles(1);
			
		}
            
        pair<string,int> P = mol.checkValencyError();
        if ((P.second == 1 || P.second==2 || P.second==3 || P.second==6))
        {
            HasValencyMismatch = true;
			ErrorInfo = pair<string,string>(P.first,mol.moleculestring());
        }
		
		//cout<<mol.moleculestring()<<endl;
        prod_mol.push_back(mol);  
		if (mol.isaromaticmolecule())isOneProductAromatic = true;
		totalFinalAromaticRings+= mol.NumberOfAromaticRings();
        
    }

	if (HasValencyMismatch)
	{
		//the following if condition checks if there was an aromatic ring broken in the process of reaction to see if there indeed was 
		// a valency mismatch, it was caused because upon broken aromaticitiy resulted in a redistribution of electrons that was meaningless
		//this is possible for example in Furan where protonation cannot happen on the cc bond farthest from oxygen because redistribution of electrons is not possible without having two atoms in the molecule with charge
		//this check however does not catch all cases - the best way is to somehow find out if the atom with mismatch in valency was previously aromatic or not
		
		if (totalFinalAromaticRings< totalInitialAromaticRings)
		{
			//cout<<"reaction rule "<<Rt.getRuleName()<<" appears to break aromaticity, but the resultant molecule has valency mismatches, so I am skipping this particular reaction!"<<endl;
		}
		else
		{
			throw ErrorInfo;
		}
	}
	else
	{
		
		generated_rxn G;
		string reactionMFstring ="";
		int Num_reactants = M.size();
		if (OnlyOneReactantUsed)Num_reactants =1;

		for (int i=0;i<Num_reactants;i++)
		{
			string * molecptr;
			molecptr= StringRegistry::getStringPointer(M.at(i).moleculestring());
			G.add_reactants(molecptr);
			
			int counter = 1;
			for (int z = 0; z < Rt.FragmentCopyIndex.size(); z++)
			{
				
				if (i==Rt.FragmentCopyIndex.at(z))
				{
					G.add_reactants(molecptr);
					counter++;
				}
			}

			if (counter > 1)
			{
				char RxnCount[5];
				sprintf(RxnCount, "%d", counter);
				if (counter>1)reactionMFstring+=RxnCount[0];
				if (counter>9)reactionMFstring+=RxnCount[1];
				if (counter>99)reactionMFstring+=RxnCount[2];
				reactionMFstring+= " "+M.at(i).GetMF();
			}
			else reactionMFstring+=M.at(i).GetMF();

			if (i!=Num_reactants-1) reactionMFstring+=" + ";

		}

		reactionMFstring+=" ---> ";
		
		for (int i=0;i<prod_mol.size();i++)
		{
			string * pdctptr = StringRegistry::getStringPointer(prod_mol[i].moleculestring());
			G.add_products(pdctptr);
			reactionMFstring+=prod_mol[i].GetMF();
			if (i!=prod_mol.size()-1)reactionMFstring+=" + ";
			
			vector<int> AtomsFromOthers;
			AtomsFromOthers.resize(1+Rt.FragmentCopyIndex.size(),0);
			for (int k =0;k<prod.second.at(i).size();k++)
			{
				int atom_index = prod.second.at(i).at(k);
				
				if (atom_index >=ReactantSizes[0])
				{
					int incrementalSize = ReactantSizes[0];
					for (int z = 1;z < ReactantSizes.size();z++)
					{
						
						incrementalSize+=ReactantSizes[z];
						if (atom_index<incrementalSize)
							AtomsFromOthers[z-1]++; //NOTE: this will not lead to vector subscript out of range because ReactantSizes is always set one more than AtomsFromOthers!
					}
				}
			}

			G.setAtomsFromOthers(i,AtomsFromOthers);
		}

		for (int i=0;i<Hyd_products.size();i++)
		{
	        
			string * Hptr = StringRegistry::getStringPointer(Hyd_products[i].first);
			G.add_products(Hptr);
			if (Hyd_products[i].first.compare(0,1,"[")==0 && Hyd_products[i].first.compare("[HH]")!=0)
			{
				reactionMFstring+=Hyd_products[i].first.substr(1,Hyd_products[i].first.length()-2);
			}
			else if (Hyd_products[i].first.compare("[HH]")==0)reactionMFstring+="H2";
			else reactionMFstring+=Hyd_products[i].first;
			if (i!=Hyd_products.size()-1)reactionMFstring+=" + ";
			vector<int> AtomsFromOthers;
			AtomsFromOthers.resize(1+Rt.FragmentCopyIndex.size(),0);
			if (Hyd_products[i].second>=ReactantSizes[0])
			{
				int incrementalSize = ReactantSizes[0];
				for (int z = 1;z < ReactantSizes.size();z++)
				{
					incrementalSize+=ReactantSizes[z];
					if (Hyd_products[i].second<incrementalSize)
						AtomsFromOthers[z]++;
				}
				//Note that for hydrogen molecule, we only take info of one of the donors. the other reactant, if any, is an equal donor.
				//I need to check, but I suspect that this will lead to the parent molecule of the H atom that was first mentioned in the bond formation instruction 
				//to be included as the parent (closer) molecule for Hydrogens. 
			}
			G.setAtomsFromOthers(prod_mol.size()+i,AtomsFromOthers);
		}
	    
		G.set_reactionMF(reactionMFstring);

		int actual_freq = reaction_freq;//calculating the actual frequency. Need to divide by n if they are all equal.
		
		//if (G.ProductsWithSameRadicalType())actual_freq = actual_freq/2;
		if (actual_freq>1)actual_freq = actual_freq/G.IdenticalProductsFactor();
		G.set_frequency(actual_freq);
		if (OnlyOneReactantUsed)G.setIntramolecularity(true);
		else G.setIntramolecularity(false);
		
		gen_rxns.push_back(G);
		//cout<<G.reactionstring()<<endl;
	}
	
    
        
}
void Reaction::findmatch(int i)
{
	Patternmatch reactantmatch(M.at(i),Rt.reactant_pattern[i],0);
	reactantmatch.unique_matches();
	reactpattlist.push_back(reactantmatch);
	
						
}


bool Reaction::perform_changes(Atomcontainer &A,vector<pair<string,int> >& HP)
{
	
	
	vector<int>Hdisconnectlabel;
	for (int i=0;i<Rt.bondchanges.size();i++)
	{
		
		int fatomlabel=Rt.bondchanges[i].first_atom_label;
		int satomlabel=Rt.bondchanges[i].second_atom_label;

		int fatom=A.findatomwithlabel(fatomlabel);
		int satom=A.findatomwithlabel(satomlabel);

		if (fatom>=0 && satom>=0)
		{
			if (Rt.bondchanges[i].set_bond_order==0)
				A.breakbond(fatom,satom);
			if (Rt.bondchanges[i].set_bond_order==1)
			{
				if (A.find_BO(fatom,satom)==-1)	A.formbond(fatom,satom);
				else return false;
			}
			if (Rt.bondchanges[i].set_bond_order==-1)
			{
			
				A.changeBO(fatom,satom,Rt.bondchanges[i].change_bond_order);
				
			}
			
		}
		else if (fatom==-1 && satom==-1)//I am assuming that if the labels are -1, they are definitely Hydrogens! 
		{
			
			
			if (Rt.bondchanges[i].set_bond_order==1)
			{	
				//HP.push_back("H");
				int parentatom1 = A.findatomwithHlabel(fatomlabel);
				int parentatom2 = A.findatomwithHlabel(satomlabel);

				HP.push_back(pair<string,int> ("[HH]", parentatom1));
				
				for (int k=Hdisconnectlabel.size()-1;k>=0;k--)
				{
					if ((Hdisconnectlabel[k]==fatomlabel) || (Hdisconnectlabel[k]==satomlabel))
						Hdisconnectlabel.erase(Hdisconnectlabel.begin()+k);
				}

				//Hydrogen counts are already set when the bond is broken. 

				
			}
		}
		else
		{
			if (Rt.bondchanges[i].set_bond_order==0)
			{
				
				if (fatom==-1)
				{
					Hdisconnectlabel.push_back(Rt.bondchanges[i].first_atom_label);
					A.setHydrogens(satom,(A.getHydrogens(satom)-1));
					A.removeHlabel(satom,fatomlabel);
					
				}
				else
				{
					Hdisconnectlabel.push_back(Rt.bondchanges[i].second_atom_label);
					A.setHydrogens(fatom,(A.getHydrogens(fatom)-1));
					A.removeHlabel(fatom,satomlabel);
				}
				
			}

			if (Rt.bondchanges[i].set_bond_order==1)
			{
				for (int k=Hdisconnectlabel.size()-1;k>=0;k--)
				{
					if ((Hdisconnectlabel[k]==fatomlabel) || (Hdisconnectlabel[k]==satomlabel))
					Hdisconnectlabel.erase(Hdisconnectlabel.begin()+k);
				}
				if (fatom==-1)
				{
					A.setHydrogens(satom,(A.getHydrogens(satom)+1));
					A.addHlabel(satom,fatomlabel);
				}
				else
				{
					A.setHydrogens(fatom,(A.getHydrogens(fatom)+1));
					A.addHlabel(fatom, satomlabel);
				}
				
			
			}
			if (Rt.bondchanges[i].set_bond_order==-1)
			{
				if (Rt.bondchanges[i].change_bond_order==4)//if you had a Hydrogen, it has to be single bonded! so forming a partial bond will entail increase in BO of 4.
				{
					int Hcontaining_atom;
					int heavy_atom;
					
					if (fatom==-1)
					{
						A.setHydrogens(satom,(A.getHydrogens(satom)-1));
						A.removeHlabel(satom,fatomlabel);
						Hcontaining_atom =A.findatomwithHlabel(fatomlabel);
						heavy_atom = satom;
					}
					else
					{
						A.setHydrogens(fatom,(A.getHydrogens(fatom)-1));
						A.removeHlabel(fatom,satomlabel);				
						Hcontaining_atom = A.findatomwithHlabel(satomlabel);
						heavy_atom = fatom;
					}

					A.addAtom("H","H",0,'H',0);
					A.setlabel(A.getsize()-1,fatomlabel);
					if (Hcontaining_atom!=-1)
					{
						A.formbond(A.getsize()-1,Hcontaining_atom);
						A.setHydrogens(Hcontaining_atom,A.getHydrogens(Hcontaining_atom)-1);
					}
					A.formbond(A.getsize()-1,heavy_atom);
					A.changeBO(A.getsize()-1,heavy_atom,4);
					
					
				}
				
			}

						

		}

		
	}

	map<int,string>::iterator it;
	it=Rt.mod_atomtype.begin();
	for (;it!=Rt.mod_atomtype.end();it++)
	{
		int atomnumber=A.findatomwithlabel((*it).first);
		if (atomnumber>=0)
		{
			string c=(*it).second;
			//if the atom is a generic atom like $ or & or X, then we need to find the actual element and replace the wildcards with its symbol
			//This is needed for reactions where the transformations are like & -> &+ indicating that a heteroatom picks up the charge. 
			if ((c.compare(0,1,"$")==0) || (c.compare(0,1,"&")==0) ||(c.compare(0,1,"X")==0))
			{
				string c2="";
				c2=A.getatom(atomnumber)->get_atomtype_name()[0];//take the first character! We assume that atoms are only single
				c[0]=c2[0];
			}
			string d=A.getatom(atomnumber)->get_atom_symbol();
			A.setatomtypename(atomnumber,c);

			//accounting for the correct isotope is done in updatesquarebrackets! check!
			
			A.setatomsymbol(atomnumber,c);
			A.setInitialAtomproperties(atomnumber);
			A.setatomvalency(atomnumber);
		}
		if (atomnumber==-1)
		{
			for (int i=0;i<Hdisconnectlabel.size();i++)
			{
				if ((*it).first==Hdisconnectlabel[i])
				{
					string Hsymbol ="";
					Hsymbol = (*it).second;
					if (Hsymbol.length()>1)Hsymbol="["+Hsymbol+"]";
					
					HP.push_back(pair<string,int> (Hsymbol, A.findatomwithHlabel(Hdisconnectlabel[i])));
				}
			}
		}			
	}

	return true;

}

void Reaction::GenerateRxnWithFragmentCopies(Atomcontainer A, int i, int Numb_rxns)
{
	int CurrentFragmentCopyIndex = Rt.FragmentCopyIndex.at(i);
	Atomcontainer B2=(Atomcontainer) M[CurrentFragmentCopyIndex];
	int Initialsize = A.getsize();
	A.merge(B2);
	for (int j=0;j<reactpattlist[CurrentFragmentCopyIndex].number_of_unique_matches();j++)
	{
		Atomcontainer A2 = A;
		
		for (int k=0;k<reactpattlist[CurrentFragmentCopyIndex].get_unique_matches(j).size();k++)
		{
			int atomnumber = (reactpattlist[CurrentFragmentCopyIndex].get_unique_matches(j))[k];
			int atomlabel = Rt.reactant_pattern[CurrentFragmentCopyIndex].getlabel(k);
						
			vector<int> atomHlabel = Rt.reactant_pattern[CurrentFragmentCopyIndex].getHlabel(k);
			
			map<int,int>::iterator it;
											
			it = Rt.FragmentCopyLabels.at(i).find(atomlabel);
			atomlabel= (*it).second;
			A2.setlabel(atomnumber+Initialsize, atomlabel);
			A2.setHlabel(atomnumber+Initialsize, atomHlabel);
			for (it=Rt.FragmentCopyLabels.at(i).begin();it!=Rt.FragmentCopyLabels.at(i).end();it++)
				A2.changeHlabel((*it).first,(*it).second);
		}
		int Number_reactions = Numb_rxns*reactpattlist[Rt.FragmentCopyIndex.at(i)].getMatchFrequency(j)*reactpattlist[Rt.FragmentCopyIndex.at(i)].H_factor_unique_match(j);
		if (i < Rt.FragmentCopyIndex.size()-1)
			GenerateRxnWithFragmentCopies(A2,i+1,Number_reactions);
		else
		{
			vector<pair<string,int> >Hyd_products;
			vector<Molecule> prod_mol;
			perform_changes(A2,Hyd_products);
			addproducts(A2,prod_mol,Hyd_products, Number_reactions);

		}
					
	}

}



generated_rxn Reaction::get_generated_rxns(int i)
{
	return gen_rxns[i];
}


int Reaction::number_rxns_generated()
{
	return gen_rxns.size();
}

bool Reactiontype::isIntraMolecularAlso()
{
	return IntraMolecularRxnAlso;
}

bool Reactiontype::isIntraMolecularOnly()
{
	return IntraMolecularRxnOnly;
}




