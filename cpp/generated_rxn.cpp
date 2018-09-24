#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <map>
#include <vector>
//#include <utility>
#include <deque>
//#include <cctype>
//using namespace std;

//#include <stdio.h>
//#include <stdlib.h>

//#include "common.h"
#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "singleatom.h"
//#include "compositeatom.h"
//#include "atomcontainer.h"
#include "molecule.h"
#include "patternmatch.h"
#include "substructure.h"
//#include "reactiontype.h"
//#include "generated_rxn.h"
#include "generated_rxn.h"

using std::string; 
//using std::pair; 
using std::vector; 
using std::map; 
using std::deque; 
using std::cout; using std::endl; 



generated_rxn::generated_rxn()
{
	reactants.resize(0);
	products.resize(0);
	frequency = 1;
	RuleIndex = -1;
	IsRxnIntramolecular = false;
}

void generated_rxn::add_reactants(string* S)
{
	reactants.push_back(S);
}

void generated_rxn::add_products(string* S)
{
	products.push_back(S);
}

string* generated_rxn::get_reactants(int i)
{
	return reactants[i];
}

string* generated_rxn::get_products(int i)
{
	return products[i];
}
string generated_rxn::reactionstring()
{
	string rstring="";
	for (int i=0;i<reactants.size();i++)
	{
		rstring+=*get_reactants(i);
		if (i!=reactants.size()-1)
			rstring+=".";
	}
	rstring+=">>";
	
	for (int i=0;i<products.size();i++)
	{
		rstring+=*get_products(i);
		if (i!=products.size()-1)
			rstring+=".";
	}
	return rstring;
}
int generated_rxn::number_pdcts()
{
	return products.size();
}
int generated_rxn::number_reactants()
{
	return reactants.size();
}

int generated_rxn::get_frequency()
{
	return frequency;
}

void generated_rxn::set_frequency(int i)
{
	frequency = i;
}

void generated_rxn::set_reactionMF(std::string s)
{
	reactionMF = s;
}

string generated_rxn::get_reactionMF()
{
	return reactionMF;
}

bool generated_rxn::IsReactant(string* S)
{
	bool result = false;
	for (int i=0;i<reactants.size();i++)
	{
		if ((*S)==(*(reactants[i])))
		{
			result = true;
			break;
		}
	}
	return result;
}

bool generated_rxn::IsProduct(string* S)
{
	bool result = false;
	for (int i=0;i<products.size();i++)
	{
		if ((*S)==(*products[i]))
		{
			result = true;
			break;
		}
	}
	return result;
}

void generated_rxn::set_rule(int i)
{
	RuleIndex = i;
}

int generated_rxn::get_rule()
{
	return RuleIndex;
}

int generated_rxn::occurence(string* Sptr)
{
	return abs(NetOccurence(Sptr));
}
	
int generated_rxn::NetOccurence(std::string * Sptr)
{
	return (countSpecies(products,Sptr)- countSpecies(reactants,Sptr));
}

int generated_rxn::getdeltaN()
{
	return (reactants.size()-products.size());
}


int generated_rxn::countSpecies(deque<string*>& deq, string* Sptr)
{
	int result = 0;
	for (int i=0;i<deq.size();i++)
	{
		if ((*Sptr)==(*(deq[i])))
		{
			result++;
			
		}
	}
	return result;
}

int generated_rxn::NetPatternDiff(string S)
{
	Substructure Sub(S, patternsize(S));

	int netDiff = 0;

	for (int i =0;i<products.size();i++)
	{
		Molecule mol (*products[i], moleculesize(*products[i]));
		mol.unique_smiles();
		Patternmatch P(mol,Sub,0);

		netDiff+=P.GetDistinctMatches();
	}

	for (int i = 0;i<reactants.size();i++)
	{
		Molecule mol (*reactants[i], moleculesize(*reactants[i]));
		mol.unique_smiles();
		Patternmatch P(mol,Sub,0);

		netDiff-=P.GetDistinctMatches();
	}

	return netDiff;
}

		

int generated_rxn::getReactantStoich(string* Sptr)
{
	return countSpecies(reactants, Sptr);
}

int generated_rxn::getProductStoich(string* Sptr)
{
	return countSpecies(products,Sptr);
}



void generated_rxn::setAtomsFromOthers(int a, vector<int> b)
{
    AtomsFromOtherReactants[a]=b;
}

vector<int> generated_rxn::getAtomsFromOthers(int a)
{
    vector<int> v;
    if (AtomsFromOtherReactants.count(a)>0)
        return AtomsFromOtherReactants[a];
    else return v;
}

void generated_rxn::setParentMolecule(std::map<int,int> m)
{
    ParentMolecule = m;
	
}

int generated_rxn::getParentMolecule(int a)
{
    if (ParentMolecule.count(a)>0)
        return ParentMolecule[a];
    else 
	{
		cout<<"oh no, returning 0!"<<endl;
		return 0;
	}
}
int generated_rxn::getParentMolecule(string* S)
{
    int k=0;
    for (int i=0;i<products.size();i++)
    {
        if ((*S)==(*products[i]))
        {
            k=i;
            break;
        }
    }
    return getParentMolecule(k);
}

vector<string*> generated_rxn::getDaughterMolecules(string* S)
{
	vector<string*> daughters;
	for (int i=0;i<reactants.size();i++)
	{
		if ((*S)==(*reactants[i]))
		{
			map<int,int>::iterator it;
			for (it=ParentMolecule.begin();it!=ParentMolecule.end();it++)
			{
				if (it->second==i)
					daughters.push_back(products[it->first]);
			}
			break;
		}
	}
	return daughters;
}
				

bool generated_rxn::isSameNetReaction(generated_rxn & G)
{
	bool returnvalue = true;
	if (number_reactants()!=G.number_reactants() || number_pdcts()!=G.number_pdcts())
	{
		returnvalue = false;
	}
	else
	{
		for (int i=0;i<number_reactants();i++)
		{
			if (!G.IsReactant(get_reactants(i)))
			{
				returnvalue = false;
				break;
			}
		}
		if (returnvalue)
		{
			for (int i=0;i<number_pdcts();i++)
			{
				if (!G.IsProduct(get_products(i)))
				{
					returnvalue = false;
					break;
				}
			}
		}
	}
	return returnvalue;
}

bool generated_rxn::isReverseReaction(generated_rxn& G)
{
	bool returnvalue = true;

	

	if (number_reactants()!=G.number_pdcts() || number_pdcts()!=G.number_reactants())
	{
		returnvalue = false;
	}
	else
	{
		for (int i=0;i<number_pdcts();i++)
		{
			if (!G.IsReactant(get_products(i)))
			{
				returnvalue = false;
				break;
			}
		}
		if (returnvalue)
		{
			for (int i=0;i<number_reactants();i++)
			{
				if (!G.IsProduct(get_reactants(i)))
				{
					returnvalue = false;
					break;
				}
			}
		}
	}
	
	return returnvalue;
}

bool generated_rxn::isReactionIntramolecular()
{
	return IsRxnIntramolecular;
}

void generated_rxn::setIntramolecularity(bool b)
{
	IsRxnIntramolecular = b;
}
	

int generated_rxn::netMassDiff()
{
	int netMassDiff = 0;

	for (int i =0;i<reactants.size();i++)
	{
		Molecule mol(*reactants[i], moleculesize(*reactants[i]));
		netMassDiff+=mol.MolecularWeight();
	}

	for (int i=0;i<products.size();i++)
	{
		Molecule mol(*products[i],moleculesize(*products[i]));
		netMassDiff-=mol.MolecularWeight();
	}

	return netMassDiff;
}

bool generated_rxn::ProductsWithSameRadicalType()
{
	//I check the count of each type of radical atomtype in the reactants and products and take the difference. 
	//in the new radicals formed, I check if there is one which is a multiple of 2, if yes, return true, else false; 
	//TODO: not complete, for e.g. does not consider the case when the product leads to four new radicals such that there are two new radical types (I don't know if I have to divide by 4 in this case? )
	
	map<string, int> radicalAtomTypes;
	for (int i=0;i<products.size();i++)
	{
		Molecule mol(*products[i],moleculesize(*products[i]));
		for (int i =0;i<mol.getsize(); i ++)
		{
			if (mol.getatom(i)->get_up()==1)
			{
				if (radicalAtomTypes.count(mol.getatomtype(i))>0)
					radicalAtomTypes[mol.getatomtype(i)]+=1;
				else radicalAtomTypes[mol.getatomtype(i)]=1;
			}
		}

	}

	for (int i=0;i<reactants.size();i++)
	{
		Molecule mol(*reactants[i],moleculesize(*reactants[i]));
		for (int i =0;i<mol.getsize(); i ++)
		{
			if (mol.getatom(i)->get_up()==1)
			{
				if (radicalAtomTypes.count(mol.getatomtype(i))>0)
					radicalAtomTypes[mol.getatomtype(i)]-=1;
				else radicalAtomTypes[mol.getatomtype(i)]=-1;
			}
		}

	}

	//now check if there are net new radicals and in multiples of 2. 

	for (map<string,int>::iterator it = radicalAtomTypes.begin(); it!=radicalAtomTypes.end();it++)
	{
		if (it->second>0 && it->second%2 ==0)
			return true;
	}
	return false;

}

int generated_rxn::IdenticalProductsFactor()
{
	map<string,int> productCounts;
	for (int i=0;i<products.size();i++)
	{	
		if (productCounts.count(*products[i])>0)
			productCounts[*products[i]]+=1;
		else productCounts[*products[i]]=1;
	}

	int prodFactor = 1;

	//first the map should be even-sized and all entries should be equal.
    if (productCounts.size()%2==0)
	{
		int value = productCounts.begin()->second;
		for(map<string,int>::iterator it = productCounts.begin(); it!=productCounts.end();it++)
		{
			if (it->second==value)
				prodFactor*= it->second;
			else 
			{
				prodFactor = 1;
				break;
			}
		}
	}

	return prodFactor;
}



//end of generated_rxns implementation

