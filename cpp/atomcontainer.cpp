#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <map>
#include <vector>
//#include <cctype>
#include <queue>
#include <algorithm>
#include <set>
#include <utility>

//#include <stdio.h>
//#include <stdlib.h>
#include <cmath>

//#include "common.h"
#include "additionalfunc.h"
#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
#include "atomcontainer.h"
//#include "molecule.h"
//#include "patternmatch.h"
//#include "substructure.h"
//#include "reactiontype.h"

//using namespace std;
using std::fabs; 
using std::map; using std::vector; using std::string; using std::set; 
using std::pair; using std::queue; using std::sort;
using std::cout; using std::endl;

//start of Atomcontainer implementation
Atomcontainer::Atomcontainer(int m)//constructor of Atomcontainer class
{
	atoms.reserve(m);//reserving a size of m
	Adjacency.resize(m);//Adjacency list
	BO.resize(m);;//Bond Order
	NN.resize(m);//Number of neighbors
	value.resize(m);//value of each atom (used to evaluate rank)
	rank.resize(m);//rank
	classRank.resize(m);//Class Rank
	Hydrogens.resize(m, 0);//the number of Hydrogens
	size=m;//size of the container
	label.resize(m, -1);//label of an atom. -1 indicates no labelling
	aromatic_atoms.clear();
	H_label.resize(m);
	
}

Atomcontainer::Atomcontainer()
{
	Adjacency.resize(0);
	BO.resize(0);
	size=0;
	NN.resize(0);
	value.resize(0);
	label.resize(0);
	rank.resize(0);
	classRank.resize(0);
	Hydrogens.resize(0);
	aromatic_atoms.clear();
	H_label.resize(0);
}

Atomcontainer::~Atomcontainer()
{
	
	/*vector<Atom*>::iterator it;
	for (it=atoms.begin();it!=atoms.end();it++)
	{
		delete *it;
	}*/

	for (int i=0;i<atoms.size();i++)
	{
		delete atoms[i];
	}
	
	
	atoms.clear();
	BO.clear();
	Adjacency.clear();
	
	
	NN.clear();
	value.clear();
	label.clear();
	
	rank.clear();
	
	Hydrogens.clear();
	classRank.clear();
	aromatic_atoms.clear();
	H_label.clear();
	//cout<<"cleared things up"<<endl;
}

Atomcontainer::Atomcontainer(const Atomcontainer &a)
{
	for (int i=0;i<a.size;i++)
	{
		atoms.push_back(a.atoms[i]->clone());
		atoms[i]->set_nature(a.atoms[i]->get_nature());
	}
	
	Adjacency=a.Adjacency;
	BO=a.BO;
	size=a.size;
	NN=a.NN;
	value=a.value;
	label=a.label;
	rank=a.rank;
	classRank =a.classRank;
	Hydrogens=a.Hydrogens;
	aromatic_atoms = a.aromatic_atoms;
	H_label = a.H_label;

	Allrings = a.Allrings;
}

Atomcontainer& Atomcontainer::operator =(const Atomcontainer &a)
{
	if (this != &a)
	{
		for (int i=0;i<atoms.size();i++)
		{
			delete atoms[i];
		}
		atoms.clear();
		
		for (int i=0;i<a.size;i++)
		{
			atoms.push_back(a.atoms[i]->clone());
			atoms[i]->set_nature(a.atoms[i]->get_nature());
		}
		Adjacency=a.Adjacency;
		BO=a.BO;
		size=a.size;
		NN=a.NN;
		value=a.value;
		label=a.label;
		rank=a.rank;
		classRank = a.classRank;
		Hydrogens=a.Hydrogens;
		aromatic_atoms=a.aromatic_atoms;
		H_label = a.H_label;
		Allrings = a.Allrings;
		
	}
	return *this;
}


void Atomcontainer::setrank(vector<int>& rank_vector)//setting rank
{
	rank=rank_vector;
}
int Atomcontainer::getBplusNBatomcount(int i) const
{
	return Adjacency[i].size();
}

int Atomcontainer::getrank(int r) const //retrieving rank
{
	return rank[r];
}
int Atomcontainer::getclassRank(int r) const
{
	return classRank[r];
}

/* -- removed - new one below
void Atomcontainer::breakTies()
{
	set<int> AtomRanks;
	for (int i=0;i<size;i++)
	{
		AtomRanks.insert(rank[i]);
	}
	set<int>::iterator it;
	for (it=AtomRanks.rbegin();it!=AtomRanks.rend();it++)
	{
		int count =0;
		vector<int>Equalranklist;
		for (int i=0;i<size;i++)
		{
			
			if (rank[i]==(*it))
			{
				count++;
				if (count>1)Equalranklist.push_back(i);
			}
		}
		if (count>1)
		{
			for (int i=0;i<Equalranklist.size();i++)
			{
				rank[Equalranklist[i]]+=i+1;
			}
			//cout<<"break ties called on "<<*it<<endl;
			break;
		}
	}
	//cout<<"done"<<endl;
}*/

//modified by SR - Sep 11 2013; following J. Chem. Inf. Comput. Sci. 2003, 43, 735-742
void Atomcontainer::breakTies()
{
	set<int> AtomRanks;
	for (int i=0;i<size;i++)
	{
		AtomRanks.insert(rank[i]);
	}
	set<int>::iterator it;
	//find the smallest rank that is tied
	vector<int> AtomsWithSmallestTiedRank;
	int SmallestTiedRank = 0;
	int count =0;
	for (it=AtomRanks.begin();it!=AtomRanks.end();it++)
	{
		AtomsWithSmallestTiedRank.clear();
		SmallestTiedRank = 0;
		count = 0;
		
		for (int i=0;i<size;i++)
		{
			
			if (rank[i]==(*it))
			{
				count++;
				AtomsWithSmallestTiedRank.push_back(i);
				SmallestTiedRank = rank[i];
				
			}
		}
		if (count>1)break;
	}
	if (count>1)
	{
	    //cout<<"tie broken "<<SmallestTiedRank<<endl;
		for (int i=0;i<size;i++)
		{
			if (rank[i]>=SmallestTiedRank && i!= AtomsWithSmallestTiedRank[0]) //increment every rank higher than smallest tied rank by 1, except the first atom with the smallest tied rank
				rank[i]++;
		}
		
	}
	
	//cout<<"done"<<endl;
}
int Atomcontainer::distinctRankCount() const
{
	int uniquerankcount =1;
	for (int i=1;i<size;i++)
	{
		int marker=0;
		for (int j=0;j<i;j++)
		{
			if (rank[j]==rank[i])
			{
				marker=1;
				break;
			}
		}
		if (marker==0)uniquerankcount++;
	}

	return uniquerankcount;
}


void Atomcontainer::evaluate_rank()
{
	// This function evaluates the rank of each of the atoms in the atomcontainer based on value
	// Note that by this method if two atoms are of rank 1, say, the next lower rank will be 3 and not 2. 
	
	for(int i=0;i<size;i++)
	{
		rank[i]=1;
	}
	for (int i=0;i< size;i++)
	{	for (int j=0;j<=i;j++)
		{	
			if (fabs(value[j]-value[i])>0.1)
			{
				if ((value[j]-value[i])>0.2)rank[i]++;//this set of commands ensures the above mentioned feature of the function
				else rank[j]++;
			}
		}
	}
}

void Atomcontainer::EvaluateInitialValue()
{
	vector <string> compositeatomlist;
	for (int i=0;i<size;i++)
	{
		if (atoms[i]->get_nature()==1)
		{
			compositeatomlist.push_back(atoms[i]->get_element_name());
			//cout<<atoms[i]->get_element_name()<<endl;
		}
	}
	sort(compositeatomlist.begin(),compositeatomlist.end());
	//for (int i =0;i<compositeatomlist.size();i++)
//	{
//		cout<<compositeatomlist[i]<<"  ";
//	}
//	cout<<endl;
		
	for (int i=0;i<size;i++)
	{

		if (atoms[i]->get_nature()==1)
		{
			for (int k=0;k<compositeatomlist.size();k++)
			{
				if (compositeatomlist[k]==atoms[i]->get_element_name())
				{
					value[i]+=(200+k)*1.0e10;
					break;
				}
			}
		}
			
		if (atoms[i]->get_nature()==0)
		{
					
			if (atoms[i]->get_isotope_number()>0)// calculating the appropriate mass number accounting for isotopes
				value[i]=atoms[i]->get_isotope_mass_number()*1.0e9;
			else
				value[i]=atoms[i]->get_n_mass_number()*1.0e9;

			value[i]=value[i]+atoms[i]->get_atomic_number()*1.0e8;
					
			if (isaromatic(i))
				value[i]=value[i]+1.0e7;
		}

		value[i]=value[i]+ atoms[i]->get_valency()*1.0e6+NN[i]*1.0e5;// Using Nearest neighbor count and valency for the initial value
		
		if (atoms[i]->get_charge()>=1)// using charge (if charge>=1, then charge sign is 1, and abs value is 1)
			value[i]=value[i]+1.0e3 + atoms[i]->get_charge()*1.0e2;
		if (atoms[i]->get_charge()<0)// if charge is negative, then charge sign is 2 and the abs value is 1)
			value[i]=value[i]+2.0e3 + (atoms[i]->get_charge()*1.0e2);
		// if charge is zero, no calculation is required anyways
		//Next we count the double bonds and add it
			
		value[i]=value[i]+dbcount(i)*1.0e1;// using double bond count in value calculation
		value[i]=value[i]+Hydrogens[i];// Finally, adding the number of Hydrogens to get the initial value
	}
}

void Atomcontainer::EvaluateInitialValue(int z)
{
	vector <string> compositeatomlist;
	for (int i=0;i<size;i++)
	{
		if (atoms[i]->get_nature()==1)
		{
			compositeatomlist.push_back(atoms[i]->get_element_name());
			cout<<atoms[i]->get_element_name()<<endl;
		}
	}
	sort(compositeatomlist.begin(),compositeatomlist.end());
	for (int i =0;i<compositeatomlist.size();i++)
	{
		cout<<compositeatomlist[i]<<"  ";
	}
	cout<<endl;
		
	for (int i=0;i<size;i++)
	{

		if (atoms[i]->get_nature()==1)
		{
			for (int k=0;k<compositeatomlist.size();k++)
			{
				if (compositeatomlist[k]==atoms[i]->get_element_name())
				{
					value[i]+=(200+k)*1.0e10;
					break;
				}
			}
		}
			
		if (atoms[i]->get_nature()==0)
		{
					
			if (atoms[i]->get_isotope_number()>0)// calculating the appropriate mass number accounting for isotopes
				value[i]=atoms[i]->get_isotope_mass_number()*1.0e9;
			else
				value[i]=atoms[i]->get_n_mass_number()*1.0e9;

			value[i]=value[i]+atoms[i]->get_atomic_number()*1.0e8;
					
			if (isaromatic(i))
				value[i]=value[i]+1.0e7;
		}

		value[i]=value[i]+ atoms[i]->get_valency()*1.0e6+NN[i]*1.0e5;// Using Nearest neighbor count and valency for the initial value
		
		if (atoms[i]->get_charge()>=1)// using charge (if charge>=1, then charge sign is 1, and abs value is 1)
			value[i]=value[i]+1.0e3 + atoms[i]->get_charge()*1.0e2;
		if (atoms[i]->get_charge()<0)// if charge is negative, then charge sign is 2 and the abs value is 1)
			value[i]=value[i]+2.0e3 + (atoms[i]->get_charge()*1.0e2);
		// if charge is zero, no calculation is required anyways
		//Next we count the double bonds and add it
			
		value[i]=value[i]+dbcount(i)*1.0e1;// using double bond count in value calculation
		value[i]=value[i]+Hydrogens[i];// Finally, adding the number of Hydrogens to get the initial value
	}
}

bool Atomcontainer::isaromatic(int n) const
{
	int flag=0;
	for (int i=0;i<aromatic_atoms.size();i++)
	{
		if (aromatic_atoms[i]==n)
		{
			flag=1;
			break;
		}
	}
	if (flag==1)
		return true;
	else 
		return false;
}

void Atomcontainer::EvaluateValue()
{
	for (int i=0;i<size;i++)
	{
		value[i]=1;			
		for (int j=0;j<NN[i];j++)
		{
			value[i]=value[i]*prime(rank[Adjacency[i][j]]-1);
			if (BO[i][j]>1 && BO[i][j]<4)
			{
				value[i]=value[i]*prime(rank[Adjacency[i][j]]-1);
			}
			if (BO[i][j]==3)
			{
				value[i]=value[i]*prime(rank[Adjacency[i][j]]-1);
			}
			
			//evaluating value as product of primes corresponding to rank of neighbours 
		}
	}
}


vector<set<int> > Atomcontainer::EvaluateAtomClasses()
{
	vector<set<int> > AC;
	vector<bool>IsCounted;
	IsCounted.resize(size,false);

	for (int i=0;i<size;i++)
	{
		
		if (!IsCounted[i])
		{
			set<int> NewSet;
			NewSet.insert(i);
			for (int j=i+1;j<size;j++)
			{
				if (rank[i]==rank[j])
				{
					NewSet.insert(j);
					IsCounted[j] = true;
				}
			}
			AC.push_back(NewSet);
		}
	}

	return AC;
}

bool Atomcontainer::IsEqualClasses(std::vector<set<int> > A, std::vector<set<int> > B) const
{
	int counter =0;
	if (A.size()==B.size())
	{

		for (int i=0;i<A.size();i++)
		{
			for (int j=0;j<B.size();j++)
			{
				if (A[i]==B[j])
				{
					counter++;
					break;
				}
			}
		}
		if (counter==A.size())return true;
		else return false;
	}
	else return false;
}
			

void Atomcontainer::reverse_ranks()
{
	for (int i=0;i<size;i++)
	{
		rank[i]=size+1-rank[i];
	}
}

void Atomcontainer::push_to_first(int i,int j)// this function pushes an element of hte adjacency list into the first spot
{
	int swap1,swap2;
	swap1=Adjacency[i][j];
	swap2=BO[i][j];
	for (int k=j-1;k>=0;k--)
	{
		Adjacency[i][k+1]=Adjacency[i][k];
		BO[i][k+1]=BO[i][k];
	}
	Adjacency[i][0]=swap1;
	BO[i][0]=swap2;
}


void Atomcontainer::sort_adjacency()
{
	//this function sorts the adjacency list in descending order of ranks
	// and also puts the double bonds in the first. this ensures atoms with higher ranks are given higher preference in the DFS
	// and atoms which are connected by a double bond are given the highest preference
	int swap1, swap2;
	
	for (int i=0;i<size;i++)
	{
		for (int j=1;j<NN[i];j++)
		{
	
			for (int k=0;k<NN[i]-j;k++)
			{
	
				if (rank[Adjacency[i][k]]>rank[Adjacency[i][k+1]])
				{
					//swapping part of the Rank sort algorithm
					swap1=Adjacency[i][k];
					swap2=BO[i][k];
					Adjacency[i][k]=Adjacency[i][k+1];
					BO[i][k]=BO[i][k+1];
					Adjacency[i][k+1]=swap1;
					BO[i][k+1]=swap2;
				}
			}
		}
	
		//the loop for transfering the double bonded atoms to the first. Note that the following code ensures that, in case
		// of multiple double bonds, the atoms are in turn sorted according to rank as the loop is done twice in case of two
		// or more double bonds! 
		int db_count=0;
		for (int j=1;j<NN[i];j++)
		{
			if (BO[i][j]==2)
			{
				
				push_to_first(i,j);
				db_count++;
			}
		}
		if (db_count>1)
		{
			for (int j=1;j<db_count;j++)
			{
				push_to_first(i,j);
			}
		}
	}
}

int Atomcontainer::get_adjacency(int i,int j) const//implementation of adjacency list
{
	return Adjacency[i][j];
}
int Atomcontainer::get_BO(int i,int j) const//implementation of BO
{
	return BO[i][j];
}
int Atomcontainer::getsize() const
{
	return size;
}
int Atomcontainer::getNN(int n) const
{
	return NN[n];
}

void Atomcontainer::setBO(int i, int j, int k)
{
	for (int a=0;a<NN[i];a++)
	{
		if (Adjacency[i][a]==j)
		{
			BO[i][a]=k;
			break;
		}
	}
	for (int a=0;a<NN[j];a++)
	{
		if (Adjacency[j][a]==i)
		{
			BO[j][a]=k;
			break;
		}
	}
	
}
int Atomcontainer::dbcount(int n) const//returns the number of double bonds
{
	int dbcount;
	dbcount=0;

	
	for (int i=0;i<Adjacency[n].size();i++)
	{
		if (BO[n][i]==2)dbcount++;
	}

	return dbcount;
}

int Atomcontainer::tpcount(int n) const//returns the number of triple bonds
{
	int tpcount;
	tpcount=0;
	for (int i=0;i<NN[n];i++)
	{
		if (BO[n][i]==3)tpcount++;
	}

	return tpcount;
}

int Atomcontainer::find_BO(int a,int b) const
{
	int value=-1;
	if (a>=0 && a < Adjacency.size())
	{
		for (int i=0;i<NN[a];i++)
		{
			if (Adjacency[a][i]==b)
			{
				value=get_BO(a,i);
				break;
			}
		}
	}
	return value;
}

string Atomcontainer::getatomtype(int i) const
{
	string s;
	s=atoms[i]->get_atomtype_name();
	return s;
}

void Atomcontainer::erase(int i)
{
	delete atoms[i];
	atoms.erase(atoms.begin()+i);
	NN.erase(NN.begin()+i);
	value.erase(value.begin()+i);
	Adjacency.erase(Adjacency.begin()+i);
	BO.erase(BO.begin()+i);
	rank.erase(rank.begin()+i);
	label.erase(label.begin()+i);
	Hydrogens.erase(Hydrogens.begin()+i);
	classRank.erase(classRank.begin()+i);
	H_label.erase(H_label.begin()+i);

	for (int j=0;j<aromatic_atoms.size();j++)
	{
		if (aromatic_atoms[j]==i)
		{
			aromatic_atoms.erase(aromatic_atoms.begin()+j);
			break;
		}
	}

	size--;
	
	for (int j=0;j<size;j++)
	{
		for (int k=0;k<Adjacency[j].size();k++)
		{
			if (Adjacency[j][k]==i)Adjacency[j].erase(Adjacency[j].begin()+k);
		}
	}

	for (int j=0;j<size;j++)
	{
		for (int k=0;k<Adjacency[j].size();k++)
		{
			if (Adjacency[j][k]>i)
				Adjacency[j][k]--;
		}
	}

	for (int j=0;j<aromatic_atoms.size();j++)
	{
		if (aromatic_atoms[j]>i)aromatic_atoms[j]--;
	}
	
}

void Atomcontainer::print_adjacency_list() const
{
	cout<<"printing adjacency list"<<endl;
	for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		
		cout<<i<<" ("<<NN[i]<<")"<<"("<<atoms[i]->get_valency()<<")"<<" ("<<atoms[i]->get_atom_symbol()<<") ("<<Hydrogens[i]<<")"<<rank[i]<<"  "<<classRank[i]<<"--->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}
}

Atom* Atomcontainer::getatom(int i) const
{
	return atoms[i];
}

float Atomcontainer::getvalue(int i) const
{
	return value[i];
}
int Atomcontainer::getHydrogens(int i) const
{
	return Hydrogens[i];
}

void Atomcontainer::merge(const Atomcontainer & A)
{
	int m=size;
	size=size+A.size;
	
	Adjacency.reserve(size);
	BO.reserve(size);
	NN.reserve(size);
	value.reserve(size);
	rank.reserve(size);
	classRank.reserve(size);
	label.reserve(size);
	Hydrogens.reserve(size);
	H_label.reserve(size);

	for (int i=0;i<A.size;i++)
		{
			atoms.push_back(A.atoms[i]->clone());
			atoms[i+m]->set_nature(A.atoms[i]->get_nature());
		}

	
	Adjacency.insert(Adjacency.end(),A.Adjacency.begin(),A.Adjacency.end());
	BO.insert(BO.end(),A.BO.begin(),A.BO.end());
	NN.insert(NN.end(),A.NN.begin(),A.NN.end());
	value.insert(value.end(),A.value.begin(),A.value.end());
	rank.insert(rank.end(),A.rank.begin(),A.rank.end());
	label.insert(label.end(),A.label.begin(),A.label.end());
	Hydrogens.insert(Hydrogens.end(),A.Hydrogens.begin(),A.Hydrogens.end());
	H_label.insert(H_label.end(),A.H_label.begin(),A.H_label.end());
	
	for (int i=m;i<size;i++)
	{
		for (int j=0;j<Adjacency[i].size();j++)
		{
			Adjacency[i][j]=Adjacency[i][j]+m;
		}
	}
	for (int i=0;i<A.aromatic_atoms.size();i++)
	{
		aromatic_atoms.push_back(A.aromatic_atoms[i]+m);
	}

	for (int i=0;i<size;i++)
	{
		NN[i]=Adjacency[i].size();
	}

}

void Atomcontainer::resetlabel(int m)
{
	for (int i=m;i<size;i++)
	{
		label[i]=0;
	}
}
	

int Atomcontainer::findatomwithlabel(int j) const
{
	int m=-1;
	for (int i=0;i<size;i++)
	{
		if (label[i]==j) 
		{
			m=i;
			break;
		}
	}
	return m;
}
int Atomcontainer::getlabel(int i) const
{
	return label[i];
}
void Atomcontainer::setlabel(int i,int j)
{
	label[i]=j;
}

void Atomcontainer::formbond(int i,int j)
{
	Adjacency[i].push_back(j);
	Adjacency[j].push_back(i);
	BO[i].push_back(1);
	BO[j].push_back(1);
	NN[i]++;
	NN[j]++;
}

void Atomcontainer::breakbond(int i,int j)
{
	for (int k=0;k<Adjacency[i].size();k++)
	{
		if (Adjacency[i][k]==j)
		{
			Adjacency[i].erase(Adjacency[i].begin()+k);
			BO[i].erase(BO[i].begin()+k);
			break;
		}
	}
	for (int k=0;k<Adjacency[j].size();k++)
	{
		if (Adjacency[j][k]==i)
		{
			Adjacency[j].erase(Adjacency[j].begin()+k);
			BO[j].erase(BO[j].begin()+k);
			break;
		}
	}
	NN[i]--;
	NN[j]--;
}

void Atomcontainer::changeBO(int i, int j, int k)
{
	for (int l=0;l<Adjacency[i].size();l++)
	{
		if (Adjacency[i][l]==j)
		{
			BO[i][l]=BO[i][l]+k;
			break;
		}
	}
	for (int l=0;l<Adjacency[j].size();l++)
	{
		if (Adjacency[j][l]==i)
		{
			BO[j][l]=BO[j][l]+k;
		}
	}
}

pair< vector <Atomcontainer>, vector< vector<int> > > Atomcontainer::connectedcomponents() const
{
	vector<Atomcontainer> A;
	vector <int> color;
	queue <int>connected_atoms;
	color.resize(size,-1);
	int counter=0;
	vector< vector<int> > components;
	for (int i=0;i<size;i++)
	{
		if (color[i]==-1)
		{
			components.resize(counter+1);
			color[i]=counter;
			connected_atoms.push(i);
			components[counter].push_back(i);
			while(!connected_atoms.empty())
			{
				int j=connected_atoms.front();
				for (int k=0;k<Adjacency[j].size();k++)
				{
					if (color[Adjacency[j][k]]==-1)
					{
						color[Adjacency[j][k]]=counter;
						connected_atoms.push(Adjacency[j][k]);
						components[counter].push_back(Adjacency[j][k]);
					}
				}
				connected_atoms.pop();
			}
			counter++;
		}
	}
	
	for (int i=0;i<counter;i++)
	{
		A.push_back(form_atomcontainer(components[i]));
	}
	return pair <vector<Atomcontainer>, vector<vector <int> > >(A, components);
}

Atomcontainer Atomcontainer::form_atomcontainer(vector<int> P) const
{
	Atomcontainer A(P.size());
	vector<int>::iterator it;
	

	for (int i=0;i<P.size();i++)
	{
		
		A.atoms.push_back(atoms[P[i]]->clone());
		A.atoms[i]->set_nature(atoms[P[i]]->get_nature());
		for (int j=0;j<Adjacency[P[i]].size();j++)
		{
			int value=Adjacency[P[i]][j];
			it=find(P.begin(),P.end(),value);
			if (it!=P.end())
			{
				
				A.Adjacency[i].push_back(it-P.begin());
			}
			A.BO[i].push_back(BO[P[i]][j]);				
		}
		
	}
	for (int i=0;i<A.getsize();i++)
	{
		A.NN[i]=A.Adjacency[i].size();
		A.Hydrogens[i]=Hydrogens[P[i]];
		/*Will not update label, H_label, value, rank, classrank*/
		
		
	}
	for (int i=0;i<P.size();i++)
	{
		if (isaromatic(P[i]))A.aromatic_atoms.push_back(i);
	}
	
	return A;
}

void Atomcontainer::setatomtypename(int i, std::string S)
{
	atoms[i]->reset_properties();
	atoms[i]->set_atomtype_name(S);
}
void Atomcontainer::setatomsymbol(int i,string S)
{
	atoms[i]->set_atom_symbol(S);
}
void Atomcontainer::setInitialAtomproperties(int i)
{
	atoms[i]->set_initial_properties();
}

void Atomcontainer::setatomvalency(int i)
{
	atoms[i]->set_valency();
}


int Atomcontainer::AromaticBondCount(int i) const
{
	int counter=0;
	for (int k=0;k<BO[i].size();k++)
	{
		if (BO[i][k]==4)counter++;
	}
	return counter;
}

bool Atomcontainer::IsAdjacentAtomWithBO(int n, string S, int i) const
{
	bool answer = false;
	for (int j=0;j<NN[n];j++)
	{
		if (atoms[Adjacency[n][j]]->get_element_name()==S)
		{
			if (BO[n][j]==i)
			{
				answer=true;
				break;
			}
		}
	}
	return answer;
}

int Atomcontainer::getBatomcount(int i) const
{
	int counter = 0;
	for (int j=0;j<Adjacency[i].size();j++)
	{
		if (BO[i][j]!=5)counter++;
	}
	return counter;
}

bool Atomcontainer::HasNBinteractions() const
{
	bool returnvalue = false;
	for (int i=0;i<size;i++)
	{
		if (getBatomcount(i)<getBplusNBatomcount(i))
		{
			returnvalue = true;
			break;
		}
	}
	return returnvalue;
}

vector<int> Atomcontainer::getAromaticAtoms() const
{
	return aromatic_atoms;
}

vector<int> Atomcontainer::getHlabel(int i) const
{
	return H_label[i];
}
void Atomcontainer::setHlabel(int i, std::vector<int> P)
{
	H_label[i]=P;
}

int Atomcontainer::findatomwithHlabel(int i) const
{
	int m =-1;
	for (int j=0;j<size;j++)
	{
		bool found = false;
		for (int k=0;k<H_label[j].size();k++)
		{
			if (H_label[j][k]==i)
			{
				m=j;
				found = true;
				break;
			}
		}
		if (found)break;
	}
	return m;
}
void Atomcontainer::changeHlabel(int i, int j)
{
	for (int k=0;k<size;k++)
	{
		bool found = false;
		for (int l=0;l<H_label[k].size();l++)
		{
			if (H_label[k][l]==i)
			{
				H_label[k][l]=j;
				found =  true;
				break;
			}
		}
		if (found)break;
	}
}

void Atomcontainer::removeHlabel(int i, int j)
{
	for (int k=0;k<H_label[i].size();k++)
	{
		if (H_label[i][k]==j)
		{
			H_label[i].erase(H_label[i].begin()+k);
		}
	}
}

void Atomcontainer::addHlabel(int i, int j)
{
	H_label[i].push_back(j);
}
	
void Atomcontainer::setHydrogens(int i, int j)
{
	Hydrogens[i]=j;
}

int Atomcontainer::getElectronicHashValue(int i) const
{
	int elecHash = 1;
	if (atoms[i]->get_charge()!=0 && atoms[i]->get_up()==0)
	{
		int charge = atoms[i]->get_charge();
		int prime = 1;
		if (charge>0)prime =2;
		else
		{
			prime =3;
			charge =-charge;
		}
		if (atoms[i]->get_atomtype_name().compare("C*")==0)
		{
			prime =5;
		}
		for (int j=0;j<charge;j++)
		{
			elecHash=elecHash*prime;
		}
	}

	if (atoms[i]->get_up()==1)
	{
		if (atoms[i]->get_charge()==1) elecHash=elecHash*11;
		else elecHash=elecHash*7;
	}

	return elecHash;
}

int Atomcontainer::getTotalElectronicValue() const
{
	int totalvalue = 1;
	for (int i=0;i<size;i++)
	{
		totalvalue=totalvalue*getElectronicHashValue(i);
	}
	return totalvalue;
}

void Atomcontainer::addAtom(std::string At, std::string As, int iso, char e, int n)
{
	if (n==0)
	{
		atoms.push_back(new SingleAtom(At,As,iso,e));
	}
	else
		atoms.push_back(new CompositeAtom(At, As));

	setInitialAtomproperties(atoms.size()-1);
	atoms.back()->set_valency();

	size+=1;
	Adjacency.resize(size);//Adjacency list
	BO.resize(size);;//Bond Order
	NN.resize(size);//Number of neighbors
	value.resize(size);//value of each atom (used to evaluate rank)
	rank.resize(size);//rank
	classRank.resize(size);//Class Rank
	Hydrogens.resize(size, 0);//the number of Hydrogens
	label.resize(size, -1);//label of an atom. -1 indicates no labelling
	
	H_label.resize(size);
}

int Atomcontainer::getGroupHash(int i)const
{
	int GroupHash = AtomValueForGA(atoms[i]->get_atomtype_name(), isaromatic(i));	
	
	//cout<<GroupHash;
	for (int j=0;j<getNN(i);j++)
	{
		int adjacentAtom = get_adjacency(i,j);
		int NNValue=AtomValueForGA(getatomtype(adjacentAtom), isaromatic(adjacentAtom));
		int bovalue = get_BO(i,j);
		for (int k=0;k<bovalue;k++)
		{
			//cout<<" (x "<<NNValue<<") "<<getatomtype(adjacentAtom)<<" "<< GroupHash;
			GroupHash=GroupHash*NNValue;
		}
	}
	//cout<<endl;

	return GroupHash;
}
int Atomcontainer::AtomValueForGA(std::string S, bool aromatic) const
{
	int factor =1;
	if (aromatic)factor = 2;
	if (S.compare(0,1,"C")==0)return 2*factor;//returning a factor of 2 will be enough. this will multiply the hash by 2*2^4 = 32
	else if (S.compare(0,1,"H")==0)return 3*factor; //since we have an upper limit on the valency of atom and each connection must be specified - we are fine.
	else if (S.compare(0,1,"O")==0)return 5*factor;
	else if (S.compare(0,1,"N")==0)return 7*factor;
	else if (S.compare(0,1,"S")==0)return 11*factor;
	else if (S.compare(0,1,"P")==0)return 13*factor;
	else 
	{
		int temp = CompositeAtomsRegistry::getIndexOfAtom("{"+ S + "}");
		if (temp!=-1)
			return (prime(6+temp)*factor);
		else return 0;
	}
}
	
int Atomcontainer::getNNElecHash(int i) const
{
	int Hash = 1;
	
	for (int j=0;j<getNN(i);j++)
	{
		Hash=Hash*getElectronicHashValue(get_adjacency(i,j));

		//Hash=Hash*getElectronicHashValue(NNIndex)*(int)pow(11.0,dbcount(NNIndex))*(int)pow(13.0,tpcount(NNIndex));//11 and 13 because getElectronicHashValue goes up to 7
	}
	
	return Hash;
}

/* TO DO: BROKEN
int Atomcontainer::getNNBondCount(int i, int k) const
{
	int counter = 0;

	for (int j=0;j<getNN(i);j++)
	{
		
		
		
		if (get_BO(i,j)!=k)
			counter+=dbcount(get_adjacency(i,j));
	}
	return counter;
}*/

int Atomcontainer::getNNDoubleBondCount(int i) const
{
	int counter = 0;

	for (int j=0;j<getNN(i);j++)
	{
		counter+=dbcount(get_adjacency(i,j));
		if (get_BO(i,j)==2)
			counter--;//not counting any double bonds between atom i and j
	}
	return counter;
}

int Atomcontainer::getNNTripleBondCount(int i) const
{
	int counter = 0;

	for (int j=0;j<getNN(i);j++)
	{
		counter+=tpcount(get_adjacency(i,j));
		if (get_BO(i,j)==3)
			counter--;//not counting any triple bonds between atom i and j
	}
	return counter;
}


	
int Atomcontainer::RingCountOfAtom(int i) const
{
	return Allrings.CountAtom(i);
}
		

map<int,int> Atomcontainer::getClassesFreqMapNeighboringAtom(int i) const
{
	map<int,int> UniqueClassesFreqMap;
	for (int j=0;j<NN[i];j++)
	{
		int classOfNeighbor = classRank[Adjacency[i][j]];
		if (UniqueClassesFreqMap.count(classOfNeighbor)==0)
			UniqueClassesFreqMap.insert(pair<int,int>(classOfNeighbor,1));
		else UniqueClassesFreqMap[classOfNeighbor]++;	
	}

	return UniqueClassesFreqMap;
}

string Atomcontainer::AtomCenteredGroupForGA(int i) const
{
	string S = "";

	S+=getatomtype(i);
	int NNDoubleBonds = getNNDoubleBondCount(i);
	int NNTripleBonds = getNNTripleBondCount(i);
	for (int j=0;j<getNN(i);j++)
	{
		S+="(";
		int bondOrder = get_BO(i,j);
		

		if (bondOrder==2)
		{
			S+="=";
			NNDoubleBonds--;
		}
		if (bondOrder==3)
		{
			S+="#";
			NNTripleBonds--;
		}

		S+=getatomtype(get_adjacency(i,j));

		for (int k=1;k<NNDoubleBonds;k++)
			S+="d";
		for (int k=1;k<NNTripleBonds;k++)
			S+="t";
		S+=")";
	}
	return S;
}

int Atomcontainer::getNNDoubleTripleBondFactors(int i) const
{
	int NNDoubleBonds = getNNDoubleBondCount(i);
	int NNTripleBonds = getNNTripleBondCount(i);
	return (int)(pow(11.0,NNDoubleBonds)*pow(13.0,NNTripleBonds));//11 and 13 because NNElecHash takes prime numbers up to 7
}

set<string> Atomcontainer::GetElements() const 
{
	set<string> elements;
	for (int i =0;i<size;i++)
	{
		elements.insert(atoms[i]->get_element_name());
	}
	return elements;
}

