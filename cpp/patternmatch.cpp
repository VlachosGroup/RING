#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
//#include <map>
#include <set>
#include <vector>
#include <utility>
//using namespace std;
//
//#include <stdio.h>
//#include <stdlib.h>
//
//
//#include "common.h"
#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "singleatom.h"
//#include "compositeatom.h"
//#include "atomcontainer.h"
//#include "molecule.h"
//#include "substructure.h"
#include "patternmatch.h"

//Start of implementation of Patternmatch class
using std::vector; using std::pair; using std::string; using std::set;
using std::cout; using std::endl; 

Patternmatch::Patternmatch()
{
	Molsize= 0;
	Subsize = 0;
	Hfactor.clear();
}
Patternmatch::Patternmatch(const Molecule &Mol,const Substructure &Sub, int patternflag, int FirstAtomIndex)//implementation of constructor
{
	Mo=&Mol;
	Su=&Sub;
	
	Subsize=Sub.size;//setting the size of substructure
	Molsize=Mol.size;//setting the size of the molecule
	for (int i=0;i<Sub.size;i++)//creating the empty rows in M
	{
		M.push_back(vector <int>(Mol.size));	
		
	}
	InitializeAndSetMa(FirstAtomIndex);
	
	if (patternflag==0)
	{
		find_matches(Mol,Sub, patternflag);//finds all the matches
		calcHfactor();
	}
	else find_matches(Mol,Sub, 1);
	
}

Patternmatch::Patternmatch(const Molecule & Mol, const Substructure & Sub, int patternflag)
{
	Mo= &Mol;
	Su= &Sub;
	
	Subsize=Sub.size;//setting the size of substructure
	Molsize=Mol.size;//setting the size of the molecule
	for (int i=0;i<Sub.size;i++)//creating the empty rows in M
	{
		M.push_back(vector <int>(Mol.size));	
		
	}
	InitializeAndSetMa(-1);
	if (patternflag==0)
	{
		find_matches(Mol,Sub, patternflag);//finds all the matches
		calcHfactor();
	}
	else find_matches(Mol,Sub, 1);
	
}




void Patternmatch::InitializeAndSetMa(int FirstAtomIndex)
{
	Ma.resize(Su->size);//creating the firs dimension of size Subsize
	for (int i=0;i<Su->size;i++)
	{
		Ma[i].resize(Su->size);//creating the second dimension of size Subsize

		for (int j=0;j<Su->size;j++)
			Ma[i][j].resize(Mo->size);//creating the third dimension of Ma with size Molsize
	}
	
	// intializing the Ma matrix. the loop checks if the atomtypes match, NN, dbcount and tpcount of the Substructure
	// is less than that of the Molecule
	for (int i=0;i<Su->size;i++)
	{
		string c;
		c=Su->getatomtype(i);
		//cout <<i<<" "<<c<<endl;
		for (int j=0;j<Mo->size;j++)
		{
			
			if (i==0 && FirstAtomIndex>=0)
			{
				if (j==FirstAtomIndex)
				{
					//cout <<"First atom assigned"<<endl;
					Ma[0][i][j]=1;
					break;
				}
			}
			else
			{
				string d;
				d=Mo->getatomtype(j);
				//cout<<d<<endl;
				if (Su->atoms[i]->get_nature()==Mo->atoms[j]->get_nature())//are they the same type of atoms - Single or composite?
				{
					
					if (IsSameCharacteristrics(i,j))// do they have the same characteristics - aromatic, cyclic or allylic? 
					{
						
						if (checkatomtypematch(c,d))//do they ahve the same atomtype? 
						{
							
							if (checkconnectivitymatch(i,j))//do they have compatible connectivity structure - Molecule should have greater than or the same number of Hydrogens, double bonds, or triple bonds?
							{
								
								if(checkatomenvironmentmatch(i,j))//do they have the same environments? 
								{
									Ma[0][i][j]=1;
									
								}
							}
						}
					}
				}
			}
		}			
	}
	
}




bool Patternmatch::checkatomtypematch(string c, string d)
{
	/*we check if strings c and d are of the same atomtype. Note that we check if the nature of the two atoms match. 
	Therefore, if c is a wildcard character ($,&,X), then both atoms are necessarily single atoms because these wildcard 
	characters are assumed to only include single atoms. So we check only the part of the atomtype after the first character to see if they also match in electronic configuration*/
	bool rvalue = false;

	if (c==d) return true;
	else if (c.compare(0,1, "$")==0)
	{
		
		string e = c.substr(1);
		string f = d.substr(1);
		if (e==f) 
		{
			
			return true;
		}
		else
		{
			
			return false;
		}
	}
	else if (c.compare(0,1,"&")==0)
	{
		if ((d.compare(0,1,"N")==0) || (d.compare(0,1,"O")==0)|| (d.compare(0,1,"P")==0)|| (d.compare(0,1,"S")==0))
		{
			string e = c.substr(1);
			string f = d.substr(1);
			if (e==f) return true;
			else return false;
		}
		else return false;
	}
	else if (c.compare(0,1,"X")==0)
	{
		if (d.compare(0,1,"H")!=0)
		{
			string e = c.substr(1);
			string f = d.substr(1);
			if (e==f) return true;
			else return false;
		}
		else return false;
	}
	else return false;
}



bool Patternmatch::checkconnectivitymatch(int i,int j)
{
	
	if( (Su->NN[i]<=Mo->NN[j]) && (Su->Hydrogens[i]<=Mo->Hydrogens[j]) && (Su->dbcount(i)<=Mo->dbcount(j)) && (Su->tpcount(i)<=Mo->tpcount(j)))
		return true;
	else return false;
}


bool Patternmatch::checkatomenvironmentmatch(int i, int j)
{
	int outerflag=1;
	env_set es = Su->AtomEnv[i];
	
	

	for (int k=0;k<es.size();k++)
	{
		env_tuple et = es.get(k);

		int ex_flag=et.get_flag();
		string St=et.get_type();
		
		int check_freq=et.get_freq();
		int comparisonop=et.get_comp_opr();
		int counter=0;
		int checkflag=0;
		if (St=="H" || St=="-H" || St=="~H")
			counter=Mo->Hydrogens[j];
		else if (St=="$" || St=="-$" || St=="~$")
			counter = Mo->Hydrogens[j] + Mo->NN[j];
		else 
		{
			string Groupstring;
			Groupstring = "'$'";
			
			string symbol = Mo->atoms[j]->get_atom_symbol(); 
			if (Mo->isaromatic(j))Groupstring+="[{0}]";
			Groupstring=Groupstring + St;
			Substructure S (Groupstring, patternsize(Groupstring));		
			counter = Patternmatch((*Mo), S, 0, j).GetDistinctMatches();
			//cout<<j<<"  "<<counter<<endl;
			
		}

		
		if (((comparisonop==0) && (counter==check_freq)) || ((comparisonop==1) &&(counter>check_freq)) || ((comparisonop==2) &&(counter>=check_freq)) || ((comparisonop==-1) &&(counter<check_freq)) || ((comparisonop==-2) &&(counter<=check_freq)))
			checkflag=1;
		
		
		if (ex_flag!=checkflag)
		{
			
			outerflag=0;
			break;
		}
	}
	

	if (outerflag==1)
	{
		
		return true;
	}
	else return false;
}

bool Patternmatch::checkatomflagmatch(int i,int j)
{
	int flag=0;
	if (((Su->atom_flag[i][0]==0) && Mo->isaromatic(j)) ||((Su->atom_flag[i][0]==1) && !Mo->isaromatic(j)))
	{
		flag=1;
	}
	else if(((Su->atom_flag[i][1]==0) && Mo->isringatom(j)) ||((Su->atom_flag[i][1]==1) &&!Mo->isringatom(j)))
	{
		flag=1;
	}
	else if (((Su->atom_flag[i][2]==0) && Mo->isallylicatom(j)) ||((Su->atom_flag[i][2]==1) && !Mo->isallylicatom(j)))
	{
		flag=1;
	}
	
	if (Su->atom_flag[i][4]>-1)
	{
		
		int rs =Mo->MaxRingSizeOfAtom(j);
		if (Su->atom_flag[i][3]==-2 && rs < Su->atom_flag[i][4])flag=1;
		else if (Su->atom_flag[i][3]==-1 && rs <= Su->atom_flag[i][4])flag=1;
		else if (Su->atom_flag[i][3]==0 && rs !=Su->atom_flag[i][4])flag=1;
		else if (Su->atom_flag[i][3]==1 && rs >= Su->atom_flag[i][4])flag =1;
		else if (Su->atom_flag[i][3]==2 && rs > Su->atom_flag[i][4])flag =1;
		
	}
	
	if (flag==0)
		return true;
	else return false;
}

bool Patternmatch::IsSameCharacteristrics(int i, int j)
{
	if (Su->atoms[i]->get_nature()==0 && Mo->atoms[j]->get_nature()==0)
	{
		

		bool a = (Su->atoms[i]->IsSymbolAliphatic() && ((Su->atom_flag[i][0]==0) || (Su->atom_flag[i][0]==-1)));
		bool b = Mo->atoms[j]->IsSymbolAliphatic();
		if ((a)|| (!a && !b))
		{
			return checkatomflagmatch(i,j);
		}
		else 
		{
			return false;
		}
	}
	else 
	{

		return checkatomflagmatch(i,j);
	}
}
	


int Patternmatch::check_zeros()//implementation part of the function that checks if any one of rows of M is zero
{
	//the function checks if any of the rows is zero by adding the each of the rows and checking if the sum were zero or not
	//if the sum was zero, the result is 1 indicating there is indeed a zero norm row.
	int result;
	result=0;
	for (int i=0;i<Subsize;i++)
	{
		int adder;
		adder=0;
		for (int j=0;j<Molsize;j++)
		{
			adder=adder+M[i][j];
		}
		if (adder==0)
		{	result=1;
			break;
		}
	}
	return result;
	
}
	
void Patternmatch::find_matches(const Molecule& Mol, const Substructure& Sub, int patterncheckflag)//implementation of the function that finds matches
{
	int* F= new int[Molsize];//creates the array F
	int* H= new int[Subsize];//creates the array H
	//myfile<<"the size of Mol is: "<<Molsize<<endl;
	for (int i=0;i<Molsize;i++)
	{	F[i]=0;}//sets the elements of F to zero
	for (int i=0;i<Subsize;i++)
	{	H[i]=-1;}//sets the elements of H to -1

	int d;// d is the row index
	d=1;
	M=Ma[d-1];//initializing M to Ma[0].
	refine(Mol,Sub);//refinement operation
	
	Ma[d-1]=M;//copying the refined M onto Ma[d-1]
	int k;//k is the column index
	k=0;
	if (check_zeros()==1)//if there are zero norm rows, then the d and k are set so that the while loop below is skipped
	{
		d=1;
		k=Molsize;
	}

	int checking_flag=0;

	while(((d>1)||(k<Molsize)) && checking_flag==0)//ensures that search is within the dimensions of M matrix
	{
		M=Ma[d-1];//setting M as Ma
		
		while(k<Molsize)//if the column pointer is less than the row length (which is equal to the number of columns or Molsize)
		{
			if ((M[d-1][k]==1) && (F[k]==0))//checking if the particular entry in M is one and its corresponding 
			{break;}						//F is zero to ensure potential atom for matching and has not yet been 
			else k++;						//matches with other atoms as yet. If not, check the next column.
		}

		if (k<Molsize)//if k is less than Molsize, indicating the possibility of a match
		{
			Ma[d-1]=M;//copying M onto Ma before refinement
			for (int j=0;j<Molsize;j++)
			{
				if (j!=k) M[d-1][j]=0;//set all other elements in the row to zero
			}
	
			refine(Mol,Sub);//refinement
			
	
			if (check_zeros()==1)//checking if there are zero norm rows. if there are any, then k is incremented and the outer while loop starts over
			{
				k=k+1;
				continue;
			}
			else//if not, indicating a match, then either more atoms need to be matched or the entire pattern has been matched
			{
				if (d<Subsize)//checking if more atoms need to be matched by comparing the value of 'd' at which match occured with Subsize
				{
					
					H[d-1]=k;//updating H, F, d,k, and Ma
					F[k]=1;
		
					d=d+1;
					k=0;
					Ma[d-1]=M;
					continue;
				}
				else
				{
					//Matching has been found
					// the matching atoms and their ranks are appended to Matches and rank_matches vectors
					if (patterncheckflag==1) checking_flag=1;
					H[d-1]=k;
					vector <int> A;
					vector <int> B;
					A.clear();
					B.clear();
					for (int j=0;j<Subsize;j++)
					{
						A.push_back(H[j]);
						
						B.push_back(Mol.getclassRank(H[j]));
					}
					Matches.push_back(A);
					rank_matches.push_back(B);
					F[k]=0;
					k++;
					
				}
			}
		}
		else //if k was greater than Molsize, then all the above steps are skipped and a backtrack is issued
		{
			if (d==1) break;//if d==1, get out of the loop
			else// if not, issue backtracking and update H,F,d,M and k appropriately (that is to set the values to that of the previous match)
			{
				if (H[d-1]>=0)
				{F[H[d-1]]=0;}
				//myfile<<"backtracking at row "<<d<<endl;
				H[d-1]=-1;
				d=d-1;
				
				M=Ma[d-1];
				//print_M();
				k=H[d-1]+1;
				F[H[d-1]]=0;
			
				continue;
			}
		}
	}
	delete[] F;
	delete[] H;

	
}

void Patternmatch::refine(const Molecule& Mol,const Substructure& Sub)//refinement function implementation
{
	//the function basically checks, for all elements of M(i,j) that is 1, if the neighbors of the atom 'i' in question in substructure
	//all have a matching with the neighbors of the atom 'j'in the Molecule.
	
	for (int i=0;i<Subsize;i++)
	{
		for (int j=0;j<Molsize;j++)
		{
			if (M[i][j]==1)
			{
				int counter;
				counter=0;
				
				
				for (int k=0;k<Sub.getNN(i);k++)
				{
					int a1,b1;
					a1=Sub.get_adjacency(i,k);
					b1=Sub.get_BO(i,k);
					vector <int> flag;//this flag vector stores 0 or 1 corresponding to each of the neighbors of the particular atom - this keeps track if a given neighbor has already been matched by a previous fragment atom.
					//this variable keeps a status check to avoid multiple mappings of different fragment atom to the same molecule atom.
					flag.resize(Mol.getNN(j));

					for (int l=0;l<Mol.getNN(j);l++)
					{
						int a2,b2;
						
						a2=Mol.get_adjacency(j,l);
						
						b2=Mol.get_BO(j,l);
						
						if ((M[a1][a2]==1) && (flag[l]==0) && ((b1==b2) || (b1==6) ||((b1==7) && (b2>1 && b2<4))))
						{
							
							if ((Sub.isringbondcheck(i,a1)==-1)||((Sub.isringbondcheck(i,a1)==0) && !Mol.isringbond(j,a2)) || ((Sub.isringbondcheck(i,a1)==1) && Mol.isringbond(j,a2)))
							{
								
								
								
								counter++;//flag counter is increased when the neighbor of the atom i matches with atleast one of the neighbors of j.
								flag[l]=1;
								break; //  
							}
							
						}
					}
					
				}
				if (counter<Sub.getNN(i))//if flag is less than number of nearest neighbors of i, then some atom(s) were not matched 
				{M[i][j]=0;}
			}
		}
	}
}
void Patternmatch::list_matches()// list all the Matches
{
	int n_matches;
	n_matches=(int) Matches.size();//size of the vector Matches
	cout<<"The number of matches is   "<<n_matches<<endl;
	
	for (int i=0;i<n_matches;i++)
	{
		cout<<"the match # "<<i+1<<" is "<<endl;
		for (int j=0;j<Subsize;j++)
		{
			cout<<Matches[i][j]<<"("<<rank_matches[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}
	
}

int Patternmatch::GetDistinctMatches()
{
	/*int counter = 0;
	vector<bool> MatchNotToBeConsidered;
	MatchNotToBeConsidered.resize(rank_matches.size(),false);
	for (int i=0;i<rank_matches.size(); i++)
	{
		if (!MatchNotToBeConsidered[i])
		{
			counter++;
			for (int j=i+1;j<rank_matches.size();j++)
			{
				if (atomsetmatch(i,j))MatchNotToBeConsidered[j]=true;
			}
		}
	}
	return counter;*/

	return GetDistinctAllAndRingMatches().first;
}
			


void Patternmatch::unique_matches()//lists the unique matches
{
	int n_matches;
	n_matches=(int) rank_matches.size();
	int* value=new int[n_matches];
	int* flag=new int[n_matches];//this keeps track if a match has been considered earlier while calculating 
	vector<bool>IsAllRanksSame;
	IsAllRanksSame.resize(n_matches,true);
	for (int i=0;i<n_matches;i++)
	{
		value[i]=1;
		flag[i]=0;
		for (int j=0;j<Subsize;j++)
		{
			value[i]=value[i]*prime(rank_matches[i][j]);//calculates a value for each of the match based on the rank of the atoms of the match
			if (j>0)
				if (rank_matches[i][j-1]!=rank_matches[i][j])IsAllRanksSame[i]=false;
		}	
	}
	
	//by virtue of using Prime numbers and unique rank for calculating the value, each unique match has a unique value
	//the loop below finds out the total number of matches with the same value.
	for (int i=0;i<n_matches;i++)
	{
		if (flag[i]==0)
		{
			unique_matcheslist.push_back(i);
			int counter=1;
			int eq_count=1;
			for (int j=i+1;j<n_matches;j++)
			{
				if (value[j]==value[i])
				{
					flag[j]=1;
					counter++;
					if (atomsetmatch(i,j) && !IsAllRanksSame[i])
					{
						eq_count++;
						unique_matcheslist.push_back(j);
					}
				}	
			}
			counter=counter/(eq_count);
			int k=unique_matches_frequency.size();
			for (int j=k;j<unique_matcheslist.size();j++)
			{
				unique_matches_frequency.push_back(counter);
			}		
			
			
		}
	}
	
	delete[] value;
	delete[] flag;

}

void Patternmatch::print_M()// prints out M
{
	for (int i=0;i<Subsize;i++)
	{
		for (int j=0;j<Molsize;j++)
		{
			cout<<M[i][j]<<"  ";
		}
		cout<<endl;
	}
}

void Patternmatch::print_Ma(int d)// prints out Ma[d]
{
	
	for (int i=0;i<Subsize;i++)
	{
		for (int j=0;j<Molsize;j++)
		{
			cout<<Ma[d][i][j]<<"  ";
		}
		cout<<endl;
	}
}

void Patternmatch::calcHfactor()
{
	Hfactor.clear();
	for (int a=0;a<Matches.size();a++)
	{
		int counter=1;
		
		for (int k=0;k<Subsize;k++)
		{
			if (Su->Hydrogens[k]>0)
			{	
				int Mol_H = (Mo)->Hydrogen_counter(Matches[a][k]).first;
				int Sub_H = Su->Hydrogens[k];
				counter = counter*factorial(Mol_H)/(factorial(Sub_H)*factorial(Mol_H-Sub_H));
				//counter=counter*factorial((Mo)->Hydrogen_counter(Matches[a][k]).first)/factorial((Su)->Hydrogens[k]);
			}
		}
		Hfactor.push_back(counter);
	}
}	


int Patternmatch::H_factor(int a)
{
	return Hfactor[a];
}

int Patternmatch::H_factor_unique_match(int a)
{
	return Hfactor[unique_matcheslist[a]];
}

int Patternmatch::getMatchFrequency(int a)
{
	/*int freq = 0;
	for (int i=0;i<unique_matcheslist.size();i++)
	{
		if (unique_matcheslist[i]==a)
		{
			freq = unique_matches_frequency[i];
			break;
		}
	}
	return freq;*/

	if (a < unique_matcheslist.size())return unique_matches_frequency[a];
	else return 0;
}
bool Patternmatch::atomsetmatch(int i,int j)
{
	int flag=0;
	for (int k=0; k<Subsize;k++)
	{
		for (int l=0;l<Subsize;l++)
		{
			if (Matches[i][k]==Matches[j][l])
			{	flag++;
				break;
			}
		}
	}
	if (flag<Subsize)
		return false;
	else
		return true;
}

int Patternmatch::number_of_matches()
{
	return Matches.size();
}

int Patternmatch::number_of_unique_matches()
{
	return unique_matcheslist.size();
}
vector<int> Patternmatch::get_unique_matches(int i)
{
	return Matches[unique_matcheslist[i]];
}

vector<int> Patternmatch::get_match(int i)
{
	return Matches[i];
}

set<int> Patternmatch::AtomsCoveredBySubstr(int a)
{
	set<int> AlreadyCoveredAtoms;
	for (int i=0;i<Matches.size();i++)
	{
		AlreadyCoveredAtoms.insert(Matches[i][a]);
	}
	return AlreadyCoveredAtoms;
}

bool Patternmatch::IsMatchWithinRing(int p)
{
	if (p<Matches.size())
	{
		for (int i=0;i<Matches[p].size();i++)
		{
			if (!Mo->isringatom(Matches[p][i]))
				return false;
		}
		return true;
	}
	else return false;
}

pair<int,int> Patternmatch::GetDistinctAllAndRingMatches()
{
	int DistinctAll = 0;
	int DistinctRing = 0;

	for (int i=0;i<Matches.size();i++)
	{
		bool isDistinct = true;
		for (int j=0;j<i;j++)
		{
			if (atomsetmatch(i,j))
			{
				isDistinct = false;
				break;
			}
		}
		if (isDistinct)
		{
			DistinctAll++;
			if (IsMatchWithinRing(i))DistinctRing++;
		}
	}
	return pair<int,int> (DistinctAll,DistinctRing);
}


