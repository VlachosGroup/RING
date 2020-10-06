#include <iostream>
#include <fstream>
#include <string>

#include <cstring>
#include <sstream>
#include <map>
#include <vector>
#include <cctype>
#include <algorithm>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Classheader.h"
#include "AdditionalFunctions.h"


//start of Molecule implementation





Molecule::Molecule(string stringname, int m)
		:Atomcontainer(m)// Constructor of Substructure class
{
	smilesstring=stringname;//assigning smilesstring as stringname
	int r;
	r=-1;//atom counters assignment
	Triplet* ring=new Triplet[100];//ring Triplet declaration
	
	int ringcounter;
	ringcounter=0;//ring counter initialization
	int i;
	i=0;
	string atomtypename;
	atomtypename="";//string that will store the value of the atomtype_name 
	string atomsymbol;//string that will store the value of the atomtype_symbol
	atomsymbol="";

	int parent;
	parent=-1;//variable that has the information of parent of the new atom to be added
	int bond_ident;
	bond_ident=1;//gives bond_identification. 1 is single, 2 is double
	int* branch=new int[m];//initialization of array that stores branch nodes
	int branchcounter;//counts the number of branches. Increases by 1 when a new branch forms, decreases when it ends
	branchcounter=0;//initializing the branchcounter to zero
	int isotopevalue=0;//0 is normal, else stores the actual isotope atomic mass.

	vector < string > ringnumber;//stores as string the number index of hte ring in question.
	bool IsCompositeAtomFound=false;
	bool InsideSquareBrackets=false;
	bool ExplicitSingleBond=false;
	bool ExplicitDoubleBond=false;
	int NumberOfRings=0;//keeps track of the number of Rings the atom is involved in. 
	bool IsAfterPercentageSign=false;
	bool IsBeforeElementinSquareBraces = true;
	while(i<smilesstring.length())//each character of the smilesstring is read and processed one by one
	{
		string c;
		
	
		c=smilesstring[i];// c stores the character in position i.
	
		
		if (c=="[")
		{
			//cout<<InsideSquareBrackets<<endl;
			//cout<<"got into square brackets"<<endl;
			InsideSquareBrackets = true;
			//cout<<InsideSquareBrackets<<endl;
			atomsymbol+="[";
			IsBeforeElementinSquareBraces = true;
			i++;
			c=smilesstring[i];
		}
		

		
		if (InsideSquareBrackets && atomsymbol.length()>=2 && smilesstring.compare(i,1,"H")==0 && !isalpha(smilesstring[i-1]) && smilesstring.compare(i-1,1,"}")!=0)
		{
			//this if statement is for cases like [n+H] or [{HTA}+H], etc
			atomsymbol+=c;
			i++;
			Hydrogens[r]=1;
			cout<<smilesstring<<endl;
			
		}
		c=smilesstring[i];
		if (c=="]")
		{
			InsideSquareBrackets = false;
			
			atomsymbol+="]";
		}
		
		if(isalpha(smilesstring[i]) || (c=="{") ) //checking if the character is an alphabet
		{
					
			NumberOfRings=0;
			IsAfterPercentageSign=false;
			if (c=="{")
			{
				IsCompositeAtomFound=true;
				atomsymbol+="{";
				i++;
				c=smilesstring[i];
				while (c!="}")
				{
					atomtypename+=smilesstring[i];
					atomsymbol+=smilesstring[i];
					i++;
					c=smilesstring[i];
				}
				//cout<<"i th character is "<<smilesstring[i]<<endl;
				atomsymbol+="}";
			}

			if (isalpha(smilesstring[i]))
			{			
				atomsymbol+=c;			
				if (c == "D")//if char is D, Deuterium, then atomtype still is H. and isotope flag is set to 1.
				{
					atomtypename="H";
					isotopevalue=2;
				}
				else 				
				{	
					atomtypename=toupper(smilesstring[i]);
				}
				
			}
			
			r++;//incrementing atomcounter
			if (IsBeforeElementinSquareBraces) IsBeforeElementinSquareBraces = false;
			if (InsideSquareBrackets && smilesstring.compare(i+1,1,"H")==0)
			
			{//todo - make an explicit character comparison!
				
				atomsymbol+=smilesstring[i+1];
				i++;
				Hydrogens[r]=1;//you can ignore the character and add 1 Hydrogen to Hydrogens because there
				//is atleast one Hydrogen attached to the atom corresponding to r. Example [nH] has one aromatic
				//nitrogen attached to 1 Hydrogen while [CH3] is a methyl group with 3 Hydrogens
				if (smilesstring.compare("[{M}]")==0)cout<<"here "<<smilesstring<<endl;
			}
			


				
			if (parent>=0)//if the parent value is non-negative (implying the atom is not first), 
				//then update adjacency list of both the daughter and parent atoms
			{
					//cout<<"adding atoms "<<r<<endl;
					//cout<<"the size of atomcontainer is "<<Adjacency.size()<<endl;
					//Adjacency[r].position[NN[r]]=parent;
				
				Adjacency[r].push_back(parent);
				//cout<<"reached1"<<endl;
				//BO[r].position[NN[r]]=bond_ident;
				BO[r].push_back(bond_ident);
				//cout<<"reached2"<<endl;
				//Adjacency[parent].position[NN[parent]]=r;
				Adjacency[parent].push_back(r);
				//cout<<"reached3"<<endl;
				//BO[parent].position[NN[parent]]=bond_ident;
				BO[parent].push_back(bond_ident);
				//cout<<"reached"<<endl;
									
				bond_ident=1;
				string parentsymbol;
				parentsymbol=atoms[parent]->get_atom_symbol();
				if ((islower(atomsymbol[0])|| (atomsymbol.length()>1 && islower(atomsymbol[1]))) && islower(parentsymbol[0]))
				{
					if (!ExplicitSingleBond && !ExplicitDoubleBond)
					{
						BO[r][BO[r].size()-1]=4;
						BO[parent][BO[parent].size()-1]=4;
					}
					ExplicitSingleBond = false; 
					ExplicitDoubleBond = false;
				}
					
			}
			parent=r;//now for subsequent atom addition, r is the parent
			
		}
		if (c=="%")
		{
			NumberOfRings++;
			ringnumber.push_back("");
			ringcounter++;
			//cout<<"adding ring counter here"<<endl;
			IsAfterPercentageSign=true;
		}
		//c=smilesstring[i];
		if (isdigit(smilesstring[i]))//checking if the character is a digit
		{
			//cout<<"is a digit"<<endl;
			//cout<<"i is "<<i<<" "<<smilesstring<<"  "<<smilesstring[i]<<endl;

			//cout<<InsideSquareBrackets<<endl;
			if (InsideSquareBrackets)
			{
				//if inside square brackts, they refer to Hydrogens or to isotopes
				if (!IsBeforeElementinSquareBraces)
				{
					Hydrogens[r]=atoi(c.c_str());
					atomsymbol+=smilesstring[i];
					
				}
				else
				{
					isotopevalue = isotopevalue*10 + atoi(c.c_str());
				}
				
				
				
			}
			else
			{
				if (!IsAfterPercentageSign)
				{
					//cout<<"adding ring counter "<<endl;
					ringcounter++;
					NumberOfRings++;
					ringnumber.push_back("");
				}
				ringnumber[NumberOfRings-1]+=smilesstring[i];//add this character to ringnumber
				if (bond_ident==2)ringnumber[NumberOfRings-1]+="=";
				if (bond_ident==3)ringnumber[NumberOfRings-1]+="#";
				if (bond_ident==5)ringnumber[NumberOfRings-1]+="_";
				bond_ident = 1;
				
			}				
		}

		

		if (c=="(")//if open brackets are identified, then increment branch counter and update branch array with the parent atom
		{
			branchcounter++;
			branch[branchcounter-1]=parent;
		
		}
		if (c==")")//if the closed brackets are identified, decrease branch counter and assign the branch point as parent for subsequent addition
		{
			parent=branch[branchcounter-1];
			branchcounter--;
		
		}
		if ((c=="+")||(c==":")||(c=="*")||(c=="."))//checking atomtypes and updating the atomtypename and atomsymbol appropriately
		{
			atomsymbol=atomsymbol+smilesstring[i];
			atomtypename=atomtypename+smilesstring[i];
		
		}

		if (c=="-")
		{
			if (InsideSquareBrackets)
			{
				atomsymbol+=smilesstring[i];
				atomtypename+=smilesstring[i];
			}
			else
			{
				bond_ident=1;
				ExplicitSingleBond=true;
			}
					
		}

		if (c=="=")//identifying double bonds
		{
			bond_ident=2;
			ExplicitDoubleBond=true;
	
		}
		if (c=="#")//identifying triple bonds
		{
			bond_ident=3;
		
		}
		if (c=="_")//identifying nonbonded interactions
		{
			bond_ident=5;
		}
		/*if (c=="'")//identifying isotopes
		{
			isotopenumber=1;
			atomsymbol=atomsymbol+c;
			//cout<<"isotope found"<<endl;
		
		}*/
		/*when the next character is either one of the element characters or if the current position is end of string
		  it is understood that new atom needs to be initialized with the current values. The following codes ensure that*/
				
		string c2;
		if (i<smilesstring.length()-1)
			c2= smilesstring[i+1];

		if ((((i<smilesstring.length()-1 && isalpha(smilesstring[i+1]))|| (c2=="[") || (c2=="{"))&& !InsideSquareBrackets) ||(i==(smilesstring.length()-1)))
		{
			if(r>=0)//r is initially set to -1
			{
				if (IsCompositeAtomFound)
				{
					
					atoms.push_back(new CompositeAtom(atomtypename, atomsymbol));
				}
				else
				{
					char elementname=atomtypename[0];
					
					atoms.push_back(new SingleAtom(atomtypename,atomsymbol,isotopevalue,elementname));//instantiating an object of Atom class
					//cout<<"the atom type name is "<<atomtypename<<endl;
					//cout<<"get atom symbol gives "<<atoms[r]->get_atom_symbol()<<endl;

					
					isotopevalue=0;
					
				}
				atomtypename="";
				atomsymbol="";
				//if (r==0)cout<<"1. the first atom is "<<atoms[0]->get_atomtype_name()<<"   "<<atoms[0]->get_up()<<endl;
				setInitialAtomproperties(r);//setting properties of the atomtype
				//if (r==0)cout<<"2. the first atom is "<<atoms[0]->get_atomtype_name()<<"   "<<atoms[0]->get_up()<<endl;
				atoms[r]->set_valency();//setting valency of the atom (atomtype)
				//if (r==0)cout<<"3. the first atom is "<<atoms[0]->get_atomtype_name()<<"   "<<atoms[0]->get_up()<<endl;
				IsCompositeAtomFound=false;
				
				int l;
				
				for (int k=0;k<ringnumber.size();k++)
				{
					string ringvalue="";
					int bo=1;
					for (int rv=0;rv<ringnumber[k].length();rv++)
					{
						
						if (isdigit(ringnumber[k][rv]))ringvalue+=ringnumber[k][rv];
						if (ispunct(ringnumber[k][rv]))
						{
							string bovalue;
							bovalue= ringnumber[k][rv];
							if (bovalue=="=")bo=2;
							if (bovalue=="#")bo=3;
							if (bovalue=="_")bo=5;
						}
					}
					l=atoi(ringvalue.c_str());//conversion of ringnumber value from string to integer
					
					ring[l-1].second=ring[l-1].first;
					ring[l-1].first=r;//updating the ring pair
					ring[l-1].third=bo;
					
					
				}

	//			cout<<r<<endl;
			}
			ringnumber.clear();
			
		}

		i++;
	}

	//cout<<"before "<<atoms[0]->get_atomtype_name()<<"  "<<atoms[0]->get_up()<<"  "<<atoms[0]->get_valency()<<" "<<atoms[0]->get_charge()<<endl;
	
	//ring bonds need to be added and are done below. Note that ringcounter was incremented each time a number was detected
	//For each ring, ringcounter doublecounts. Hence the need to divide the final value by 2.
	ringcounter=ringcounter/2;

	for (int j=0;j<ringcounter;j++)
	{
		int j2;
		j2=ring[j].first;
		int j3;
		j3=ring[j].second;
		if (!atoms[j2]->IsSymbolAliphatic() && !atoms[j3]->IsSymbolAliphatic())ring[j].third=4;
		Adjacency[j2].push_back(j3);
		BO[j2].push_back(ring[j].third);

		Adjacency[j3].push_back(j2);
		BO[j3].push_back(ring[j].third);
		/*
		Adjacency[j2].position[NN[j2]]=j3;
		BO[j2].position[NN[j2]]=1;

		Adjacency[j3].position[NN[j3]]=j2;
		BO[j3].position[NN[j3]]=1;
		*/
		
			
	}


	for (int i=0;i<size;i++)
	{
		NN[i]=Adjacency[i].size();
	}
	
	remove_Hydrogens();
	aromatic_atoms.clear();
	aromatic_rings.clear();
	update_aromaticity_details();
	
	/*cout<<"the initial adjacency matrix is"<<endl;
	for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		cout<<i<<" ("<<NN[i]<<")  "<<atoms[i]->get_atom_symbol()<<" ("<<Hydrogens[i]<<")"<<" ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}*/
	delete [] ring;
	delete [] branch;
	
	Readjustproperties();//Readjusts valency to account for Higher oxidation states and Hydrogens
	calculateHydrogens();
	find_all_rings();
	Allrings.sizesort();
	find_aromatic_rings();
	update_aromaticity_details();
	//cout<<"going to update square braces"<<endl;
	UpdateSquareBraces();
	//cout<<"update done"<<endl;
	EvaluateMF();
	find_allylic_atoms();	
	/*cout<<"the initial adjacency matrix is"<<endl;
	for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		cout<<i<<" ("<<NN[i]<<")  "<<atoms[i]->get_atom_symbol()<<" ("<<Hydrogens[i]<<")"<<" ("<<atoms[i]->get_up()<<")"<<" ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}*/
	//cout<<atoms[0]->get_atomtype_name()<<"  "<<atoms[0]->get_charge()<<"  "<<atoms[0]->get_up()<<endl;
	//cout<<"done"<<endl;
	
	
}


Molecule::Molecule(Atomcontainer& AtC)
{
	smilesstring="";
	size=AtC.getsize();
	Adjacency.resize(size);
	BO.resize(size);
	NN.resize(size);
	value.resize(size);
	Hydrogens.resize(size);
	label.resize(size);
	rank.resize(size);
	classRank.resize(size);
	aromatic_atoms.resize(size);
	H_label.resize(size);
	atoms.clear();

	
	for (int i=0;i<size;i++)
	{
		//Atom* at = new Atom();
		//at = AtC.getatom(i);
		//atoms.push_back(at);
		//atoms.push_back(AtC.atoms[i]->clone());
		atoms.push_back(AtC.getatom(i)->clone());
		atoms[i]->set_nature(AtC.getatom(i)->get_nature());
		NN[i]=AtC.getNN(i);
		//value[i]=AtC.getvalue(i);
		label[i]=AtC.getlabel(i);
		rank[i]=AtC.getrank(i);
		classRank[i]=AtC.getclassRank(i);
		Hydrogens[i]=AtC.getHydrogens(i);
		
		H_label[i] = AtC.getHlabel(i);

		
		for (int j=0;j<NN[i];j++)
		{
			Adjacency[i].push_back(AtC.get_adjacency(i,j));
			BO[i].push_back(AtC.get_BO(i,j));
		}
		//delete at;
		 
	}
	aromatic_atoms=AtC.getAromaticAtoms();
	
	
	for (int i=0;i<size;i++)
	{
		setInitialAtomproperties(i);
		atoms[i]->set_valency();
		atoms[i]->DropRingIdentifier();

		if ((atoms[i]->get_atom_symbol()=="C=") || (atoms[i]->get_atom_symbol()=="O="))
			cout<<"i need to fix this"<<"  "<<atoms[i]->get_atom_symbol()<<endl;

		int aromatic_bond_count=0;
		
		for (int j=0;j<NN[i];j++)
		{
			if (BO[i][j]==4)
				aromatic_bond_count++;
		}
		if (aromatic_bond_count==0 || aromatic_bond_count ==1)//if the atom is not aromatic, then check if the atom was initially aromatic and if yes, update atomsybol!
		{	
			atoms[i]->MakeSymbolAliphatic();
			
		}
		
				
	}
	remove_Hydrogens();
	find_all_rings();
	Allrings.sizesort();

	vector<int>BrokenAromaticAtomsList;
	

	for (int i=0;i<aromatic_atoms.size();i++)
	{
		if (!isringatom(aromatic_atoms[i]))atoms[aromatic_atoms[i]]->MakeSymbolAliphatic();
		
		if (isringatom(aromatic_atoms[i]) && atoms[aromatic_atoms[i]]->IsSymbolAliphatic())
		{
			BrokenAromaticAtomsList.push_back(aromatic_atoms[i]);
		}
		
	}
	
	for (int i=0;i<Allrings.size();i++)
	{
		
		bool HasBrokenAtom = false;
		for (int j=0;j<Allrings.get(i).size();j++)
		{
			
			for (int k=0;k<BrokenAromaticAtomsList.size();k++)
			{
				if (Allrings.get(i).get(j)==BrokenAromaticAtomsList[k])
				{
					HasBrokenAtom = true;
					break;
				}
			}
			if (HasBrokenAtom)break;
		}
		if (HasBrokenAtom)PerceiveKekule(Allrings.get(i));
	}

	/*cout<<"the initial adjacency matrix is"<<endl;
	for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		cout<<i<<" ("<<NN[i]<<")  "<<atoms[i]->get_atom_symbol()<<" ("<<Hydrogens[i]<<")"<<" ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}*/
	

	aromatic_atoms.clear();
	aromatic_rings.clear();
	update_aromaticity_details();
	Readjustproperties();//Readjusts valency to account for Higher oxidation states and Hydrogens
	//calculateHydrogens();
	
	find_aromatic_rings();
	update_aromaticity_details();
	
	UpdateSquareBraces();
	EvaluateMF();
	//cout<<"going to do allylic atoms"<<endl;
	find_allylic_atoms();

	/*for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		cout<<i<<" ("<<NN[i]<<")  "<<atoms[i]->get_atom_symbol()<<" ("<<Hydrogens[i]<<")"<<" ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}*/
	
}


void Molecule::find_all_rings()
{
	//this function finds the set of all rings in a molecule
	//cout<<"size of atomcontainer is"<<size<<endl;
	vector <vector <int> > container;
	vector <int> row;
	
	vector <Path> all_paths;
	for (int i=0;i<size;i++)
	{
		row.clear();
		row.push_back(i);
		for (int j=0;j<NN[i];j++)
		{
			row.push_back(Adjacency[i][j]);
		}
		container.push_back(row);
		
	}
	//cout<<"container size is"<<container.size()<<endl;
	//check for leaves
	int leaves_counter;
	leaves_counter=0;
	for (int i=0;i<size;i++)
	{
		if (NN[i]==1) 
		{
			leaves_counter=1;
			break;
		}
	}
	//cout<<leaves_counter<<endl;
	vector <int> remove_atoms;
	while (leaves_counter!=0)
	{
		leaves_counter=0;
		
		for (int i=0;i<container.size();i++)
		{
			if (container[i].size()==2)
			{
				int j=container[i][1];
				for (int k=0;k<container[j].size()-1;k++)
				{
					if (container[j][k+1]==i)
					{
						container[j].erase(container[j].begin()+k+1);
						container[i].erase(container[i].begin()+1);
						//cout<<"here"<<endl;
						break;		
					}
											
				}
				remove_atoms.push_back(i);	
				if (container[j].size()<2)
				{
					remove_atoms.push_back(j);
				}
				
			}
		}
		
		//cout<<"I can reach"<<endl;
		
		//cout<<" reach"<<endl;
		for (int i=0;i<container.size();i++)
		{
			if (container[i].size()==2)
			{
				leaves_counter=1;
				break;
			}
		}
		//cout<<"REACH"<<endl;
	}
	for (int i=remove_atoms.size()-1;i>=0;i--)
	{
		for (int j=0;j<container.size();j++)
		{
			if (container[j][0]==remove_atoms[i])
			{
				container.erase(container.begin()+j);
				break;
			}
		}
	}

	remove_atoms.clear();
	//cout<<"the number of atoms is"<<container.size()<<endl;
	/*for (int i=0;i<container.size();i++)
	{
		cout<<container[i][0]<<"--> ";
		
		cout<<endl;
	}*/
	for (int i=0;i<container.size();i++)
	{
		Path A;
		for (int j=0;j<container[i].size()-1;j++)
		{
			if (container[i][j+1]>container[i][0])
			{
				A.clear();
				A.add(container[i][0]);
				A.add(container[i][j+1]);
				all_paths.push_back(A);
			}
		}
	}
	//cout<<"size of all paths is: "<<all_paths.size()<<endl;
	/*for (int i=0;i<all_paths.size();i++)
	{
		all_paths[i].print();
	}*/


	
	while (container.size()>0)
	{
		int Min_NN=100;
		int d_atom;
		for (int i=0;i<container.size();i++)
		{
			if (container[i].size()-1<=Min_NN)
			{
				d_atom=i;
				Min_NN=container[i].size();
			}
		}
		int del_atom=container[d_atom][0];
		
		
		vector <int> remove_paths;
		
		for (int i=0;i<all_paths.size();i++)
		{
			if ((all_paths[i].getfirst()==del_atom) && (all_paths[i].getlast()==del_atom))
			{
				
				all_paths[i].remove(0);
				Allrings.add(all_paths[i]);
				remove_paths.push_back(i);
				
			}
			else if ((all_paths[i].getfirst()==del_atom) || (all_paths[i].getlast()==del_atom))
			{
				for (int j=i+1;j<all_paths.size();j++)
				{
					if ((all_paths[j].getfirst()==del_atom) || (all_paths[j].getlast()==del_atom))
					{
						
						Path P1=all_paths[i];
						Path P2=all_paths[j];
						if (P1.intersection(P2)==0)
						{
							P1.join(P2,del_atom);
							all_paths.push_back(P1);
						}

					}
				}
				remove_paths.push_back(i);

			}
			

		}
		
		for (int i=remove_paths.size()-1;i>=0;i--)
		{
			all_paths.erase(all_paths.begin()+remove_paths[i]);
		}
	
			
		container.erase(container.begin()+d_atom);
		
		//cout<<"done"<<endl;
		//cout<<all_paths.size()<<endl;
		/*for (int i=0;i<all_paths.size();i++)
		{
			all_paths[i].print();
		}*/
		//cout<<"done"<<endl;
	}
}


void Molecule::find_aromatic_rings()
{
	if (Allrings.size()>0)
		Allrings.sizesort();
	Path P;
	//cout<<"number of rings is "<<Allrings.size()<<endl;
	for (int i=0;i<Allrings.size();i++)
	{
		P=Allrings.get(i);
		//P.print();
		int flag=0;
		int non_sp2_atomcount=0;
		for (int j=0;j<P.size();j++)
		{
			if(isaromatic(P.get(j)))
			{
				flag++;
				//cout<<"yes!"<<endl;
			}
		}
		if (flag==P.size())
		{
			//PerceiveKekule(P);
			aromatic_rings.push_back(i);
		}
		else
		{
			int pi_elec;
			pi_elec=0;
			//cout<<"size of P is "<<P.size()<<endl;
			//P.print();
			if (!RingWithNBInteractions(P))
			{
				for (int j=0;j<P.size();j++)
				{
					string c;
					//cout <<j<<endl;
					int Catom=P.get(j);				
				
					c=atoms[Catom]->get_atomtype_name();
					
					/*
					There are six steps in aromaticity check
					1. First check for composite or single atom - composite atom are assumed to be non aromatic
					2. Check for carbon and its cases 
					3. Check for heteroatom and its cases
					4. Check for the case of more than one double bond - not aromatic
					5. Check for quartenary and higher atom - not aromatic
					6. Check for mono radicals - uncharged ones without a double bond adds one pi electron				
					*/
					int db_in_ring = AdjacentBondInRingWithOrder(P, Catom, 2);
					

					if (atoms[Catom]->get_nature()==1)
					{
						non_sp2_atomcount++;
					}
					else//if the atom is SingleAtom
					{
						
					
						if (c.compare(0,1,"C")==0)//if the atom is a carbon
						{
							if (c=="C")
							{
								
								if (!atoms[Catom]->IsSymbolAliphatic())pi_elec++;//if it is neutral aromatic carbon, add one to pi electron count
								if (dbcount(Catom)==0)non_sp2_atomcount++;//if the neutral carbon has no double bonds, its not sp2!
								if (db_in_ring==1)pi_elec++;//if the neutral carbon has one double bond in the ring, then add one to pi electron count
								else 
								{
									if (db_in_ring==0 &&dbcount(Catom)==1 && IsAdjacentAtomWithBO(Catom, "C", 2))non_sp2_atomcount++;//NOTE-not really non sp2 in this case, but non aromatic nevertheless
								}
							}
							if (c=="C-")
							{
								if (!atoms[Catom]->IsSymbolAliphatic())pi_elec+=2;//if it is an aromatic carbanion, add 2 electrons
								if (dbcount(Catom)==0)pi_elec+=2;//if it is a carbanion with no double bonds, then again add 2 electrons corresponding to the lone pair						
							}
						}					
						else//if not a carbon
						{
							if (atoms[Catom]->get_lp()>=1)
							{
								if (c=="O")pi_elec+=2;
								if (c=="S")  
								{
									if (!atoms[Catom]->IsSymbolAliphatic())pi_elec+=2;
									else
									{
										if (IsAdjacentAtomWithBO(Catom, "O", 2))non_sp2_atomcount++;
										else pi_elec+=2;
									}
								}
								if ((c=="N") || (c=="P"))
								{
									if (isaromatic(Catom))
										pi_elec+=(Hydrogens[Catom]+1);
									else
									{
										if (db_in_ring==1)pi_elec++;
										if (db_in_ring==0)pi_elec+=2;
										if (Hydrogens[Catom]>1)non_sp2_atomcount++;
									}									
								}
							}				
	
							if (((c=="N+") || (c=="P+")))
							{
								if (isaromatic(Catom))pi_elec++;
								else
								{
									if (db_in_ring==1)
										pi_elec++;
									else non_sp2_atomcount++;
								}
							}
							if ((c=="S+"))
							{	
								if (isaromatic(Catom))pi_elec++;
								else
								{
									if ((db_in_ring==1) && !(IsAdjacentAtomWithBO(Catom, "O", 2)))pi_elec++;
									else non_sp2_atomcount++;
								}
							}
						}
					}
					if (db_in_ring==2)non_sp2_atomcount++;
					if ((getBatomcount(Catom)+ Hydrogens[Catom])>=4)non_sp2_atomcount++;
					if ((atoms[Catom]->get_up()==1)&& (atoms[Catom]->get_charge()==0) && (dbcount(Catom)==0))pi_elec++;
					
				}
				//cout<<"the number of pi electrons is "<<pi_elec<<endl;
				if (((pi_elec%4)==2) && (non_sp2_atomcount==0)&& (pi_elec>2))
				{
					aromatic_rings.push_back(i);
					for (int j=0;j<P.size();j++)
					{
						if (!(isaromatic(P.get(j))))aromatic_atoms.push_back(P.get(j));					
					}	
					
				}				
			}			
		}		
	}
}


void Molecule::PerceiveKekule(Path P)
{
	vector<int> PiElectronCount;
	vector<int> BondOrder;
	PiElectronCount.resize(P.size());
	BondOrder.resize(P.size(), -1);
	int totalcount=0;
	for (int j=0;j<P.size();j++)
	{
		string c;
		atoms[P.get(j)]->MakeSymbolAliphatic();
					
		int Catom=P.get(j);				
				
		c=atoms[Catom]->get_atomtype_name();
		if (AdjacentBondInRingWithOrder(P, Catom, 1)==0)
		{
		
			if (c.compare("C")==0)
			{
				if (dbcount(Catom)==0)PiElectronCount[j]=1;
				else PiElectronCount[j]=0;
			}
			else if ((c.compare("C-")==0) || (c.compare("O")==0) || (c.compare("S")==0))PiElectronCount[j]=2;
			else if ((atoms[Catom]->get_up()==1)&& (atoms[Catom]->get_charge()==0) && (dbcount(Catom)==0))PiElectronCount[j]=1;
			else if ((c.compare("N")==0) || (c.compare("P")==0))
			{
				PiElectronCount[j] = Hydrogens[Catom]+1;
			}
			else if ((c.compare("N+")==0) || (c.compare("P+")==0)||(c.compare("O+")==0) || (c.compare("S+")==0))PiElectronCount[j]=1;
			else PiElectronCount[j]=0;
			totalcount+=PiElectronCount[j];
		}
	}
	
   // P.print();
	bool fixedBO = false;
	for (int i=0;i<P.size();i++)
	{
		/*if ((PiElectronCount[i]==0) || (PiElectronCount[i]==2) || (atoms[P.get(i)]->get_up()==1))
		{
			if (BondOrder[i]==-1)BondOrder[i]=1;
			BondOrder[(i+P.size()-1)%P.size()]=1;
			fixedBO = true;
			
		}*/
		int totalbonds = getBatomcount(P.get(i))+dbcount(P.get(i))+2*tpcount(P.get(i))+Hydrogens[P.get(i)];
		if (atoms[P.get(i)]->get_valency()==totalbonds ||(atoms[P.get(i)]->get_valency()==totalbonds-2) ||(atoms[P.get(i)]->get_valency()==totalbonds-4))
		{
			if (BondOrder[i]==-1)
			{
				BondOrder[i]=1;
				BondOrder[(i+P.size()-1)%P.size()]=1;
				fixedBO = true;
			}
		}
		int bo = find_BO(P.get(i),P.get((i+1)%P.size()));
		if (bo!=-1 && bo!=4)BondOrder[i]=bo;

	}

	int j=P.size();
	if (fixedBO)
	{
		for (int i=0;i<P.size();i++)
		{
			if (BondOrder[i]>=1)
			{
				j=i+P.size()+1;
				break;
			}
		}
	}
	for (int k=j;k<j+P.size();k++)
	{
		if (BondOrder[k%P.size()]==-1)
		{
			if (BondOrder[(k+1)%P.size()]==2 || BondOrder[(k-1)%P.size()]==2)BondOrder[k%P.size()]=1;
			else BondOrder[k%P.size()]=2;
		}
	}
	
	//for (int i=0;i<BondOrder.size();i++)
	//	cout<<BondOrder[i]<<"  ";
	//cout<<endl;

	for (int i=0;i<BondOrder.size();i++)
	{
		//cout<<P.get(i)<<"  "<<P.get((i+1)%P.size())<<"  "<<BondOrder[i]<<endl;
		if (isaromaticbond(P.get(i), P.get((i+1)%P.size())))
			setBO(P.get(i),P.get((i+1)%P.size()),BondOrder[i]);
	}

	/*cout<<"the initial adjacency matrix is"<<endl;
	for (int i=0;i<size;i++)//printing out the adjacency list 
	{
		cout<<i<<" ("<<NN[i]<<")  "<<atoms[i]->get_atom_symbol()<<" ("<<Hydrogens[i]<<")"<<" ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"  ";
		}
		cout<<endl;
	}*/
	
}

bool Molecule::isallylicatom(int n) const
{
	/*set<int>::iterator it;
	it = allylic_atoms.find(n);
	if (it!=allylic_atoms.end())return true;
	else return false;*/
	if (allylic_atoms.count(n)==0)return false;
	else return true;
}


void Molecule::update_aromaticity_details()
{
	for (int i=0;i<size;i++)
	{
		
		string c;
		
		c=atoms[i]->get_atom_symbol();
		
		if (isaromatic(i) || (!atoms[i]->IsSymbolAliphatic())) 
		{
			if (!(isaromatic(i)))
			{
				
				aromatic_atoms.push_back(i);
			}
			else atoms[i]->MakeSymbolAromatic();
			
			for (int j=0;j<NN[i];j++)
			{		
				//cout<<"the adjacent atom is "<<Adjacency[i][j]<<" "<<d<<"  "<<isaromatic(Adjacency[i][j])<<" "<<isaromatic(i)<<"   "<<isaromatic(2)<<endl;
				if (isaromatic(Adjacency[i][j]) || (!atoms[Adjacency[i][j]]->IsSymbolAliphatic()))
				{
					
					if (InSameRing(i, Adjacency[i][j]))BO[i][j]=4;
					if (!(isaromatic(Adjacency[i][j])))
					{
						aromatic_atoms.push_back(Adjacency[i][j]);
					}
					else atoms[Adjacency[i][j]]->MakeSymbolAromatic();

				}
			}
		}
	}
}


string Molecule::unique_smiles(int z)
{

	
	
	EvaluateInitialValue(z);	//Evaluate the initial value of all the atoms

	
	evaluate_rank();//evaluate rank after the initial value has been calculated

	for (int r =0;r<size;r++)
		cout<<rank[r]<<"  "<<value[r]<<endl;
	cout<<"after initial evaluation"<<endl;
	print_adjacency_list();
	findRanksOfAllAtoms();
	cout<<"after finding final ranks"<<endl;
	classRank=rank;

	print_adjacency_list();
	
	while(distinctRankCount()!=size)//if there are ties, then unique rank count is less than size of hte molecules)
	{
		//Now we should break ties
		cout<<"going to break ties"<<endl;
		breakTies();		
		cout<<"finished ties"<<endl;
		print_adjacency_list();
		findRanksOfAllAtoms();
		cout<<"After calculating ranks"<<endl;
		print_adjacency_list();

	}
	
	
	sort_adjacency();//sorting the adjacency list in order of ranks and double bonds
	
	
	cout<<"the sorted adjacency list"<<endl;
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<classRank[i]<<"  "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}
	//cout<<"going to dfsvisit"<<endl;
	
	dfsvisit();//call the dfsvisit function for generating the Unique SMILES

	return smilesstring;
}




/*

string Molecule::unique_smiles(int i)
{
	EvaluateInitialValue();	//Evaluate the initial value of all the atoms
	vector<int> originalrank;
	
	evaluate_rank();//evaluate rank after the initial value has been calculated
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}
	int old_unique_count;
	old_unique_count=distinctRankCount();//calculate the unique ranks count 
	
	originalrank=rank;
	vector<set<int> > PrevClasses;
	PrevClasses = EvaluateAtomClasses();
	int test_cond=1;
	vector<int>PrevRank;
	// the following piece of code tries to re-evaluate rank until final rank converges. 
	while (test_cond!=0)
	{
		PrevRank=rank;
		test_cond=0;//flag for testing rank convergence
		
		EvaluateValue();
		evaluate_rank();//re-evaluating rank
	
		ModEqRankUsingPrev(PrevRank);
		
		vector<set<int> >newClasses;
		newClasses = EvaluateAtomClasses();
		if (PrevRank!=rank)test_cond++;
		int new_unique_count;
		new_unique_count=distinctRankCount();

		if (test_cond!=0)
		{
			
			if (new_unique_count==size)test_cond=0;
		}
		if (test_cond!=0)
		{
			if(old_unique_count==new_unique_count)
			{
			
				if (IsEqualClasses(PrevClasses,newClasses))test_cond=0;				
			}
		}
		PrevClasses = newClasses;
		old_unique_count=new_unique_count;
		
	}
	cout<<""<<endl;
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}
	classRank=rank;
	
	if (distinctRankCount()!=size)//if there are ties, then unique rank count is less than size of hte molecules)
	{
		//Now we should break ties
		
		breakTies();	
		//cout<<"after first break ties"<<endl;
		for (int i=0;i<size;i++)//printing the final sorted adjacency list
		{
			cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<Hydrogens[i]<<"  ---->";
			for (int j=0;j<NN[i];j++)
			{
				cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
			}
			cout<<endl;
			//cout<<atoms[i].get_isotope_number()<<endl;
		}
		
		test_cond=1;
		while (test_cond!=0)
		{
			PrevRank=rank;
			EvaluateValue();
			evaluate_rank();
			ModEqRankUsingPrev(PrevRank);

			if (distinctRankCount()==size)test_cond=0;
			else 
			{
				breakTies();
				test_cond=1;
			}
		}
	}
	
	
	sort_adjacency();//sorting the adjacency list in order of ranks and double bonds
	
	cout<<"the sorted adjacency list"<<endl;
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}
	
	
	dfsvisit();//call the dfsvisit function for generating the Unique SMILES

	return smilesstring;
}

*/

void Molecule::findRanksOfAllAtoms()
{
	vector<int> originalrank;
	int old_unique_count;
	old_unique_count=distinctRankCount();//calculate the unique ranks count 
	
	originalrank=rank;
	vector<set<int> > PrevClasses;
	PrevClasses = EvaluateAtomClasses();
	int test_cond=1;
	vector<int>PrevRank;
	// the following piece of code tries to re-evaluate rank until final rank converges. 
	while (test_cond!=0)
	{
		PrevRank=rank;
		test_cond=0;//flag for testing rank convergence
		//cout<<"here"<<endl;
		EvaluateValue();
		evaluate_rank();//re-evaluating rank
	
		ModEqRankUsingPrev(PrevRank);
		
		vector<set<int> >newClasses;
		newClasses = EvaluateAtomClasses();
		//for (int r =0;r<PrevRank.size();r++)
		//	cout<<PrevRank[r]<<"  "<<rank[r]<<endl;
		if (PrevRank!=rank)test_cond++;
		int new_unique_count;
		new_unique_count=distinctRankCount();

		if (test_cond!=0) 
		{
			
			if (new_unique_count==size)test_cond=0;
		}
		if (test_cond!=0)
		{
			if(old_unique_count==new_unique_count)
			{
			
				if (IsEqualClasses(PrevClasses,newClasses))test_cond=0;	
				//cout<<"equal classes"<<endl;
			}
		}
		PrevClasses = newClasses;
		old_unique_count=new_unique_count;
		
	}
}

string Molecule::unique_smiles()
{
	
	
	EvaluateInitialValue();	//Evaluate the initial value of all the atoms

	
	evaluate_rank();//evaluate rank after the initial value has been calculated

	findRanksOfAllAtoms();
	//cout<<"done"<<endl;
	classRank=rank;

	//print_adjacency_list();
	
	while(distinctRankCount()!=size)//if there are ties, then unique rank count is less than size of hte molecules)
	{
		//Now we should break ties
		//cout<<"going to break ties"<<endl;
		breakTies();		
		//cout<<"finished ties"<<endl;
		//print_adjacency_list();
		findRanksOfAllAtoms();
		//cout<<"After calculating ranks"<<endl;
		//print_adjacency_list();

	}
	
	
	sort_adjacency();//sorting the adjacency list in order of ranks and double bonds
	
	
	/*cout<<"the sorted adjacency list"<<endl;
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<classRank[i]<<"  "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}*/
	//cout<<"going to dfsvisit"<<endl;
	
	dfsvisit();//call the dfsvisit function for generating the Unique SMILES

	return smilesstring;
}
	

/*				
string Molecule::unique_smiles()
{
	
	
	EvaluateInitialValue();	//Evaluate the initial value of all the atoms
	vector<int> originalrank;
	
	evaluate_rank();//evaluate rank after the initial value has been calculated
	int old_unique_count;
	old_unique_count=distinctRankCount();//calculate the unique ranks count 
	
	originalrank=rank;
	vector<set<int> > PrevClasses;
	PrevClasses = EvaluateAtomClasses();
	int test_cond=1;
	vector<int>PrevRank;
	// the following piece of code tries to re-evaluate rank until final rank converges. 
	while (test_cond!=0)
	{
		PrevRank=rank;
		test_cond=0;//flag for testing rank convergence
		//cout<<"here"<<endl;
		EvaluateValue();
		evaluate_rank();//re-evaluating rank
	
		ModEqRankUsingPrev(PrevRank);
		
		vector<set<int> >newClasses;
		newClasses = EvaluateAtomClasses();
		if (PrevRank!=rank)test_cond++;
		int new_unique_count;
		new_unique_count=distinctRankCount();

		if (test_cond!=0)
		{
			
			if (new_unique_count==size)test_cond=0;
		}
		if (test_cond!=0)
		{
			if(old_unique_count==new_unique_count)
			{
			
				if (IsEqualClasses(PrevClasses,newClasses))test_cond=0;				
			}
		}
		PrevClasses = newClasses;
		old_unique_count=new_unique_count;
		
	}
	//cout<<"done"<<endl;
	classRank=rank;
	
	if (distinctRankCount()!=size)//if there are ties, then unique rank count is less than size of hte molecules)
	{
		//Now we should break ties
		//cout<<"going to break ties"<<endl;
		breakTies();		
		//cout<<"finished ties"<<endl;
		
		test_cond=1;
		while (test_cond!=0)
		{
			//cout<<"inside while loop"<<endl;
			PrevRank=rank;
			EvaluateValue();
			evaluate_rank();
			ModEqRankUsingPrev(PrevRank);

			//cout<<"actually in here"<<endl;

			if (distinctRankCount()==size)test_cond=0;
			else 
			{
				breakTies();
				test_cond=1;
			}
		}
	}
	
	
	sort_adjacency();//sorting the adjacency list in order of ranks and double bonds
	
	
	cout<<"the sorted adjacency list"<<endl;
	for (int i=0;i<size;i++)//printing the final sorted adjacency list
	{
		cout<<i<<" ("<<NN[i]<<") "<<getatomtype(i)<<" "<<rank[i]<<" "<<classRank[i]<<"  "<<Hydrogens[i]<<"  ---->";
		for (int j=0;j<NN[i];j++)
		{
			cout<<Adjacency[i][j]<<"("<<BO[i][j]<<")"<<"("<<rank[Adjacency[i][j]]<<")"<<"  ";
		}
		cout<<endl;
		//cout<<atoms[i].get_isotope_number()<<endl;
	}
	//cout<<"going to dfsvisit"<<endl;
	
	dfsvisit();//call the dfsvisit function for generating the Unique SMILES

	return smilesstring;
	
	
}

*/
vector<int> Molecule::RankOrder(std::vector<int> & Atomlist, std::vector<int> & PrevRank) const
{
	vector<int> ModRanks;
	ModRanks.resize(Atomlist.size(),1);

	for (int i=0;i<Atomlist.size();i++)
	{
		for (int j=0;j<=i;j++)
		{
			if (PrevRank[Atomlist[j]]>PrevRank[Atomlist[i]])ModRanks[j]++;
			if (PrevRank[Atomlist[j]]<PrevRank[Atomlist[i]])ModRanks[i]++;
		}
	}
	return ModRanks;
}

void Molecule::ModEqRankUsingPrev(vector<int>& PrevRank)
{
	vector<bool>IsChecked;
	IsChecked.resize(size,false);
	

	for (int i=0;i<size;i++)
	{
		vector<int>Atomlist;
		if (!IsChecked[i])
		{
			for (int j=i+1;j<size;j++)
			{
				if (rank[j]==rank[i])
				{
					Atomlist.push_back(j);
					IsChecked[j]=true;
					
				}
			}
			if (!Atomlist.empty())
			{
				Atomlist.push_back(i);
				vector <int>ModRanks;
				
				ModRanks = RankOrder(Atomlist,PrevRank);
				
				for (int j=0;j<Atomlist.size();j++)
				{
					rank[Atomlist[j]]+=(ModRanks[j]-1);
				}
			}
			
		}
	}
}
	



void Molecule::print_smiles()//function to print out the smilesstring
{
	cout<<"the SMILES string is"<<endl;
	cout<<smilesstring<<endl;
}

void Molecule::set_smiles(string stringname)//setting the smiles string
{
	smilesstring=stringname;
}

void Molecule::dfsvisit()//visiting the graph structure using DFS approach
{
	
	int* color = new int[size];//the array that stores the color:0-whie;1-grey and 2-black
	int* parent = new int[size];//the array that stores the information of the parent
	int nodes_to_cover;//defining a counter for number of nodes to cover and initializing it to total number of atoms
	nodes_to_cover=size;
	int Cnode=-1;//Current node
	string newsmiles;//intermediate string that holds the unique smiles as DFS proceeds
	newsmiles="";
	int* adj_nodes_visited_count =new int[size];//each atoms has a counter to note the number of adjacent atoms visited
	int ring_count;//counter for the rings
	ring_count=0;
	IntPair* ring= new IntPair[100];//Pair to store the ring atom number
	vector<int>ringBond;
	for (int i=0;i<size;i++)
	{
		color[i]=0;//coloring all nodes white initially
		adj_nodes_visited_count[i]=0;
	}
	int lowest;
	lowest=rank[0];
	for (int i=1;i<size;i++)
	{
		if (lowest>rank[i])lowest=rank[i];
	}
	for (int i=0;i<size;i++)
	{
		if (rank[i]==lowest)
		{
			Cnode=i;//Finding the first node with rank 1 and declaring it as the current node. 
			break;
		}
	}
	
	parent[Cnode]=-1;//the root node does not have a parent. Hence the value -1.
	//cout<<"getting into while loop"<<endl;
	while(nodes_to_cover!=0)
	{
		//cout<<"nodes to cover: "<<nodes_to_cover<<endl;
		//cout<<"Cnode is "<<Cnode<<" visited nodes count "<<adj_nodes_visited_count[Cnode]<<endl;
		if (color[Cnode]==0)//checking if the color of the node is white
		{

			color[Cnode]=1;//if found, color it grey
			newsmiles=newsmiles+atoms[Cnode]->get_atom_symbol();//add the symbol of the atomtype
			if (adj_nodes_visited_count[Cnode]==NN[Cnode])
			{
				nodes_to_cover--;//if the adjacent bonds are all visited already, then color the node black
				color[Cnode]=2;
				if ((NN[Cnode]==1)&&(nodes_to_cover>0))
				{
					newsmiles=newsmiles+")";//if the node was a leaf and there are more nodes to cover, add ')' to indicate end of branch
				}
			}
		}
		
		if (color[Cnode]==1)//checking if the color of the node is grey
		{

			if (adj_nodes_visited_count[Cnode]<NN[Cnode])//check if the Cnode has unvisited bonds 
			{
				if (adj_nodes_visited_count[Cnode]==NN[Cnode]-1)//if the adjacent bond is the last bond to cover
				{
					color[Cnode]=2;//then color it black anyway because in the outer loop this adjacent bond will be traversed
					nodes_to_cover--;//So color the node black! 
				}
				
				
				if (color[Adjacency[Cnode][adj_nodes_visited_count[Cnode]]]==0)//checking if the adjacent node is white
				{
					parent[Adjacency[Cnode][adj_nodes_visited_count[Cnode]]]=Cnode;//declaring the Cnode as the parent of this adjacent node
					
					
					
					if (adj_nodes_visited_count[Cnode]<(NN[Cnode]-1))
						newsmiles=newsmiles+"(";//adding a '(' to indicate a branch if the adjacent bond is not the last to cover
					if (BO[Cnode][adj_nodes_visited_count[Cnode]]==2) 
						newsmiles=newsmiles+"=";//checking and adding double bond
					if (BO[Cnode][adj_nodes_visited_count[Cnode]]==3)
						newsmiles=newsmiles+"#";//checking and adding triple bond
					if (BO[Cnode][adj_nodes_visited_count[Cnode]]==5)
						newsmiles+="_";
					if (BO[Cnode][adj_nodes_visited_count[Cnode]]==1)
					{
						if (isaromatic(Cnode) && isaromatic(Adjacency[Cnode][adj_nodes_visited_count[Cnode]]))
							newsmiles+="-";
					}

					adj_nodes_visited_count[Cnode]++;//incrementing the visited adjacent bonds of Cnode
		
					Cnode=Adjacency[Cnode][adj_nodes_visited_count[Cnode]-1];//update Cnode as this adjacent node

					for (int i=0;i<NN[Cnode];i++)
					{
						if(Adjacency[Cnode][i]==parent[Cnode])
						{
							push_to_first(Cnode,i);//transferring the parent to the first position in adjacency list
							break;
						}
					}
					adj_nodes_visited_count[Cnode]++;//incrementing the visited adjacent bonds of new Cnode. 

				}
				else if(color[Adjacency[Cnode][adj_nodes_visited_count[Cnode]]]==1)//checking if the adjacent node is grey
				{
					ring_count++;//if the node is grey, then remove the node from adjacency list and maintain it in a ring Pair
					
					ring[ring_count-1].first=Cnode;
					ring[ring_count-1].second=Adjacency[Cnode][adj_nodes_visited_count[Cnode]];
					ringBond.push_back(BO[Cnode][adj_nodes_visited_count[Cnode]]);
					
					int adj=Adjacency[Cnode][adj_nodes_visited_count[Cnode]];//storing the adjacent node number in adj
					int BOValue=BO[Cnode][adj_nodes_visited_count[Cnode]];//storing the bond order value in adj
					//adding the ring identifiers in the symbol of the atomtype
					string symbol,symbol2;
					symbol=atoms[Cnode]->get_atom_symbol();
					symbol2=atoms[adj]->get_atom_symbol();
					if (BOValue==2)symbol+="=";
					if (BOValue==5)symbol+="_";
					if (BOValue==3)symbol+="#";
					char out[5];
					//_itoa_s(ring_count,out,5,10);
					sprintf(out, "%d", ring_count);
					symbol=symbol+out[0];
					symbol2=symbol2+out[0];
					if (ring_count>=10)
					{	symbol=symbol+out[1];
						symbol2=symbol2+out[1];
					}
					if (ring_count>=100)
					{
						symbol=symbol+out[2];
						symbol=symbol+out[2];
					}
					
					atoms[Cnode]->set_atom_symbol(symbol);
					atoms[adj]->set_atom_symbol(symbol2);
					//cout<<atoms[Cnode].get_atomtype_symbol()<<endl;
								

					for(unsigned int i=adj_nodes_visited_count[Cnode];i<NN[Cnode]-1;i++)//removing the node by pushing it to the end of adjacency list and decrementing NN
					{
						Adjacency[Cnode][i]=Adjacency[Cnode][i+1];
						BO[Cnode][i]=BO[Cnode][i+1];
					}
					NN[Cnode]--;//Decrementing NN of Cnode
					
					

					//locating Cnode in adj
					for (unsigned int i=0;i<NN[adj];i++)
					{
						if (Adjacency[adj][i]==Cnode)
						{
							for (unsigned int j=i;j<NN[adj]-1;j++)//locating the position that holds Cnode and pushing it to the last
							{
								Adjacency[adj][j]=Adjacency[adj][j+1];
								BO[adj][j]=BO[adj][j+1];
							}
							break;
						}
					}

					NN[adj]--;//decrementing NN of the adj
					
					//check if adj is now to be colored black
					if (adj_nodes_visited_count[adj]==NN[adj])
					{
						color[adj]=2;
						nodes_to_cover--;
					}
				}
			}
			else//if the node was colored grey but all its bonds have been visited, then color it black
			{
				if (color[Cnode]!=2)//ensuring if the node was not colored black already, thereby preventing over counting
				{
					color[Cnode]=2;
					nodes_to_cover--;
				}
			}
		}
		if (color[Cnode]==2)//if the Cnode is colored Black and if its not the rootnode, then reassign the Cnode as the parent of the current Cnode.
		{
			if (parent[Cnode]>=0)
			{
				Cnode=parent[Cnode];
			}
		}


	}
	set_smiles(newsmiles);//set the unique smiles
	if (ring_count>0)//if the ring counter showed prescence of rings, do the DFS visit the second time
	{
		//cout<<"performing second loop of DFSvisit"<<endl;
		dfsvisit();
	}

	for (unsigned int i=0;i<ring_count;i++)//Adjacency list and nearest neighbors are updated again to inculde broken bonds.
	{
		Adjacency[ring[i].first][NN[ring[i].first]]=ring[i].second;
		Adjacency[ring[i].second][NN[ring[i].second]]=ring[i].first;
		BO[ring[i].first][NN[ring[i].first]]=ringBond[i];
		BO[ring[i].second][NN[ring[i].second]]=ringBond[i];

		NN[ring[i].first]++;
		NN[ring[i].second]++;
	}
	
	//update_aromaticity_details();
	delete[] color;
	delete[] parent;
	delete[] ring;
	delete[] adj_nodes_visited_count;
	
}

void Molecule::print_rings()
{
	//cout<<"the number of rings in Molecule is: "<<Allrings.size()<<endl;
	Allrings.print_all();
}

void Molecule::print_aromatic_rings()
{
	if (aromatic_rings.size()>0)
	{
		//cout<<"the aromatic rings are the rings"<<endl;
		for (unsigned int i=0;i<aromatic_rings.size();i++)
		{
			//cout<<aromatic_rings[i]+1<<" ";
		}
		//cout<<endl;
	}
	else
	{
		//cout<<"there are no aromatic rings in the Molecule"<<endl;
	}
}
bool Molecule::isaromaticbond(int m,int n) const
{
	
	if (isaromatic(m) && isaromatic(n) && (find_BO(m,n)==4)) 
		return true;
	else 
		return false;
}

void Molecule::remove_Hydrogens()
{
	/*This function removes trims loose hydrogens in a molecule. It only trims neutral hydrogen atoms such as hydrogens of a methyl group. The atom to which the hydrogen was attached has its Hydrogens 
	value incremented by one*/
	vector <int> H_atoms;
	for (int i=0;i<size;i++)
	{
		if (getatomtype(i)=="H")
		{
			if (getBatomcount(i)==getBplusNBatomcount(i))//atom has no NB interactions-- round about, I know!
			{
				for ( int j=0;j<size;j++)
				{
					for (unsigned int k=0;k<Adjacency[j].size();k++)
					{
						if (Adjacency[j][k]==i && BO[j][k]==1)
						{
							Adjacency[j].erase(Adjacency[j].begin()+k);
							Hydrogens[j]++;
							break;
						}
					}
				}
				H_atoms.push_back(i);
			}
		}
	}
	int last_atom = 0;
	if (H_atoms.size()==size)//if the molecule is H or HH
		last_atom++;
	//cout<<"found all Hydrogens"<<endl;
	for (int i=H_atoms.size()-1;i>=last_atom;i--)
	{
		erase(H_atoms[i]);
	}
	//cout<<"erased Hydrogens"<<endl;
	for (int i=0;i<size;i++)
	{
		NN[i]=Adjacency[i].size();
	}
}

bool Molecule::ishydrogenic() const
{
	if (size<2)
	{
		if (atoms[0]->get_atomic_number()==1)
			return true;
		else return false;
	}
	else return false;
}


IntPair Molecule::Hydrogen_counter(int i) const
{
	int D_count=0;
	
	for (int j=0;j<NN[i];j++)
	{
		if (atoms[Adjacency[i][j]]->get_atom_symbol()=="D")	D_count++;
	}
	IntPair P;
	P.first=Hydrogens[i];
	P.second=D_count+ Hydrogens[i];
	return P;
}

bool Molecule::isaromaticmolecule() const
{
	if (aromatic_rings.size()>0)
		return true;
	else
		return false;
}

bool Molecule::iscyclicmolecule() const
{
	if (Allrings.size()>0)
		return true;
	else
		return false;
}

bool Molecule::isringatom(int i) const
{
	if (Allrings.is_present(i))
		return true;
	else
		return false;
}



int Molecule::MaxRingsize() const
{
	if (Allrings.size()>0)
		return Allrings.get(0).size();
	else return 0;
}

int Molecule::SmallestRingSize() const
{
	if (Allrings.size()>0)
		return Allrings.get(Allrings.size()-1).size();
	else return 0;

}

int Molecule::MaxRingSizeOfAtom(int i) const
{
	return Allrings.SizeContainingAtom(i);
}
string Molecule::moleculestring() const
{
	return smilesstring;
}

bool Molecule::isneutral() const
{
	int flag=0;
	//cout<<smilesstring<<endl;
	for (unsigned int i=0;i<size;i++)
	{
		
		if (atoms[i]->get_charge()>0)
		{
			//cout<<"positive charge at "<<i<<endl;
			//cout<<"charge is "<<atoms[i].get_charge()<<endl;
			flag=1;

			break;
		}
	}
	
	if (flag==0) return true;
	else return false;
}
int Molecule::totalcharge() const
{
	int counter=0;
	for (unsigned int i=0;i<size;i++)
	{
		
		counter+=atoms[i]->get_charge();
	}
	return counter;
}

int Molecule::totalupElectrons() const
{
	int counter = 0;
	for (unsigned int i=0;i<size;i++)
	{
		counter+=atoms[i]->get_up();
	}

	return counter;
}

void Molecule::Readjustproperties()
{
	//this function checks for the cases of N,S,P to see if the initially set valency is less than NN. If greater than NN, 
	//it means that the atom is in its higher oxidation state. Else, it is in its lower state. 
	
	for (unsigned int i=0;i<size;i++)
	{
		atoms[i]->readjustatomproperties(getBatomcount(i), Hydrogens[i], dbcount(i), tpcount(i));
			
	}
}

void Molecule::calculateHydrogens()
{
	for (unsigned int i=0;i<size;i++)
	{
		int HydrogenCount = 0;
		HydrogenCount = atoms[i]->get_valency()-dbcount(i)-2*tpcount(i)-getBatomcount(i);
		
		
		if (atoms[i]->get_atomtype_name().compare(0,1,"C"))//if carbon
		{
			if (isaromatic(i))//if aromatic, one less count of Hydrogens. if three aromatic bonds already, then no Hydrogen possible!
			{
				HydrogenCount--;
				if (AromaticBondCount(i)==3)HydrogenCount=0;
			}
		}
		else
		{
			if (isaromatic(i))
			{
				if (Hydrogens[i]==0)HydrogenCount--;
				else HydrogenCount=Hydrogens[i];
				

			}
			
			//for aromatic systems with declared hydrogens, I assume they are correct; 
		}    // For other non aromatic cases, the higher oxidation state has already been calculated and HydrogenCount shows correct value. 

		if (HydrogenCount<0)HydrogenCount=0;//this will ensure i dont get negative value of hydrogen for aromatic oxygens.
		Hydrogens[i]=HydrogenCount;
	}
}

bool Molecule::IsNeighbor(int i, int j) const
{
	bool answer = false;
	for (unsigned int k=0;k<NN[i];k++)
	{
		if (Adjacency[i][k]==j)
		{
			answer=true;
			break;
		}
	}
	return answer;
}

bool Molecule::InSameRing(int i,int j) const
{
	
	bool answer = false;
	if (Allrings.is_present(i) && Allrings.is_present(j))
	{
		
		
		for (unsigned int k=0;k<Allrings.size();k++)
		{
			
			if ((Allrings.get(k).contains(i)) && (Allrings.get(k).contains(j)))
			{
				answer = true;
				break;
			}
		}
		
	}
	
	return answer;
}

bool Molecule::isringbond(int i, int j) const
{
	if (IsNeighbor(i,j) && (find_BO(i,j)!=5) && InSameRing(i,j))
		return true;
	else return false;
}
int Molecule::AdjacentBondInRingWithOrder(Path& P, int i, int j) const
{
	int count=0;
	for (unsigned int k=0;k<P.size();k++)
	{
		if((find_BO(i,P.get(k))==j))count++;		
	}
	return count;
}

bool Molecule::RingWithNBInteractions(Path& P) const
{
	bool nb = false;
	for (unsigned int i=0;i<P.size()-1;i++)
	{
		for (unsigned int j=i+1;j<P.size();j++)
		{
			if((find_BO(P.get(i),P.get(j))==5))
			{
				nb=true;
				break;
			}
		}
	}
	return nb;
}


string Molecule::SetWithoutHydrogenCount(string s)
{
	string newstring;
	bool InsideCurlyBraces = false;
	if (s.compare(0,1,"[")==0)
	{
		for (unsigned int i=0;i<s.length();i++)
		{
			
			if (s.compare(i,1,"{")==0)InsideCurlyBraces = true;
			if (s.compare(i,1,"}")==0)InsideCurlyBraces = false;
			
			if ((i>=2) && s.compare(i,1,"H")==0 && !InsideCurlyBraces)
			{				
				while ((s.compare(i,1,"]"))!=0)
					i++;
			}
			newstring+=s[i];
		}
		
	}
	else newstring=s;

	return newstring;
	
}

string Molecule::SetWithHydrogenCount(string s, int count)
{
	
	string newstring;
	if (count>0)
	{

		
		char HCount[5];
		
		//_itoa_s(count,HCount,5,10);
		sprintf(HCount, "%d", count);

		for (unsigned int i=0;i<s.length();i++)
		{
			if (s.compare(i,1,"]")!=0)
			{
				newstring+=s[i];
			}
			else
			{
				int index = 0;
				string toAdd="H";
				//cout<<count<<endl;
				if (count>1)
					toAdd+= HCount[0];
				
				if (s.compare(i-1,1,"H")==0)index = 1;
				else 
				{
					if (s.compare(i-2,1,"H")==0)index = 2;
				}
				newstring.resize(i-index);

				/*
				if (toAdd.length()==1)
					previous = s[i-1];
				else //ensuring that if the correct Hcount was already included, then we should ignore! 
				{
					previous+=s[i-2];
					previous+=s[i-1];
					//cout<<"previous is "<<previous<<endl;
				}
				
				if (previous.compare(toAdd)!=0)
				{
					newstring+=toAdd;
				}*/
				newstring+= toAdd;

				
				newstring+="]";
				//cout<<toAdd<<endl;
				//cout<<newstring<<endl;
			}
		}
				
	}
	else
		newstring=s;
	
	return newstring;


}

string Molecule::SetWithoutSquareBrackets(std::string s)
{
	string newstring;
	for (unsigned int i=0;i<s.length();i++)
	{
		if ((s.compare(i,1,"[")!=0) && (s.compare(i,1,"]")!=0))
			newstring+=s[i];
	}
	return newstring;
}

string Molecule::SetWithSquareBrackets(string s)
{
	if (s.compare(0,1,"[")!=0)
		return ("[" + s +"]");
	else return s;
}


				

void Molecule::UpdateSquareBraces()
{
	for (int i=0;i<size;i++)
	{
		if (atoms[i]->get_nature()==0)
		{
			string atomsymbol = atoms[i]->get_atom_symbol();
			string atomname = atoms[i]->get_atomtype_name();
			//cout<<"isotope number "<<atoms[i]->get_isotope_number()<<endl;

			
			
			if (atomname.compare(0,1,"H")==0)
			{
				if (atomname.compare("H")!=0)
				{
					atoms[i]->set_atom_symbol(SetWithSquareBrackets(atomsymbol));
				}
				else 
				{
					if (atoms.size()==1)atoms[i]->set_atom_symbol("[HH]");
				}
			}

			if ((atomname.compare(0,1,"C")==0) || (atomname.compare(0,1,"O")==0))
			{
				
				
				atomsymbol=SetWithoutHydrogenCount(atomsymbol);
				if (atomname.compare("C")==0 || atomname.compare("O")==0)
				{
					atomsymbol=SetWithoutSquareBrackets(atomsymbol);
				}
				else
					atomsymbol=SetWithSquareBrackets(atomsymbol);
						
				atoms[i]->set_atom_symbol(atomsymbol);
			}
			
				
			if ((atomname.compare(0,1,"N")==0) || (atomname.compare(0,1,"P")==0))
			{
				//cout<<"its Nitrogen!"<<endl;
				//if (isaromatic(i)) 
					//cout<<"its also aromatic"<<endl;
				//the only cases where Nitrogen needs to be 
				if (atomname.compare("N")==0|| atomname.compare("P")==0)
				{
					atomsymbol=SetWithoutHydrogenCount(atomsymbol);
					atoms[i]->set_atom_symbol(SetWithoutSquareBrackets(atomsymbol));
				}
				else
				{
					atoms[i]->set_atom_symbol(SetWithSquareBrackets(atomsymbol));
				}

				
				if ((isaromatic(i)) ||(atoms[i]->get_valency()>3 && getBatomcount(i)<=3))
				{
					//cout<<"Hydrogens is "<<Hydrogens[i]<<endl;
					if ((Hydrogens[i]>0))
					{
						//cout<<"it is aromatic and has hydrogens"<<endl;
						atomsymbol=SetWithSquareBrackets(atomsymbol);//SetWithHydrogenCount requires that square brackets are already there
						atoms[i]->set_atom_symbol(SetWithHydrogenCount(atomsymbol,Hydrogens[i]));
					}
					else atoms[i]->set_atom_symbol(SetWithoutSquareBrackets(SetWithoutHydrogenCount(atomsymbol)));
				}
				
			}
			if (atomname.compare(0,1,"S")==0)
			{
				if ((atoms[i]->get_valency()>2 && getBatomcount(i)<=2) ||(atoms[i]->get_valency()>4 && getBatomcount(i)<=4))
				{
					if ((Hydrogens[i]>0))
					{
						atomsymbol=SetWithSquareBrackets(atomsymbol);
						atoms[i]->set_atom_symbol(SetWithHydrogenCount(atomsymbol,Hydrogens[i]));
					}
					else atoms[i]->set_atom_symbol(SetWithoutSquareBrackets(SetWithoutHydrogenCount(atomsymbol)));
				}
			}

			if (atoms[i]->get_isotope_number()>0 && atoms[i]->get_isotope_number()!=atoms[i]->get_n_mass_number())
			{
				int massnumber = atoms[i]->get_isotope_number();
				char Mcount[5];
				//cout<<"massnumber is "<<massnumber<<endl;
				//_itoa_s(massnumber,Mcount,5,10);
				sprintf(Mcount, "%d", massnumber);
				string isotopestring ="";
				//cout<<Mcount[0]<<endl;
				//cout<<Mcount[1]<<endl;

				if (massnumber >=1)isotopestring+=Mcount[0];
				if (massnumber >=10) isotopestring+=Mcount[1];
				if (massnumber >=100)isotopestring+=Mcount[2];
				atomsymbol = SetWithoutSquareBrackets(atomsymbol);
				if (atomsymbol.compare(0,isotopestring.length(),isotopestring)!=0)
				{
					atomsymbol = isotopestring + atomsymbol;
				}
				//cout<<"atomsymbol is "<<atomsymbol<<"isotopestring is "<<isotopestring<<endl;
				atoms[i]->set_atom_symbol(SetWithSquareBrackets(atomsymbol));
				//cout<<atoms[i]->get_atom_symbol()<<endl;
			}
			
		

		}
		else
		{
			string atomsymbol = atoms[i]->get_atom_symbol();
			atomsymbol=SetWithSquareBrackets(atomsymbol);
				
			if (Hydrogens[i]>0)
				atoms[i]->set_atom_symbol(SetWithHydrogenCount(atomsymbol,Hydrogens[i]));
			else 
				atoms[i]->set_atom_symbol(SetWithoutHydrogenCount(atomsymbol));
		
		}
	}
}


void Molecule::EvaluateMF()
{
	vector<int>AtomCounters;//0-C, 1- H, 2-O, 3-N, 4-S, 5-P
	//vector<int>IsotopeAtomCounters;//order same as above
	AtomCounters.resize(6);
	//IsotopeAtomCounters.resize(6);
	vector<string>ElementList;
	ElementList.push_back("C");
	ElementList.push_back("H");
	ElementList.push_back("O");
	ElementList.push_back("N");
	ElementList.push_back("S");
	ElementList.push_back("P");
	map<string,int>CompositesSet;
	set<string>NonNeutralAtoms;
	set<string>IsotopicAtoms;
	int HydrogenCounter = 0;
	for (int i=0;i<size;i++)
	{
		
		HydrogenCounter+=Hydrogens[i];
		if (atoms[i]->get_nature()==0)
		{
			if (atoms[i]->get_isotope_number()==0)
			{
				
				
				if (atoms[i]->get_atomtype_name().compare("C")==0 || atoms[i]->get_atomtype_name().compare("C:")==0)AtomCounters[0]++;
				else if (atoms[i]->get_atomtype_name().compare("H")==0)AtomCounters[1]++;
				else if (atoms[i]->get_atomtype_name().compare("O")==0)AtomCounters[2]++;
				else if (atoms[i]->get_atomtype_name().compare("N")==0)AtomCounters[3]++;
				else if (atoms[i]->get_atomtype_name().compare("S")==0)AtomCounters[4]++;
				else if (atoms[i]->get_atomtype_name().compare("P")==0)AtomCounters[5]++;
				else NonNeutralAtoms.insert(atoms[i]->get_atomtype_name());
			}
			else IsotopicAtoms.insert(atoms[i]->get_atom_symbol());
		}
		else
		{
			if (CompositesSet.count(atoms[i]->get_atomtype_name())==0)
				CompositesSet[atoms[i]->get_atomtype_name()]=1;
			else
				CompositesSet[atoms[i]->get_atomtype_name()]+=1;

		}
	}
	AtomCounters[1]+=HydrogenCounter;
	
	MolecularFormula ="";
	for (int i=0;i<AtomCounters.size();i++)
	{
		
		if (AtomCounters[i]>0)
		{
			char ACount[5];
			//_itoa_s(AtomCounters[i],ACount,5,10);
			sprintf(ACount, "%d", AtomCounters[i]);
			MolecularFormula+=ElementList[i];
			if (AtomCounters[i]>1)MolecularFormula+=ACount[0];
			if (AtomCounters[i]>9)MolecularFormula+=ACount[1];
			if (AtomCounters[i]>99)MolecularFormula+=ACount[2];
		}
	}
	set<string>::iterator it;
	map<string,int>::iterator m_it;
	for (it=IsotopicAtoms.begin();it!=IsotopicAtoms.end();it++)
	{
		//(*it).insert(1,"'");
		//MolecularFormula+=(*it);
		for (int k=0;k<(*it).length();k++)
		{
			MolecularFormula+=toupper((*it)[k]);
		}
	}
	for (m_it=CompositesSet.begin();m_it!=CompositesSet.end();m_it++)
	{
		MolecularFormula+=m_it->first;
		if (m_it->second>1)MolecularFormula+= IntToStr(m_it->second);
	}
	for (it=NonNeutralAtoms.begin();it!=NonNeutralAtoms.end();it++)
	{
		MolecularFormula+=(*it);
	}
}


void Molecule::find_allylic_atoms()
{
	for (int i=0;i<size;i++)
	{
		string S;
		S = atoms[i]->get_atomtype_name();
		if (S=="C" && IsAdjacentAtomWithBO(i,"C",2))
		{
			for (int j=0;j<NN[i];j++)
			{
				string S2 = atoms[Adjacency[i][j]]->get_atomtype_name();
				if ((S2=="C" || S2=="C+" || S2=="C-" || S2=="C." || S2=="C:")&& BO[i][j]==1 && dbcount(Adjacency[i][j])==0)
				{
					allylic_atoms.insert(Adjacency[i][j]);
				}
			}
		
		}

	}
		
}

string Molecule::GetMF() const
{
	return MolecularFormula;
}

string Molecule::getMFwithoutAtoms(set<string> S) const
{
	string MolFormulaWithoutAtoms="";
	bool isRestrictedAtom = false;
	for (int i =0;i<MolecularFormula.length();i++)
	{
		string current = MolecularFormula.substr(i,1);
		
		if (S.count("{" + current + "}")==0)
		{
			if (isalpha(MolecularFormula[i]) || (isdigit(MolecularFormula[i]) && !isRestrictedAtom))
			{
				MolFormulaWithoutAtoms+=MolecularFormula[i];
				isRestrictedAtom = false;
			}
		}
		else
		{
			isRestrictedAtom = true;
		}

	}
	return MolFormulaWithoutAtoms;
}		
		

pair<string,int> Molecule::checkValencyError() const
{

	pair<string,int> P;
	P.first="";
	P.second=0;
	for (int i=0;i<size;i++)
	{
		int count = getBatomcount(i)+ Hydrogens[i]+ dbcount(i)+ 2*tpcount(i);
		//cout<<getBatomcount(i)<<" "<<Hydrogens[i]<<" "<<dbcount(i)<<"  "<<tpcount(i)<<endl;
		int valency = atoms[i]->get_valency();
		//cout<<"atomtype is "<<atoms[i]->get_atomtype_name()<<" "<<count<<" "<<valency<<endl;
		if (count!=valency)
		{
			if (!isaromatic(i))
			{
				if (atoms[i]->get_atomtype_name().compare("N+")==0 && count>valency)
					P = pair<string,int> (atoms[i]->get_atomtype_name(),4);
				else
				{
					P = pair<string, int>(atoms[i]->get_atomtype_name(),1);//throws an error that a nonaroamtic atom has a bond count-valency mismatch! 
					//cout<<count<<" "<<valency<<" "<<dbcount(i)<<" "<<getBatomcount(i)<<endl;
				}
				break;
			}
			else
			{
				if (atoms[i]->get_atomtype_name().compare(0,1,"C")==0)
				{
					if (count!=valency-1)
						P = pair<string,int> (atoms[i]->get_atomtype_name(),6);
				}
				else
				{
					if (Hydrogens[i]!=0)
					{
						if (atoms[i]->get_atomtype_name().compare("N+")==0 && count>valency)
							P = pair<string,int> (atoms[i]->get_atomtype_name(),4);
						else
							P = pair<string, int>(atoms[i]->get_atomtype_name(),2);//throws an error that the hydrogen count of an aromatic atom is wrongly given! 
						break;
					}
					else
					{
						if (atoms[i]->get_atomtype_name().compare("N+")==0 && count>valency-1)
							P = pair<string,int> (atoms[i]->get_atomtype_name(),5);
						else
							P= pair<string, int>(atoms[i]->get_atomtype_name(),3);//throws an error that an aromatic atom has a bond count (inclusive of aromatic bonds)-valency mismatch
						break;
					}
				}	
			}
		}
	}
	return P;
}

pair<int,int> Molecule::calcBranchandDistanceValue() const
{
	int brvalue = 1;
	int distance = 1;
	vector<vector<int> > path;
	path.resize(size);
	vector<int> ringatoms;
	vector<int> farthestDaughterleaf;
	//Initialize path details

	//finding shortest path length between two nodes
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<size;j++)
		{
			if (i==j)path[i].push_back(0);
			int bovalue = find_BO(i,j);
			if (bovalue>0) path[i].push_back(1);
			if (bovalue<0) path[i].push_back(1000000);//large number
		}
	}
	for (int k=0;k<size;k++)
	{
		for (int i=0;i<size;i++)
		{
			for (int j=0;j<size;j++)
			{
				int a = path[i][k] + path[k][j];
				if (path[i][j]>a) path[i][j]=a;
			}
		}
	}

	//finding the leaves 
	vector<int> leaves;
	for (int i=0;i<size;i++)
	{
		if (Adjacency[i].size() ==1)leaves.push_back(i);
	}


	if (!iscyclicmolecule()) //if acyclic
	{
		
		for (int k=0;k<leaves.size();k++)
		{
			for (int j=k+1;j<leaves.size();j++)
			{
				brvalue = brvalue * path[leaves[k]][leaves[j]];//calculate brvalue as product of distance between any two leaves
			}
		}
	}
	else
	{
		//if cyclic - we look at branch points along the ring - so we find which is the closest ring atom to that leaf and then seek distances between these parent ring atoms. 
		for (int i=0;i<size;i++)
		{
			if (isringatom(i))ringatoms.push_back(i);
		}
		
		vector<int> ClosestRingAtomOfEachLeaf; //stores the closest ring atom of each leaf - indexed according to the order in leaves
		//find the closest ring atom for each leaf.
		for (int k=0;k<leaves.size();k++)
		{
			int closestringatom = ringatoms.front();

			for (int i = 0;i<ringatoms.size();i++)
			{
				if (path[leaves[k]][ringatoms[i]]< path[leaves[k]][closestringatom])
				{
					closestringatom =i;
				}
			}
			ClosestRingAtomOfEachLeaf.push_back(closestringatom);
		}

		//distance value is calculated as product of distance between each two closestringatom

		for (int k=0;k<ClosestRingAtomOfEachLeaf.size();k++)
		{
			for (int j=k+1;j<ClosestRingAtomOfEachLeaf.size();j++)
			{
				distance = distance * path[ClosestRingAtomOfEachLeaf[k]][ClosestRingAtomOfEachLeaf[j]];
			}
		}

		//this calculates the branch value as the product of distances of each leaf from its closest ring atom

		for (int i =0;i<ClosestRingAtomOfEachLeaf.size();i++)
		{
			brvalue = brvalue * path[leaves[i]][ClosestRingAtomOfEachLeaf[i]];
		}
	}



		

	return pair<int,int>(brvalue,distance);
}

int Molecule::totalDoubleBonds() const
{
	int total =0;
	for (int i=0;i< size;i++)
	{
		total+=dbcount(i);
	}

	return total/2;
}

int Molecule::totalTripleBonds() const
{
	int total = 0;
	for (int i=0; i< size; i++)
	{
		total+= tpcount(i);
	}
	return total/2;
}

int Molecule::totalAromaticBonds() const
{
	int total = 0;
	for (int i=0; i< size; i++)
	{
		total+= AromaticBondCount(i);
	}
	return total/2;
}
	

bool Molecule::ishydrocarbon() const
{
	bool ishydrocarbon = true;

	for (int i=0;i<size;i++)
	{
		if (atoms[i]->get_atomtype_name().compare(0,1,"C")!=0)
		{
			ishydrocarbon = false;
			break;
		}
	}
	return ishydrocarbon;
}


bool Molecule::isparaffinic() const
{
	return (!iscyclicmolecule() && totalDoubleBonds()==0 && totalTripleBonds()==0 && ishydrocarbon());
}

bool Molecule::isolefinic() const
{
	return (!iscyclicmolecule() && totalDoubleBonds()>0 && totalTripleBonds()==0 && ishydrocarbon());
}

bool Molecule::isNaphthenic() const
{
	return (ishydrocarbon() && iscyclicmolecule() && !isaromaticmolecule());
}

bool Molecule::isIntermediate() const
{
	for (int i=0;i<size;i++)
	{
		//if molecule has a charge, unpaired electron or has non-bonded interactions, then it is an intermediate
		if (atoms[i]->get_charge()!=0 || atoms[i]->get_up()!=0 ||HasNBinteractions())
			return true;
	}
	return false;
}

bool Molecule::HasCompositeAtoms() const
{

	for (int i=0;i<size;i++)
	{
		if (atoms[i]->get_nature()==1)
		{
			return true;
		}
	}
	return false;
}

bool Molecule::HasBridgeHead() const
{
	bool value = false;
	for (int i=0;i<Allrings.size()-1;i++)
	{
		Path P1 = Allrings.get(i);
		int counter =0;
		for (int j=i+1 ;j<Allrings.size();j++)
		{
			Path P2 = Allrings.get(j);

			for (int k=0;k<P2.size();k++)
			{
				if (P1.contains(P2.get(k)))counter++;
			}

			if (counter<P1.size() && counter<P2.size() && counter>2)
			{
				value = true;
				break;
			}
		}
		if (value)break;
	}
	return value;
}

int Molecule::NumberOfAromaticRings() const
{
	return aromatic_rings.size();
}

int Molecule::NumberOfRings() const
{
	return Allrings.size();
}

int Molecule::totalAtomsIncludingH() const
{
	int totalCount = 0;

	for (int i=0;i<atoms.size();i++)
	{
		totalCount+=getHydrogens(i)+1;
	}
	return totalCount;
}

int Molecule::totalAtomsOfType(string S) const
{
	int totalCount = 0;

	for (int i=0;i<atoms.size();i++)
	{
		if (atoms[i]->get_atomtype_name().compare(S)==0)totalCount++;
	}
	return totalCount;
}

int Molecule::MolecularWeight() const
{
	int weight = 0;

	for (int i=0;i<atoms.size();i++)
	{
		if (atoms[i]->get_isotope_number()==0)
			weight+=atoms[i]->get_n_mass_number();
		else weight+=atoms[i]->get_isotope_mass_number();
		weight+=getHydrogens(i);
	}

	return weight;
}

bool Molecule::isAtomChiral(int a) const
{
	if (a<0 || a>=size)
	{
		cout<<"incorrect value of atom index "<<a<<" used to test for chirality"<<endl;
		return false;
	}
	else
	{
		set<int> distinctClasses;
		//TODO: checking only for neutral Carbon here - need to see it we have to check for other elements as well
		if (getatomtype(a)=="C")
		{
			for (int i=0;i<NN[a];i++)
				distinctClasses.insert(getclassRank(Adjacency[a][i]));
		
			if (getHydrogens(a)<=1 && distinctClasses.size()==4)
				return true;
			else return false;
		}
		else return false;
	}
}

int Molecule::NumberChiralAtoms() const
{
	int numChiralAtoms = 0;
	for (int i=0;i<size;i++)
		if (isAtomChiral(i))numChiralAtoms++;

	return numChiralAtoms;

}

int Molecule::OptIsomers() const
{
	map<int,int> ClassAndChiralAtomNumberMap;//stores the number of chiral atoms of a particular topological class 
	int NumEnantiomers = 1;
	for (int i =0; i<size;i++)
	{
		if (isAtomChiral(i))
		{
			NumEnantiomers = 1;//setting it to 1 because it's atleast 1(re-setting each time in the loop, that's OK --update: why? )
			if (ClassAndChiralAtomNumberMap.count(classRank[i])>0)
				ClassAndChiralAtomNumberMap[classRank[i]]+=1;
			else
				ClassAndChiralAtomNumberMap[classRank[i]]=1;
		}
	}
	
	for (map<int,int>::iterator it = ClassAndChiralAtomNumberMap.begin(); it!=ClassAndChiralAtomNumberMap.end();it++)
	{
		if (it->second==1)
			NumEnantiomers=NumEnantiomers*2;
		else if (it->second%2==0)//if even
			NumEnantiomers=NumEnantiomers*(2*it->second);//if even it's 2n
		else
			NumEnantiomers=NumEnantiomers*(2*it->second-2);//if odd it's 2n-2
	}
	return NumEnantiomers;
}


int Molecule::NumberOfLeaves() const 
{
	int leaves = 0;
	for (int i =0;i<size;i++)
	{
		if (NN[i]==1)leaves++;
	}
	return leaves;
}
	


	




	





		


			
	

			



		

				



		



				









		


//end of implementation of Molecule class





