#include <iostream>
#include <string>
#include <vector>

#include "substructure.h"

using std::string; using std::vector;
using std::cout; using std::endl;


//Substructure Implementation

Substructure::Substructure(string fragment_string, int m)
		:Atomcontainer(m)
{
	//NEED TO CLEANUP CODE
	fragmentstring=fragment_string;
	string atomtypename;
	string atomsymbol;
	atomtypename="";
	atomsymbol="";
	int parent=-1;
	int r=-1;
	int bond_ident=1;
	AtomEnv.resize(m);
	atom_flag.resize(m);
	int NumberOfRings = 0;
	int ringcounter = 0;
	bool IsAfterRingLabelSign = false;
	Triplet* ring=new Triplet[100];//ring Triplet declaration
	vector < string > ringnumber;
	for (int i=0;i<m;i++)
	{
		atom_flag[i].resize(5,-1);
	}
	vector <int>branch;
	int H_flag=0;
	
	H_label.resize(m);
	int ringflag=0;
	int notflag=0;
	bool IsCompositeAtomFound =false;
	bool ExplicitSingleBond=false;
	bool ExplicitDoubleBond=false;
	bool InsideSingleQuotes = false;

	for (unsigned int i=0;i<fragmentstring.length();i++)
	{
		string c;
		c=fragmentstring[i];
		//cout<<c<<endl;
		if (c=="'")
		{
			atomsymbol+="'";
			if (!InsideSingleQuotes)
			{
				InsideSingleQuotes = true;
				
				i++;
				c=fragmentstring[i];
			}
			else InsideSingleQuotes = false;
		}
		
		if ((isalpha(fragmentstring[i]) && c!="R") || c=="$" || c=="&" || c=="{")
		{
			
			NumberOfRings=0;
			IsAfterRingLabelSign = false;
			if (c=="{")
			{
				IsCompositeAtomFound=true;
				atomsymbol+="{";
				i++;
				c=fragmentstring[i];
				while (c!="}")
				{
					atomtypename+=fragmentstring[i];
					atomsymbol+=fragmentstring[i];
					i++;
					c=fragmentstring[i];
				}
				atomsymbol+="}";
			}
			else
			{
				atomsymbol=c;
				atomtypename=toupper(fragmentstring[i]);
			}

			
			if ((parent>=0) && (c!="H"))
			{
				H_flag=0;
				r++;
				
				Adjacency[r].push_back(parent);
				Adjacency[parent].push_back(r);
				BO[r].push_back(bond_ident);
				BO[parent].push_back(bond_ident);

				bond_ident=1;
				if (ringflag==1)
				{
					Triplet Ringinfo;
					Ringinfo.first=parent;
					Ringinfo.second=r;
					ringflag=0;
					if (notflag==1)
					{
						Ringinfo.third=0;
					}
					else Ringinfo.third=1;
					notflag=0;
					ringbondDefs.push_back(Ringinfo);
				}

				string parentsymbol;
				parentsymbol=atoms[parent]->get_atom_symbol();
				if ((islower(atomsymbol[0])|| (atomsymbol.length()>1 && islower(atomsymbol[1])) || atom_flag[r][0]==1) && (islower(parentsymbol[0]) || (atom_flag[parent][0]==1)))
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
			
			if (c=="H" && parent>=0)
			{
				H_flag=1;
				Hydrogens[parent]++;
				
			}
			
			if (parent==-1)
			{
				parent=0;
				r=0;
			}
			else
			{
				if (H_flag!=1) parent=r;
			}
		}

		if (c=="R")
		{
			NumberOfRings++;
			ringnumber.push_back("");
			ringcounter++;
			IsAfterRingLabelSign=true;
		}


		if (isdigit(fragmentstring[i]))
		{
			if (H_flag==0 && !IsAfterRingLabelSign)
			{
				if (label[r]>=0)label[r]=label[r]*10+atoi(c.c_str());
				else label[r]=atoi(c.c_str());
			}
			else if (IsAfterRingLabelSign)
			{
				ringnumber[NumberOfRings-1]+=fragmentstring[i];//add this character to ringnumber
				if (bond_ident==2)ringnumber[NumberOfRings-1]+="=";
				if (bond_ident==3)ringnumber[NumberOfRings-1]+="#";
				if (bond_ident==5)ringnumber[NumberOfRings-1]+="_";
				if (bond_ident==4)ringnumber[NumberOfRings-1]+=":";
			}
			else
				H_label[parent].push_back(atoi(c.c_str()));

		}
		

		if (c=="[")
		{
			i++;
			c=fragmentstring[i];
			if (c=="{")
			{			
				i++;
				c=fragmentstring[i];
				int atomflagcheck=1;
				bool ringsizeDef = false;
				int ringsizeOpr = 0;
				int ringsizeValue = 0;
				while(c!="}")
				{
				
					
					if (c=="s")ringsizeDef = true;
					
					if (c=="!") atomflagcheck=0;
					if (isdigit(fragmentstring[i]) && !ringsizeDef)
					{
						atom_flag[r][atoi(c.c_str())]=atomflagcheck;
					}

					
					if (c=="<" && ringsizeDef)
					{
						ringsizeOpr = 1;
					}
					if (c==">" && ringsizeDef)ringsizeOpr = -1;
					if (c=="=" && ringsizeDef)
					{
						if (ringsizeOpr<0)ringsizeOpr = -2;
						if (ringsizeOpr>0)ringsizeOpr = 2;
					}
					if (ringsizeDef && isdigit(fragmentstring[i]))
					{
						ringsizeValue = atoi(c.c_str());
						atom_flag[r][3]=ringsizeOpr;
						atom_flag[r][4]=ringsizeValue;
					}		

					i++;
					c=fragmentstring[i];
				}
				while (c!="]")
				{
					i++;
					c=fragmentstring[i];
				}
			}
			else
			{
				env_tuple et;
				et.set_flag(1);
				et.set_freq(1);
				et.set_comp_opr(0);
				
				int DepthSqBr=1;
				
				while(c!="]")
				{
					if (c=="|")
					{					
						
						i++;
						c=fragmentstring[i];
						string typestring;
						typestring="";
						while (!(c=="|" && DepthSqBr==1))
						{
							
							typestring+=c;
							
							i++;
							c=fragmentstring[i];
							if (c=="[")DepthSqBr++;
							if (c=="]")DepthSqBr--;
						}
						et.set_type(typestring);
						//cout<<typestring<<endl;
						
					}
																			
					if (c=="!")
					{
						et.set_flag(0);
					}
					if (isdigit(fragmentstring[i]))
					{
						et.set_freq(atoi(c.c_str()));
					}
					if (c==">")
					{
						et.set_comp_opr(1);
					}
					if (c=="<")
					{
						et.set_comp_opr(-1);
					}
					if (c=="=")
					{
						
						if (et.get_comp_opr()==1)
						{
							et.set_comp_opr(2);
						}
						else if (et.get_comp_opr()==-1)
						{
							et.set_comp_opr(-2);
						}
						else et.set_comp_opr(0);
						
						
					}
					i++;
					
					c=fragmentstring[i];
				}
				
				AtomEnv[r].add(et);
				
			}
		}
		if (c=="(")
		{
			branch.push_back(parent);
		}
		if (c==")")
		{
			parent=branch.back();
			branch.pop_back();
		}
		if ((c=="+")||(c=="*")||(c=="."))//checking atomtypes and updating the atomtypename and atomsymbol appropriately
		{
			atomsymbol=atomsymbol+fragmentstring[i];
			atomtypename=atomtypename+fragmentstring[i];
		
		}
		if (c==":")
		{
			if (InsideSingleQuotes)
			{
				atomsymbol=atomsymbol+fragmentstring[i];
				atomtypename=atomtypename+fragmentstring[i];
			}
			else
			{
				bond_ident=4;
				ExplicitSingleBond = false;
				ExplicitDoubleBond = false;
			}
		}

		if (c=="-")
		{
			if (InsideSingleQuotes)
			{
				atomsymbol=atomsymbol+fragmentstring[i];
				atomtypename=atomtypename+fragmentstring[i];
			}
			else
			{
				bond_ident=1;
				ExplicitSingleBond = true;
			}
		}
	
		
		if (c=="=")//identifying double bonds
		{
			bond_ident=2;
			ExplicitDoubleBond = true;
	
		}
		if (c=="#")//identifying triple bonds
		{
			bond_ident=3;
		
		}
		if (c=="~")
		{
			bond_ident=6;
		}
		if (c=="%")
		{
			bond_ident=7;
		}
		if (c=="@")
		{
			ringflag=1;
		}
		if (c=="!")
		{
			notflag=1;
		}
		if (c=="_")
		{
			bond_ident = 5;
		}
		

		
		string cc;
		if (i<fragmentstring.length()-1)
			cc=fragmentstring[i+1];

		if (((i<fragmentstring.length()-1 && isalpha(fragmentstring[i+1]) && cc!="R") || (cc=="{") || (cc=="'") ||(cc=="$") || (cc=="&") || (i==(fragmentstring.length()-1)))&& !InsideSingleQuotes)
		{
			if (((r>=0) && atomtypename!="H")||(atoms.size()==0 && atomtypename=="H"))
			{
				//cout<<r<<" adding atom "<<atomtypename<<endl;
				if (IsCompositeAtomFound)
				{
					atoms.push_back(new CompositeAtom(atomtypename, atomsymbol));
				}
				else
				{
					char elementname=atomtypename[0];
					
					atoms.push_back(new SingleAtom(atomtypename,atomsymbol,0,elementname));//instantiating an object of Atom class
									
				}
				
				setInitialAtomproperties(r);//setting properties of the atomtype
				atoms[r]->set_valency();//setting valency of the atom (atomtype)
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
							if (bovalue==":")bo=4;
						}
					}
					
					l=atoi(ringvalue.c_str());//conversion of ringnumber value from string to integer
					
					
					ring[l-1].second=ring[l-1].first;
					ring[l-1].first=r;//updating the ring pair
					ring[l-1].third=bo;
					
					
					
				}
			}
			atomtypename="";
			atomsymbol="";
			ringnumber.clear();
		}
	}

	

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
		Triplet Ringinfo;
		Ringinfo.first=j2;
		Ringinfo.second=j3;
		Ringinfo.third = 1;
		ringbondDefs.push_back(Ringinfo);

			
	}

	for (int i=0;i<size;i++)
	{
		NN[i]=Adjacency[i].size();
	}

	for (int i=0;i<size;i++)//if atom_flag corresponding to aromaticitiy is 1, then add this to aromatic_atoms list
	{
		if (atom_flag[i][0]==1)
			aromatic_atoms.push_back(i);
	}



	delete[] ring;

}

string Substructure::getstring()
{
	return fragmentstring;
}


int Substructure::isringbondcheck(int i, int j) const
{
	int returnvalue=-1;
	for (int k=0;k<ringbondDefs.size();k++)
	{
		if (((ringbondDefs[k].first==i)&&(ringbondDefs[k].second==j)) ||((ringbondDefs[k].second==i) &&(ringbondDefs[k].first==j)))
		{
			returnvalue=ringbondDefs[k].third;
			break;
		}
	}
	return returnvalue;
}

void Substructure::printstring()
{
	cout<<fragmentstring<<endl;
}
	
