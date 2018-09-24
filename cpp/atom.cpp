#include <iostream>
#include <string>

#include "atom.h"

using std::cout; using std::endl;
using std::string;

//implementation of CompositeAtom class

CompositeAtom::CompositeAtom(string name, string symbol)
				:Atom(1)
{
	atom_name =name;
	atom_symbol=symbol;
	composite_element_name="";
	string c;
	for (int i=1;i<symbol.length();i++)
	{
		c = symbol[i];
		if (c !="}" && c!="{" && c!="]" && c!="H")composite_element_name+=symbol[i];
	}
	valency=0;
	charge=0;
}

string CompositeAtom::get_element_name()
{
	return composite_element_name;
}


string CompositeAtom::get_atom_symbol()
{
	return atom_symbol;
}



void CompositeAtom::set_atom_symbol(string S)
{
	atom_symbol = S;
}

void CompositeAtom::set_charge(int i)
{
	charge=i;
}

int CompositeAtom::get_charge()
{
	return charge;
}

int CompositeAtom::get_valency()
{
	return valency;
}
void CompositeAtom::set_valency()
{
	valency=6;
}

void CompositeAtom::set_atomtype_name(std::string S)
{
	atom_name=S;
}
string CompositeAtom::get_atomtype_name()
{
	return atom_name;
}
void CompositeAtom::reset_properties()
{
	valency=charge=0;
}


void CompositeAtom::set_initial_properties()
{
	for (int i=1;i<atom_name.length();i++)
	{
		string c;
		c=atom_name[i];
		if (c=="+")
		{
			charge+=1;
		}
		if (c=="-")
		{
			charge=charge-1;
		}
		if (c=="*")
		{
			charge=1;
		}
	}
}

void CompositeAtom::set_valency_value(int i)
{
	valency=i;
}
void CompositeAtom::readjustatomproperties(int i, int h, int d, int t)
{
	set_valency_value(i+h+d+2*t);
}

int CompositeAtom::get_lp()
{
	return 0;
}

int CompositeAtom::get_up()
{
	return 0;
}

void CompositeAtom::DropRingIdentifier()
{
	
	bool InsideSqBrackets = false;
	bool InsideCurlyBrackets=false;
	for (int i=0;i<atom_symbol.length();i++)
	{
		string S;
		S=atom_symbol[i];
		if (S=="{")InsideCurlyBrackets=true;
		if (S=="}")InsideCurlyBrackets=false;
		if (S=="[")InsideSqBrackets=true;
		if (S=="]")InsideSqBrackets=false;
		
		if ((!InsideSqBrackets)&& (!InsideCurlyBrackets) &&((S=="%") ||(S=="=") || (S=="#") || (S=="_")|| isdigit(atom_symbol[i])))
		{
			atom_symbol.resize(i);
			cout<<S<<"  "<<atom_symbol<<endl;
			
			break;
		}
	}
}



//end of implementation of CompositeAtom

//implementation of SingleAtom class

SingleAtom::SingleAtom(string& name, string& symbol, int number, char elementname)
	:Atomtype(name,elementname), Atom(0)//constructor of class SingleAtom
{
	atom_symbol=symbol;
	isotope_value=number;
	isotope_mass_number = isotope_value;
}


void SingleAtom::set_atom_symbol(string stringname)//setting atom_symbol
{
	atom_symbol=stringname;
}

string SingleAtom::get_atom_symbol()//retrieving atomtype symbol
{
	return atom_symbol;
}

string SingleAtom::get_element_name()//retrieve element name
{
	return n_symbol;
}

void SingleAtom::set_isotope_number(int number)//setting isotope number (flag)
{
	isotope_value=number;
}

int SingleAtom::get_isotope_number()//retrieving isotope number
{
	return isotope_value;
}

void SingleAtom::set_valency_value(int i)
{
	valency=i;
}

void SingleAtom::readjustatomproperties(int n, int h, int d, int t)
{
	if ((n+h+d+2*t)>valency)
	{
		for (int j=0;j<n_states-1;j++)
		{
			valency=valency+2;
			lp--;
			if (valency>=(n+h+d+2*t))break;
		}
	}
}
void SingleAtom::DropRingIdentifier()
{
	
	
	bool InsideSqBrackets = false;
	for (int i=0;i<atom_symbol.length();i++)
	{
		string S;
		S=atom_symbol[i];
		
		if (S=="[")InsideSqBrackets=true;
		if (S=="]")InsideSqBrackets=false;
		
		
		if ((!InsideSqBrackets)&&((S=="%") ||(S=="=") || (S=="#") || (S=="_")|| isdigit(atom_symbol[i])))
		{
			atom_symbol.resize(i);
			
			break;
		}
	}
}

void SingleAtom::reset_properties()
{
	lp=0;
	up=0;
	charge=0;
}

void SingleAtom::set_valency() // setting valency of the SingleAtom
{

	valency=n_valency-charge-(lp-n_lp)*2-(up-n_up);
	
	if (atomtype_name.length()>1)
	{
		string c;
		c=atomtype_name[1];	
		
		if (c=="*")
		{
			valency=valency+2;
			
		}
	}
	if ((atomtype_name=="$") || (atomtype_name=="%") || (atomtype_name=="X"))
		valency=6;

}

int SingleAtom::get_valency()//function for retrieving the valency
{
	return valency;
}

int SingleAtom::get_charge()//function for retrieving the charge
{
	return charge;
}

void SingleAtom::set_initial_properties()//function for setting the properties
{
	lp=n_lp;
	up=n_up;
	charge=0;
	for (int i=1;i<atomtype_name.length();i++)
	{
		string c;
		c=atomtype_name[i];
		if (c=="+")
		{
			if (lp!=0)
			{	
				lp=lp-1;
				
			}
			charge+=1;
		}
		if (c=="-")
		{
			lp=lp+1;
			charge=charge-1;
		}
		if (c=="*")
		{
			charge=1;
		}
		if (c==".")
		{
			up=1;
		}
		if (c==":")
		{
			lp+=1;
		}
	}
}

void SingleAtom::set_atomtype_name(string stringname)//function for setting atomtype_name
{
	atomtype_name="";
	atomtype_name=stringname;
}

string SingleAtom::get_atomtype_name()//function for retrieving teh atomtype_name
{
	return atomtype_name;
}
 
int SingleAtom::get_up()
{

	return up;
}
int SingleAtom::get_lp()
{
	return lp;
}

int SingleAtom::get_atomic_number()//retrieving the atomic number
{
	return atomic_number;
}
int SingleAtom::get_n_mass_number()//retrieving the normal mass number
{
	return n_mass_number;
}
int SingleAtom::get_isotope_mass_number()//retrieving the isotope mass number
{
	return isotope_mass_number;
}

void SingleAtom::set_charge(int i)
{
	charge = i;
}


//end of SingleAtomclass implementation

//start of Atom implementation

Atom::Atom(int i)
{
	nature =i;
}

int Atom::get_nature()
{
	return nature;
}

void Atom::set_nature(int i)
{
	nature = i;
}


void Atom::MakeSymbolAliphatic()
{
	if (nature ==0)
	{
		string symbol = get_atom_symbol();
		for (int i=0;i<symbol.length();i++)
		{
			if (islower(symbol[i]))
			{
				symbol[i]=toupper(symbol[i]);
				break;
			}
		}
		set_atom_symbol(symbol);
	}

}

void Atom::MakeSymbolAromatic()
{
	if (nature ==0)
	{
		string symbol = get_atom_symbol();
		for (int i=0;i<symbol.length();i++)
		{
			if (isupper(symbol[i]) && (symbol.compare(i,1,"H")!=0))
			{
				symbol[i]=tolower(symbol[i]);
				break;
			}
		}
		set_atom_symbol(symbol);
	}
}

bool Atom::IsSymbolAliphatic()
{
	bool answer = true;
	if (nature ==0)
	{
		string symbol = get_atom_symbol();
		for (int i=0;i<symbol.length();i++)
		{
			if (islower(symbol[i]))
			{
				answer =false;
				break;
			}
		}
	}
	return answer;
}

//End of Atom implementation

