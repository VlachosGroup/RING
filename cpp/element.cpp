#include "element.h"

using namespace std;

//implementations of Element class
Element::Element(char elementname)
{
	atomic_number=n_mass_number=isotope_mass_number=0;
	n_symbol=isotope_symbol="";
	n_valency=n_lp=n_up=0;
	set_element_properties(elementname);
}

void Element::set_element_properties(char elementname)
{
	switch(elementname)
	{
	case 'C':
	case 'c':
		atomic_number=6;
		n_mass_number=12;
		isotope_mass_number=13;
		n_symbol="C";
		isotope_symbol="C'";
		n_valency=4;
		n_lp=0;
		n_up=0;
		n_states=1;
		break;

	case 'o':	
	case 'O':
		atomic_number=8;
		n_mass_number=16;
		isotope_mass_number=18;
		n_symbol="O";
		isotope_symbol="O'";
		n_valency=2;
		n_lp=2;
		n_up=0;
		n_states=1;
		break;

	case 'n':
	case 'N':
		atomic_number=7;
		n_mass_number=14;
		isotope_mass_number=15;
		n_symbol="N";
		isotope_symbol="N'";
		n_valency=3;
		n_lp=1;
		n_up=0;
		n_states=2;
		break;

	case 'H':
		atomic_number=1;
		n_mass_number=1;
		isotope_mass_number=2;
		n_symbol="H";
		isotope_symbol="D";
		n_valency=1;
		n_lp=0;
		n_up=0;
		n_states=1;
		break;

	case 'S':
	case 's':
		atomic_number=16;
		n_mass_number=32;
		isotope_mass_number=34;
		n_symbol="S";
		isotope_symbol="S'";
		n_valency=2;
		n_lp=2;
		n_up=0;
		n_states=3;
		break;

	case 'P':
	case 'p':
		atomic_number=15;
		n_mass_number=31;
		isotope_mass_number=31;
		n_symbol="P";
		isotope_symbol="P'";
		n_valency=3;
		n_lp=1;
		n_up=0;
		n_states=2;
		break;
	}
}
//End of Element implementation


//Implementation of Atomtype
Atomtype::Atomtype(string stringname, char elementname)
		:Element(elementname)	
{
	atomtype_name = stringname;
	lp = 0;
	up = 0;
	charge = 0;
}
//end of implementation of Atomtype


