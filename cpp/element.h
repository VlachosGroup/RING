#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>

class Element
{
	protected:
		int atomic_number;//atomic number
		int n_mass_number;//Normal Mass number
		int isotope_mass_number;//Isotope Mass number
        std::string n_symbol;//normal element symbol
        std::string isotope_symbol;//isotope  symbol
		int n_valency;//normal valency
		int n_lp;//normal lone pairs
		int n_up;//normal unpaired electrons
		int n_states;//the number of oxidation states (1 for C,H,O; 2 for N,P; 3 for S)
	public:
		Element(char);//Constructor of element class
		Element();//default constructor
		
		void set_element_properties(char);//sets element properties
};


class Atomtype:public Element 
{
	protected:
        std::string atomtype_name;
		int lp; //the actual no: of lone pairs
		int up; //the actual no: of unpaired electrons
		int charge; //actual charge
		int valency; //actual valency
	public:
		Atomtype(std::string, char);
		Atomtype(); 
};

#endif
