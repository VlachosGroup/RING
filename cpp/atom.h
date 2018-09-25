#ifndef SINGLEATOM_H
#define SINGLEATOM_H

#include <string>

#include "clonable.h"
#include "element.h"

class Atom:public clonable
{
	protected:
		int nature;//0- Single atom; 1- Composite atom
	public:
		Atom(int);
		Atom(){nature = 0;}
		virtual ~Atom(){}
		virtual Atom* clone(){return new Atom( *this );}
		virtual std::string get_atomtype_name(){return "";}
		virtual int get_valency(){return 0;}
		virtual std::string get_atom_symbol(){return "";}
		virtual void set_atomtype_name(std::string){}
		virtual void set_atom_symbol(std::string){}
		virtual void set_initial_properties(){}
		virtual void set_valency(){}
		virtual void set_valency_value(int){}
		virtual void set_charge(int){}
		virtual int get_lp(){return 0;}
		virtual int get_up(){return 0;}
		virtual int get_charge(){return 0;}
		int get_nature();
		void set_nature(int);
		virtual std::string get_element_name(){return"";}
		void MakeSymbolAromatic();
		void MakeSymbolAliphatic();
		bool IsSymbolAliphatic();
		virtual int get_isotope_number(){return 0;}
		virtual int get_isotope_mass_number(){return 0;}
		virtual int get_n_mass_number(){return 0;}
		virtual int get_atomic_number(){return 0;}
		virtual void reset_properties(){}
		virtual void DropRingIdentifier(){}
		virtual void readjustatomproperties(int, int,int, int){}
};


class SingleAtom:public Atomtype, public Atom//Class SingleAtom derived from Atomtype. has information on atom symbol and a value to determine the isotope
{
	protected:
        std::string atom_symbol;//atom symbol. Will also include the square braces wherever appropriate.
		int isotope_value;//0 for normal and the value for other cases.
	public:
		SingleAtom(std::string&, std::string&, int, char);//Constructor
		SingleAtom();//default constructor
		~SingleAtom(){}
		virtual Atom* clone(){ return new SingleAtom( *this );}
		void set_isotope_number(int);//sets the isotope number flag
	    void set_atom_symbol(std::string);//sets teh atom_symbol
        std::string get_atom_symbol();//gets the atom_symbol
		int get_isotope_number();//gets isotope_number 
		void set_valency_value(int);//sets the value of the valency
        std::string get_element_name();//gets the element symbol basically
		void set_valency();//set valency
		int get_valency();//get valency
		int get_charge();//get charge
		void set_initial_properties();//set properties
		void set_atomtype_name(std::string);//sets the atomtype_name
        std::string get_atomtype_name();//get atomtype_name
		int get_up();//get the number of unpaired electrons
		int get_lp();//get the number of lone pair of electrons
		int get_atomic_number();//gets atomic number
		int get_n_mass_number();//sets normal mass number
		int get_isotope_mass_number();//gets isotope mass number
		void set_charge(int);//sets atom charge
		void reset_properties();//resets all atomtype properties to zero.
		void DropRingIdentifier();//drops the ring identified off the symbol
		void readjustatomproperties(int, int,int, int);// used in Readjustproperties() in molecule
};


class CompositeAtom:public Atom
{
	protected:
        std::string composite_element_name;//this is without { and }. Just the name
        std::string atom_name;//this includes the curly braces. eg. {Pt}.
        std::string atom_symbol;//this would include the square braces too! eg. [{Pt}+]. 
		int charge;
		int valency;
	public:
		CompositeAtom(std::string, std::string);
		CompositeAtom();
		~CompositeAtom(){}
		virtual Atom* clone(){return new CompositeAtom( *this );}
		void set_atomtype_name(std::string);
        std::string get_atomtype_name();
		void set_valency();
		int get_valency();
		void set_valency_value(int);
		void set_charge(int);//sets charge
		int get_charge();
		void set_atom_symbol(std::string);
        std::string get_atom_symbol();
		void set_initial_properties();
		int get_lp();
		int get_up();
        std::string get_element_name();
		int get_isotope_number(){return 0;}
		int get_isotope_mass_number(){return 0;}
		int get_n_mass_number(){return 0;}
		int get_atomic_number(){return 0;}		
		void reset_properties();
		void DropRingIdentifier();
		void readjustatomproperties(int, int,int, int);
};


#endif
