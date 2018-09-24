#ifndef ATOMCONTAINER_H
#define ATOMCONTAINER_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <vector>
//#include <list>
#include <set>
#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;

#include "common.h"
//#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
#include "atom.h"


class Atomcontainer//Atomcontainer class. Has the information on adjacency and bond order
{
	protected:
		std::vector<Atom*> atoms;//std::vector of atoms
		std::vector< std::vector <int> > Adjacency;//Adjacency list
		std::vector < std::vector<int> > BO;//Bond order list //Note 1- single bond, 2- double, 3- triple, 4-aromatic and 5 - nonbonded
		std::vector<int> NN;//Nearest neighbor (counter) - nonbonded is also counted! 
		std::vector<float> value; //value
		std::vector<int> rank;// rank
		std::vector<int> classRank;//
		int size;// size of the container
		std::vector<int> label;//label of the atoms
		std::vector< std::vector< int > > H_label;//container to store the labels of the Hydrogens
		std::vector<int> Hydrogens;//number of Hydrogens 
		Ringset Allrings;//Stores the Set of all rings
		std::vector <int> aromatic_atoms;//std::vector containing all aromatic atoms
		void EvaluateInitialValue();//evaluates the initial value of the atoms
		void EvaluateInitialValue(int);
		void EvaluateValue();//evaluates the value of the atoms 
		void evaluate_rank();// evaluates the rank of the atoms 
		int distinctRankCount() const;//outputs the number of unique ranks (classes)
		void reverse_ranks();//reverses the order of the ranks
		void sort_adjacency();//sorts the adjacency list based on rank
		void push_to_first(int,int);//a useful manipulating function to push an element in the adjacency list to the first
		void erase(int);//erasing an atom from the atomcontainer
		Atomcontainer form_atomcontainer(std::vector<int>) const;//create an atomcontainer with a list of atoms.
		bool IsEqualClasses(std::vector<std::set<int> >, std::vector<std::set<int> >) const;//checks if two atom class classifications are equal
		void breakTies();//if two atoms have the same rank, then this function forcibly breaks the tie
		std::vector<std::set<int> > EvaluateAtomClasses();//evaluates the atom classes and puts them into sets of equivalent topological classes
		
		
	public:
		Atomcontainer(int);//constructor
		Atomcontainer();
		~Atomcontainer();
		Atomcontainer(const Atomcontainer &a);
		Atomcontainer& operator=(const Atomcontainer &a);
		void setrank(std::vector<int>& );//set rank function
		int getrank(int) const;//gets the rank 
		int getclassRank(int) const;//gets the rank of the class the atom belongs to
		int get_adjacency(int,int) const;//get the adjacency list entries
		Atom* getatom(int) const;//get the atom with the given index
		int get_BO(int,int) const;//get BO entries - this seeks the particular BO entry. the two ints are coordinates! 
		int getsize() const;//get size
		int getNN(int) const;//get NN value
		int dbcount(int) const;//gets double bond count
		int tpcount(int) const;//gets triple bond count
		int find_BO(int,int) const;//find the BO between two atoms;NOTE: gives -1 if there is no bond between the two! 
		void print_adjacency_list() const;
		bool isaromatic(int) const;//checks if atom is aromatic
        std::string getatomtype(int) const;//getting atomtype of a particular atom
		void merge(const Atomcontainer&);//merge the atomcontainer with another
		int findatomwithlabel(int) const;//find the atom with a given label - returns -1 if none is found! 
		int getlabel(int) const;//get the label of a particular atom
		void setlabel(int, int);//set the label of a particular atom to a specified value
		float getvalue(int) const;//get the value of a particular atom
		void formbond(int, int);//form bond between two specified atoms
		void breakbond(int, int);//break bond between two specified atoms
		void changeBO(int, int, int);//change BO between two specified atoms by a particular value
		void setBO(int, int, int);//Set BO between two specified atoms to a particular value
		void setatomtypename(int, std::string);//set atomtype name of a particular atom
		void setatomsymbol(int, std::string);//set atom symbol of a particular atom
		void setInitialAtomproperties(int);//set initial atom properties of specific atom
		void setatomvalency(int);//sets the valency of the atom
		void resetlabel(int);//sets all the label to zero.
		int getHydrogens(int) const;//gives the hydrogen count of the given atom
		void setHydrogens(int, int);//sets the hydroen count of a given atom to second argument's value;
		std::vector<int> getAromaticAtoms() const;//gives the std::vector of aromatic atoms of the atomcontainer
		std::vector<int> getHlabel(int) const;//gives the std::vector of hydrogen labels of a given atom
		void setHlabel(int, std::vector<int>);//sets the H label of a given atom with a std::vector of hydrogen labels.
		int findatomwithHlabel(int) const;//finds the atom that has the given argument as a Hydrogen label - finds -1 if none is found
		void changeHlabel(int, int);//change the hydrogen label given by the first argument to the value given in the second argument
		int getBplusNBatomcount(int) const;//gets the count of both bonded and nonbonded atoms that are neighboring to the given atom
		int getBatomcount(int) const;//gets the bonded atom count of the given atom
		int AromaticBondCount(int) const;//gets the count of the aromatic bonds attached to the atom
		bool HasNBinteractions() const;//checks if there are nonbonded interactions in the atomcontainer
		bool IsAdjacentAtomWithBO(int, std::string, int) const;//checks if the atom corresponding to the first argument is adjacent to an atom of the given string, with a bond order of the given third argument
        std::pair<std::vector<Atomcontainer>, std::vector< std::vector<int> > > connectedcomponents() const;//generates a std::vector of atomcontainers that are essentially the connected components of the parent atomcontainer and also the atomindices of these components
		int getElectronicHashValue(int) const;//returns positive charge ->2^charge; negative -> 3^charge; neutral ->1; radical -> 5; positive radical ->7^value; used in  
		int getTotalElectronicValue() const;//returns product of getElectronicHashValue for all the atoms
		void removeHlabel(int, int);//removes the Hlabel (second argument) from the parent atom (first argument)
		void addHlabel(int, int);//add the Hlabel (second argument) to the parent atom (first argument)
		void addAtom(std::string, std::string, int, char, int);//add a new atom -atomtypename, atomsymbol, isotopenumber, elementname, and nature
		int getGroupHash(int) const;//gets the hash value (or invariant) of atom specified that takes into consideration itself and its immediate neighbors
		int getNNElecHash(int) const;//gets nearest neighbors' electronic info as a hash value		
		//int getNNBondCount(int, int) const;//gets the number of bonds adjacent to the specified atom (first argument) of specified BOtype (Second argument), not counting the bonds to the specified atom. For eg (1,2) means finding double bond count of atoms adjacent to atom 1.
		int getNNDoubleBondCount(int) const;//gets the number of double bonds adjacent to the specified atom
		int getNNTripleBondCount(int) const;//gets the number of triple bonds adjacent to the specified atom
		int getNNDoubleTripleBondFactors(int) const;//gets a hash factor to account for double and triple bonds! 
		int RingCountOfAtom(int) const;//finds the number of rings an atom is in	
        std::map<int, int> getClassesFreqMapNeighboringAtom(int) const;//for specified atom, it creates a map of the different classes and their frequencies
        std::string AtomCenteredGroupForGA(int) const;//generates the atom-centered group of specified atom needed for group!
		int AtomValueForGA(std::string, bool) const;
        std::set<std::string> GetElements() const;
		friend class Patternmatch;
};

#endif
