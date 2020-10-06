#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <deque>
#include <algorithm>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


enum AstNodeType
{
	Undefined, OperatorAND, OperatorOR,OperatorNOT,Boolean
};

enum TokenType
{
	Error, AND, OR, NOT, EndOfText, Openbraces, Closedbraces, Size, Charge, IsAromatic, IsCyclic, MinRingSize, MaxRingSize, Fragment
};

class IntPair//IntPair is a class to store pair of elements
{
	public:
		int first;
		int second;
};
class Triplet//Triplet is a class to store element triplets
{
	public:
		int first;
		int second;
		int third;
};

class Path//Path is a collection of connected atoms
{
	protected:
		vector < int > atom_array;//array of atoms
	public:
		void add(int);//add an element
		int get(int) const;//get an element
		void reverse();//reverse the contents of the array
		void print();//print the array
		int size() const;//size 
		int getfirst() const;//the first element
		int getlast() const;//the last element
		void remove(int);//remove an element
		void join(Path,int);//join two paths
		int intersection(Path) const;//taking the intersection of two sets of atoms constituting the two paths
		void clear();//remove the contents of the array
		bool contains(int) const;// checks if a particular element is in the Path
		
};

class Ringset//collection of rings
{
	protected:
		vector < Path > ring_list;//array of paths
	public:
		void add(Path);//add a path
		Path get(int) const;//get a path
		void print(int);//print particular ring path
		void print_all();//print all rings
		void remove(int);//remove a path
		int size() const;//number of paths
		void sizesort();//sort by the size of the paths from largest to smallest
		bool is_present(int) const;
		int SizeContainingAtom(int) const;
		int CountAtom(int) const;
};

class env_tuple//triplet containing info on atom environment
{
	protected:
		
		int flag;//a flag to indicate presence or absence - 0 is absent and 1 is present
		string type;//string that defines the atomtype (in case of specifications involving single atoms) or SMARTS string (in case of environment constraints being a group)
		int freq;//frequency or the number of occurences
		int comp_opr;//comparison operator (0 for equality, 1 for greater than, -1 for less than)
		
	public:
		void set_flag(int);//set the flag value
		int get_flag();//get the flag value
		void set_type(string);//set the atomtype
		string get_type();//get the atomtype
		int get_freq();
		void set_freq(int);
		void set_comp_opr(int);
		int get_comp_opr();
		
};

class env_set//a set of environment triplets
{
	protected:
		vector<env_tuple>environment;//array of triplets
	public:
		void add(env_tuple);//add a triplet
		env_tuple get(int);//get a triplet
		void remove(int);//remove a triplet
		bool set_is_empty();//check if the array is empty
		int size();//number of triplets - the size of environment array
};

class bond_operations//contains info on the bond changes
{
	public:
		int first_atom_label;//first of the two atoms describing hte bond
		int second_atom_label;//second of the two atoms describing the bond
		int set_bond_order;//flag to set the bond order: -1 implies inactive, 0 is disconnect (break) and 1 is connect (form bond)
		int change_bond_order;//increase or decrease bond order implied by -1 and 1 for reducing the BO by 1/ increasing it by 1. 
							  // A value of 0 would denote no change in bond order. Note that the set_bond_order and change_bond_order cannot hold independent values
							  //Note that setting the bond order to 1 and change_bond_order to 1 is invalid, for example. 
};

class clonable {

	public:
		virtual ~clonable() {}
		virtual clonable* clone() = 0;
};

struct classcomp {
  bool operator() (string* a, string* b) const
  {return (*a)<(*b);}
};


struct paircomp {
	bool operator() (const pair<unsigned int, unsigned int>& lhs, const pair<unsigned int,unsigned int>& rhs)
	{
		if (lhs.first<rhs.first) return true;
		else if (rhs.first==lhs.first)
			return lhs.second<rhs.second;
		else return false;
	}
};

class Element
{
	protected:
		int atomic_number;//atomic number
		int n_mass_number;//Normal Mass number
		int isotope_mass_number;//Isotope Mass number
		string n_symbol;//normal element symbol
		string isotope_symbol;//isotope  symbol
		int n_valency;//normal valency
		int n_lp;//normal lone pairs
		int n_up;//normal unpaired electrons
		int n_states;//the number of oxidation states (1 for C,H,O; 2 for N,P; 3 for S)
	public:
		Element(char);//Constructor of element class
		Element();//default constructor
		
		void set_element_properties(char);//sets element properties
};

class Atomtype:public Element //Atomtype derived from Element. Has information on the form of the atom
{
	protected:
		string atomtype_name;// name of the atomtype. 
		int lp;//the actual no: of lone pairs
		int up;//the actual no: of unpaired electrons
		int charge;//actual charge
		int valency;//actual valency
	public:
		Atomtype(string, char);//Constructor of Atomtype
		Atomtype();//default constructor
		
};

class Atom:public clonable
{
	protected:
		int nature;//0- Single atom; 1- Composite atom
		
		
	public:
		Atom(int);
		Atom(){nature = 0;}
		virtual ~Atom(){}
		virtual Atom* clone(){return new Atom( *this );}
		virtual string get_atomtype_name(){return "";}
		virtual int get_valency(){return 0;}
		virtual string get_atom_symbol(){return "";}
		virtual void set_atomtype_name(string){}
		virtual void set_atom_symbol(string){}
		virtual void set_initial_properties(){}
		virtual void set_valency(){}
		virtual void set_valency_value(int){}
		virtual void set_charge(int){}
		virtual int get_lp(){return 0;}
		virtual int get_up(){return 0;}
		virtual int get_charge(){return 0;}
		int get_nature();
		void set_nature(int);
		virtual string get_element_name(){return"";}
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

class CompositeAtom:public Atom
{
	protected:
		string composite_element_name;//this is without { and }. Just the name
		string atom_name;//this includes the curly braces. eg. {Pt}.
		string atom_symbol;//this would include the square braces too! eg. [{Pt}+]. 
		int charge;
		int valency;
		
		
	public:
		CompositeAtom(string, string);
		CompositeAtom();
		~CompositeAtom(){}
		virtual Atom* clone(){return new CompositeAtom( *this );}
		void set_atomtype_name(string);
		string get_atomtype_name();
		void set_valency();
		int get_valency();
		void set_valency_value(int);
		void set_charge(int);//sets charge
		int get_charge();
		void set_atom_symbol(string);
		string get_atom_symbol();
		void set_initial_properties();
		int get_lp();
		int get_up();
		string get_element_name();
		int get_isotope_number(){return 0;}
		int get_isotope_mass_number(){return 0;}
		int get_n_mass_number(){return 0;}
		int get_atomic_number(){return 0;}		
		void reset_properties();
		void DropRingIdentifier();
		void readjustatomproperties(int, int,int, int);
};


class SingleAtom:public Atomtype, public Atom//Class SingleAtom derived from Atomtype. has information on atom symbol and a value to determine the isotope
{
	protected:
		string atom_symbol;//atom symbol. Will also include the square braces wherever appropriate.
		int isotope_value;//0 for normal and the value for other cases.
		
		
	public:
		
		SingleAtom(string&,string&,int,char);//Constructor
		SingleAtom();//default constructor
		~SingleAtom(){}
		virtual Atom* clone(){ return new SingleAtom( *this );}
		void set_isotope_number(int);//sets the isotope number flag
	    void set_atom_symbol(string);//sets teh atom_symbol
		string get_atom_symbol();//gets the atom_symbol
		int get_isotope_number();//gets isotope_number 
		void set_valency_value(int);//sets the value of the valency
		string get_element_name();//gets the element symbol basically
		void set_valency();//set valency
		int get_valency();//get valency
		int get_charge();//get charge
		void set_initial_properties();//set properties
		void set_atomtype_name(string);//sets the atomtype_name
		string get_atomtype_name();//get atomtype_name
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

class Atomcontainer//Atomcontainer class. Has the information on adjacency and bond order
{
	protected:
		vector<Atom*> atoms;//vector of atoms
		vector< vector <int> > Adjacency;//Adjacency list
		vector < vector<int> > BO;//Bond order list //Note 1- single bond, 2- double, 3- triple, 4-aromatic and 5 - nonbonded
		vector<int> NN;//Nearest neighbor (counter) - nonbonded is also counted! 
		vector<float> value; //value
		vector<int> rank;// rank
		vector<int> classRank;//
		int size;// size of the container
		vector<int> label;//label of the atoms
		vector< vector< int > > H_label;//container to store the labels of the Hydrogens
		vector<int> Hydrogens;//number of Hydrogens 
		Ringset Allrings;//Stores the Set of all rings
		vector <int> aromatic_atoms;//vector containing all aromatic atoms
		void EvaluateInitialValue();//evaluates the initial value of the atoms
		void EvaluateInitialValue(int);
		void EvaluateValue();//evaluates the value of the atoms 
		void evaluate_rank();// evaluates the rank of the atoms 
		int distinctRankCount() const;//outputs the number of unique ranks (classes)
		void reverse_ranks();//reverses the order of the ranks
		void sort_adjacency();//sorts the adjacency list based on rank
		void push_to_first(int,int);//a useful manipulating function to push an element in the adjacency list to the first
		void erase(int);//erasing an atom from the atomcontainer
		Atomcontainer form_atomcontainer(vector<int>) const;//create an atomcontainer with a list of atoms.
		bool IsEqualClasses(vector<set<int> >, vector<set<int> >) const;//checks if two atom class classifications are equal
		void breakTies();//if two atoms have the same rank, then this function forcibly breaks the tie
		vector<set<int> > EvaluateAtomClasses();//evaluates the atom classes and puts them into sets of equivalent topological classes
		
		
	public:
		Atomcontainer(int);//constructor
		Atomcontainer();
		~Atomcontainer();
		Atomcontainer(const Atomcontainer &a);
		Atomcontainer& operator=(const Atomcontainer &a);
		void setrank(vector<int>& );//set rank function
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
		string getatomtype(int) const;//getting atomtype of a particular atom
		void merge(const Atomcontainer&);//merge the atomcontainer with another
		int findatomwithlabel(int) const;//find the atom with a given label - returns -1 if none is found! 
		int getlabel(int) const;//get the label of a particular atom
		void setlabel(int,int);//set the label of a particular atom to a specified value
		float getvalue(int) const;//get the value of a particular atom
		void formbond(int,int);//form bond between two specified atoms
		void breakbond(int,int);//break bond between two specified atoms
		void changeBO(int,int,int);//change BO between two specified atoms by a particular value
		void setBO(int,int,int);//Set BO between two specified atoms to a particular value
		void setatomtypename(int,string);//set atomtype name of a particular atom
		void setatomsymbol(int,string);//set atom symbol of a particular atom
		void setInitialAtomproperties(int);//set initial atom properties of specific atom
		void setatomvalency(int);//sets the valency of the atom
		void resetlabel(int);//sets all the label to zero.
		int getHydrogens(int) const;//gives the hydrogen count of the given atom
		void setHydrogens(int, int);//sets the hydroen count of a given atom to second argument's value;
		vector<int> getAromaticAtoms() const;//gives the vector of aromatic atoms of the atomcontainer
		vector<int> getHlabel(int) const;//gives the vector of hydrogen labels of a given atom
		void setHlabel(int, vector<int>);//sets the H label of a given atom with a vector of hydrogen labels.
		int findatomwithHlabel(int) const;//finds the atom that has the given argument as a Hydrogen label - finds -1 if none is found
		void changeHlabel(int, int);//change the hydrogen label given by the first argument to the value given in the second argument
		int getBplusNBatomcount(int) const;//gets the count of both bonded and nonbonded atoms that are neighboring to the given atom
		int getBatomcount(int) const;//gets the bonded atom count of the given atom
		int AromaticBondCount(int) const;//gets the count of the aromatic bonds attached to the atom
		bool HasNBinteractions() const;//checks if there are nonbonded interactions in the atomcontainer
		bool IsAdjacentAtomWithBO(int, string, int) const;//checks if the atom corresponding to the first argument is adjacent to an atom of the given string, with a bond order of the given third argument
		pair<vector<Atomcontainer>, vector< vector<int> > > connectedcomponents() const;//generates a vector of atomcontainers that are essentially the connected components of the parent atomcontainer and also the atomindices of these components
		int getElectronicHashValue(int) const;//returns positive charge ->2^charge; negative -> 3^charge; neutral ->1; radical -> 5; positive radical ->7^value; used in  
		int getTotalElectronicValue() const;//returns product of getElectronicHashValue for all the atoms
		void removeHlabel(int, int);//removes the Hlabel (second argument) from the parent atom (first argument)
		void addHlabel(int,int);//add the Hlabel (second argument) to the parent atom (first argument)
		void addAtom(string, string, int, char, int);//add a new atom -atomtypename, atomsymbol, isotopenumber, elementname, and nature
		int getGroupHash(int) const;//gets the hash value (or invariant) of atom specified that takes into consideration itself and its immediate neighbors
		int getNNElecHash(int) const;//gets nearest neighbors' electronic info as a hash value		
		//int getNNBondCount(int, int) const;//gets the number of bonds adjacent to the specified atom (first argument) of specified BOtype (Second argument), not counting the bonds to the specified atom. For eg (1,2) means finding double bond count of atoms adjacent to atom 1.
		int getNNDoubleBondCount(int) const;//gets the number of double bonds adjacent to the specified atom
		int getNNTripleBondCount(int) const;//gets the number of triple bonds adjacent to the specified atom
		int getNNDoubleTripleBondFactors(int) const;//gets a hash factor to account for double and triple bonds! 
		int RingCountOfAtom(int) const;//finds the number of rings an atom is in	
		map<int,int> getClassesFreqMapNeighboringAtom(int) const;//for specified atom, it creates a map of the different classes and their frequencies
		string AtomCenteredGroupForGA(int) const;//generates the atom-centered group of specified atom needed for group!
		int AtomValueForGA(string, bool) const;
		set<string> GetElements() const;
};

class Molecule: public Atomcontainer//class Molecule derived from Atomcontainer
{
	protected:
		string smilesstring;//contains the SMILES like description of pattern
		void dfsvisit();//function for DFS traversing
		void find_all_rings();//function to find the Set of All Rings (SAR)
		void find_aromatic_rings();//find which of the SAR are aromatic rings 
		void find_allylic_atoms();//finds which of the carbon atom is allylic
		vector <int> aromatic_rings;//vector containing the indices of Allrings that are aromatic
		string MolecularFormula;//stores the molecular formula - in the order Carbon, Hydrogen,Oxygen, Nitrogen, Sulphur, and Phosphorous. Following this will be compositeatomtypes. 
		void update_aromaticity_details();//function that makes appropriate updation upon finding/ noting aromatic atoms
		set<int> allylic_atoms;//vector containing the indices of allylic atoms;
		void EvaluateMF();//evaluates the Molecular Formula of the molecule
		string SetWithoutHydrogenCount(string);//removes from the string, any Hydrogen and number following it within square brackets 
		string SetWithoutSquareBrackets(string);//removes from the string any Square Brackets
		string SetWithHydrogenCount(string, int);//adds into the string the Hydrogen count before "]" if any
		string SetWithSquareBrackets(string);//adds into the string Square Brackets. This assumes ring identifiers are not yet appended to the string
		void UpdateSquareBraces();//updates the atomsymbol to remove unwanted square braces and add in necessary ones
		void ModEqRankUsingPrev(vector<int>&);//modifies the ranks of atoms if two atoms have equal ranks in one iteration of canonical SMILES algorithm but had unequal ranks previously
		vector<int> RankOrder(vector<int>&, vector<int>&) const;
		void findRanksOfAllAtoms();
		
	public:
		Molecule(string, int);//constructor
		Molecule(Atomcontainer&);
		string moleculestring() const;
		string unique_smiles();//generate unique smiles
		string unique_smiles(int);
		void print_smiles();//print smilesstring function
		void set_smiles(string);// set the smilesstring
		void print_rings();//prints all rings
		void Readjustproperties();//readjusts properties for those that have higher oxidation state than normal
		void calculateHydrogens();//calculates the hydrogens for all the atoms
		void PerceiveKekule(Path);//Perceives the bond redistribution for a given aromatic atom path
		void print_aromatic_rings();//prints the indices of Allrings that are aromatic rings
		void remove_Hydrogens();//remove Hydrogens, if any in teh Atomcontainer
		bool isaromaticmolecule() const;//checks if the molecule is an aromatic molecules
		bool isaromaticbond(int,int) const;//checks if the bond is aromatic
		bool iscyclicmolecule() const;//checks if the molecule is cyclic
		bool isringatom(int) const;//checks if the atom is in a ring
		bool isringbond(int,int) const;//checks if bond is in a ring
		bool isallylicatom(int) const;//checks if the atom is allylic
		bool isneutral() const;//checks if the atom is neutral
		bool IsNeighbor(int, int) const;//checks if the second atom is a neighbor of the first
		bool InSameRing(int, int) const;//checks if two atoms are in the same ring
		bool RingWithNBInteractions(Path&) const;//checks if the ring has any nonbonded interactions
		bool ishydrogenic() const;//checks if the molecule is hydrogenic
		bool isparaffinic() const;//checks if the molecule is indeed a paraffin
		bool isolefinic() const;//checks if the molecule is indeed an olefin
		bool isNaphthenic() const; //checks if the molecule is Naphthenic (cyclic nonaromatic hydrocarbons)
		bool ishydrocarbon() const;//checks if the molecule is a hydrocarbon
		int totalcharge() const;//calculates the total charge
		int totalupElectrons() const;//calculates the total unpaired electron count
		int MaxRingsize() const;//Calculates the largest Ring size
		int MaxRingSizeOfAtom(int) const;//calculates the largest ring in which the atom resides
		int SmallestRingSize() const;//calculates the size of the smallest ring		
		int AdjacentBondInRingWithOrder(Path&, int, int ) const;//counts the number of adjacent bonds of a bond type (third arg) present in a ring (first arg) for an atom (second arg)
		int totalDoubleBonds() const;//calculates the total number of double bonds
		int totalTripleBonds() const;//calculates the total number of triple bonds
		int totalAromaticBonds() const; //calculates the total number of aromatic bonds
		bool HasBridgeHead() const;//checks if the molecule has a bridge head
		IntPair Hydrogen_counter(int) const;//return a pair- first being total H count, second being D count
		string GetMF() const;//get molecular formula
		string getMFwithoutAtoms(set<string>) const;
		pair<string,int> checkValencyError() const;//throws a pair <string,int> for error in valency!  
		pair<int,int> calcBranchandDistanceValue() const;//calculates a value used for lumping based on distances between leaves etc.
		bool isIntermediate() const;//checks if a molecule is an intermediate --if molecule has a charge, unpaired electron or has non-bonded interactions
		bool HasCompositeAtoms() const;//checks if the molecule has a composite atom
		int NumberOfAromaticRings() const;
		int NumberOfRings() const;
		int totalAtomsIncludingH() const;
		int totalAtomsOfType(string) const;//counts number of atoms of a given atomtype
		int MolecularWeight() const;
		bool isAtomChiral(int) const;
		int NumberChiralAtoms() const;
		int OptIsomers() const; //gives the number of optical isomers - min 1. 
		int NumberOfLeaves() const;
		
		friend class Patternmatch;//Patternmatch is a friend class
};

class Substructure: public Atomcontainer//class Substructure derived from Atomcontainer
{
	protected:
		string fragmentstring; //contains the SMARTS like description of pattern
		vector <env_set>AtomEnv;//set of env_sets - each atom could have a set of atoms or groups and environmental constraints.
		
		
		vector <vector<int> > atom_flag;
		vector <Triplet> ringbondDefs;
	public:
		Substructure (string,int);//constructor
		void printstring();//prints out the fragment string
		string getstring();//gets the fragmentstring
		
		int isringbondcheck(int,int) const;//checks if a bond is in a ring: -1 implies no check reqd, 0 implies forbid ring bond and 1 implies have ring bond
		friend class Patternmatch;//Patternmatch is a friend class
		friend class rxn_net_gen;	
		friend class Reactiontype;
};


class Patternmatch //Patternmatch class 
{
	protected:
		vector < vector<int> > Matches;//stores all the matches
		vector < vector<int> > rank_matches;//stores the rank of each of the atoms of all the matches
		vector < vector <int> > M;//the Matrix M
		vector <vector <vector <int> > > Ma;//the 3D matrix that stores all intermediate M Matrices
		void find_matches(const Molecule&,const Substructure&, int);//function that finds the matches
		void refine(const Molecule&,const Substructure&);//function for refinement
		int check_zeros();//function that checks if M has an entire row of zeros
		bool checkatomtypematch(string, string);
		bool checkconnectivitymatch(int, int);
		bool checkatomenvironmentmatch(int, int);
		bool checkatomflagmatch(int,int);
		bool IsSameCharacteristrics(int, int);
		int Subsize;//size of Substructure 
		int Molsize;//Size of the Molecule
		
		const Molecule* Mo;
		const Substructure* Su;
		vector<int> Hfactor;
		bool atomsetmatch(int,int);
		vector< int >unique_matcheslist;//list of indices of unique matches
		vector<int> unique_matches_frequency;//frequency of occurences of each
		void InitializeAndSetMa(int);
		void calcHfactor();
		Patternmatch(const Molecule&, const Substructure&,int, int);//constructor - the molecule, the substructure to be found, option - 0 implies find all matches, anything else implies stop at first match. The fourth option is an atomindex of the molecule to which the first atom of the substructure wil be matched.
		//This feature is used specifically for nested SMARTS (atom environment constraints being the presence or absence of entire groups). The index is set to -1 for no specification and a valid index value at other points. 

	public:
		Patternmatch(const Molecule&, const Substructure&,int);//Constructor - the molecule, the substructure to be found, and an option - 0 implies find all matches, anything else implies stop at 1. 
		Patternmatch();
		void list_matches();//lists all the matches in an ordered manner
		void unique_matches();//lists the unique matches - gives matches which are topologically unique - but different unique matches could have the same set of atoms. For example, matching CC in ethane would give two unique matches. To be used for reactions
		int GetDistinctMatches();//lists the distinct matches - here, any two pair of matches of the same pattern should have atleast one distinct atom in each match. To be used for structural constraints
		void print_M();//prints M
		int H_factor(int);//gets the H_factor for a match that has to be multiplied with MatchFrequency to get overall reaction frequency
		int H_factor_unique_match(int);//gets the H_factor for a particular unique match, utility same as above
		void print_Ma(int);//prints out Md which is Ma[d]
		int number_of_matches();//returns the number of matches
		int number_of_unique_matches();//returns the number of unique matches
		vector<int>get_unique_matches(int);//gets the particular unique match
		vector<int> get_match(int);//get a particular match
		int getMatchFrequency(int);//gets the number of times ith unique match (i is the argument) occurs in the molecule
		set<int> AtomsCoveredBySubstr(int);//gets the set of atoms matched by the specified atom of hte substructure from the list of all matches
		bool IsMatchWithinRing(int);//checks if a particular match is completely contained in a ring
		pair<int,int> GetDistinctAllAndRingMatches();//gives a pair of all distinct matches and those that are only within a ring
		friend class Molecule;
};	


typedef bool (*ConstrPtr)(const Molecule &);
typedef bool (*CombinedConstrPtr)(const Molecule &, const Molecule &);

class Reactiontype
{
	protected:
		vector <Substructure> reactant_pattern;//reactant pattern list
		bool IntraMolecularRxnAlso;
		bool IntraMolecularRxnOnly;
		vector <bond_operations> bondchanges;// a vector of bond_operations to store connectivity changes
		map<int,string> mod_atomtype;//a map indicating how an atomtype of an atom with the given label changes
		vector <Substructure> struct_constraints;
		vector<int> FragmentCopyIndex;//-1 for no duplicate, 0 for the first reactant, 1 for the second reactant
		vector< map<int,int> > FragmentCopyLabels;
		list<ConstrPtr> RxnConstraints;
		ConstrPtr ProductConstraints;
		CombinedConstrPtr CombinedConstraint;
		void add_bondchanges(bond_operations);
		int Cost;
		string RuleName;//name assigned to the rule
		int maxSpeciesRank;//the rank of the maximum species that can react in this rule.
		bool isSelfRxnOnly;
		bool checkRateConstant;
		double Min_k_value;
		
	public:
		Reactiontype();
		void add_reactant_pattern(Substructure);//add a Substructure to the reactant patterns
		
		//void add_constraintlist(ConstraintsPointerList);
		void add_reactantconstraint(ConstrPtr);
		void add_productconstraint(ConstrPtr);
		void add_combined_constraint(CombinedConstrPtr);
		void add_mod_atomtype(int,string);
		int get_molecularity();
		void AllowIntraMolecularRxnOnly();
		void AllowIntraMolecularRxnAlso();
		bool isIntraMolecularAlso();
		bool isIntraMolecularOnly();
		void AddReactantCopy(int, map<int,int>);//specifies which reactant is copied and what the atom maps are
		void disconnect_bond(int,int);
		void connect_bond(int,int);
		void increaseBO(int,int, int);
		void decreaseBO(int,int, int);
		void setCost(int);
		int getCost();
		void setRuleName(string);
		string getRuleName();
		int getFragmentCopyIndex(int);
		bool BreaksAromaticity(int);//finds if the atom with the given label breaks aromaticity.
		friend class Reaction;//Reaction is a friend of reactiontype
		friend class rxn_net_gen;//rxn_net_gen is a friend of reactiontype too!
		void setSpeciesRank(int);
		int getSpeciesRank();
		void AllowSelfRxnOnly();
		void setMinRateConst(double); //set's a minimum rate constant value
		double getMinRateConst();
		bool shouldCheckRateConstant();//
};


class generated_rxn
{
	protected:
		deque<string*> reactants;
		deque<string*> products;
		int frequency;
		string reactionMF;
		int RuleIndex;
		bool IsRxnIntramolecular;
		map<int, int> ParentMolecule;
		map<int, vector<int> >AtomsFromOtherReactants;//keeps track of which product has come from second and subsequent reactants
		int countSpecies(deque<string*>&, string*);//counts the number of species in the reactant/product deque
	public:
		generated_rxn();
		void add_reactants(string*);
		void add_products(string*);
		string* get_reactants(int);
		string* get_products(int);
		string reactionstring();
		bool IsReactant(string*);
		bool IsProduct (string*);
		int number_pdcts();
		int number_reactants();
		int get_frequency();
		void set_frequency(int);
		void set_reactionMF(string);
		string get_reactionMF();
		void set_rule(int);
		int get_rule();
		int occurence(string*); //gets the absolute value of the occurrence -> abs(stoich in reactant - stoich in products)
		void setAtomsFromOthers(int, vector<int> );
		void setParentMolecule(map<int, int>);
		vector<int> getAtomsFromOthers(int);
		int getParentMolecule(int);
		int getParentMolecule(string*);
		bool isSameNetReaction(generated_rxn&);
		bool isReverseReaction(generated_rxn&);
		bool isReactionIntramolecular();
		void setIntramolecularity(bool);
		vector<string*> getDaughterMolecules(string*);//gets all the molecules in the product having the given molecule as closer (parent).
		int getReactantStoich(string*);//gets stoich in reactant 
		int getProductStoich(string*);//gets stoich in product
		int NetOccurence(string*);//gives the actual value - net negative means more reactants - net positive means more products. 
		int getdeltaN();//gives the change in the number of moles. (reactants - products)
		int NetPatternDiff(string); //gives the net difference of the pattern occurrences in the reactants and products - negative means more in reactants! 
		int netMassDiff();//calculates if there is a mass difference between reactants and products

		bool ProductsWithSameRadicalType(); //checks if there are two radicals in the product of the same atomtype
		int IdenticalProductsFactor();//calculates the product of number of different products that are identical to each other. e.g. if products are P, P,and P - the factor is 3; if P P Q Q - factor is 4, but if P P Q, then it's still 1. 
		
		
		
};
	
class Reaction
{
	protected:
		vector<Molecule> M;
		Reactiontype Rt;
		vector<generated_rxn> gen_rxns;
		vector <Patternmatch> reactpattlist;
		bool GenerateIntraMolecularRxnsAlso;
		bool GenerateIntraMolecularRxnsOnly;
		bool OnlyOneReactantUsed;
		void generate_reactions();
		void findmatch(int);		
		void addproducts(Atomcontainer&, vector<Molecule>&, vector<pair<string, int> >&, int);
		bool perform_changes(Atomcontainer&,vector<pair<string,int> >&);//performs the atomtype modifications and bond order changes. Returns true if changes can be performed, false if it runs into a problem
		void GenerateRxnWithFragmentCopies(Atomcontainer, int, int);

	public:
		Reaction (vector<Molecule> &, Reactiontype &, vector<Patternmatch>& );
		generated_rxn get_generated_rxns(int);
		int number_rxns_generated();
}; 


class LumpInfo
{
	protected:
		int size;
		int HydCount;
		int DoubleBondCount;
		int TripleBondCount;
		int AromaticBondCount;
		int MaxRingSize;
		int RingCount;
		bool isCyclic;
		int BranchValue;
		int ringleavesDistance;
		string* MoleculeString;
		int rank;//rank of the molecule;
		int leaves;
		int PONAcharacteristic;//Paraffin, Olefin, Naphthenics and aromatics characteristics - To lump all olefins, paraffins and hydrocarbon aromatics of one size together.
			// - 1 don't really care or it is none of P, O, N, or A; 0-it is paraffinic, 1- it is olefinic,  2- cycloalkanes and cycloalkenes; 3-it is hydrocarbon aromatic with alkyl or no substituents; 4- hydrocaron aromatic with alkenyl substitutents 
		int PONAElectronicValue;//stores the Electronic value of POA. 
	public:
		LumpInfo(int,int,bool,int, int, string*);
		void setMoleculeString(string*);
		string getMoleculeString();
		string* getMolStringPtr();
		int getSize();
		int getHydrogens();
		int getBranchValue();
		int getringDistanceValue();
		void setBranchValue(int);
		void setringDistanceValue(int);
		int getleaves();
		void setleaves(int);
		void setPONAcharacteristic(int);
		int getPONAcharacteristic();
		int getPONAElectronicValue();
		void setPONAElectronicValue(int);
		void setRank(int);
		int getRank();
		void setBondCounts(int,int,int);//provide double, triple, aromatic
		int getDoubleBondCount();
		int getTripleBondCount();
		int getAromaticBondCount();
		int getMaxRingSize();
		int getRingCount();
		void setMaxRingSize(int);
		void setRingCount(int);
		int getScore();//This returns a score calculated as 3^Doublebonds*5^Triplebonds*7^AromaticBonds
};

class LumpingStrategy
{
	protected:
		bool toLump;
		int chainParameter;//-1 - no lumping of chains, 0-branches closest together, 1-branches farthest apart
		int ringParameter;//-1 - no lumping of rings, 0- substituents closest together, 1-substituents farthest apart
		int paraffinParameter;//-1 no paraffin lumping; 0 - all paraffins lumped to least branched, 1 - all paraffins lumped to most branched
		int olefinParameter;//-1 no olefin lumping; 0 - all olefins lumped to least branched, 1- all olefins lumped to most branched
		int aromaticsParameter;//-1 no aromatics lumping; 0 - all aromatics lumped to least branched, 1- all aromatics lumped to most branched
		int naphthenicsParameter;//-1 no cycloalkanes and cycloalkenes lumping; 1 - all lumped to least branched; 1- all lumped to most branched 
		vector<int> MoreLumpingParameter;// -1- all lumped to least branched; 1 - all lumped to most branched
		
		bool HasSetFunctionalLumpingConstraints, HasSetParaffinConstraints, HasSetOlefinConstraints, HasSetNaphthenicsConstraints, HasSetAromaticsConstraints, HasSetMoreLumpingConstraints;
				
		//molecule constraint pointers that describe structural constraints on what molecules can be lumped for initial functional lumping and subsequent lumping of paraffins, olefins, naphthenics, aromatics, and any other additional lumping
		ConstrPtr FunctionalLumpingConstrPtr;
		ConstrPtr ParaffinConstrPtr;
		ConstrPtr OlefinConstrPtr;
		ConstrPtr NaphthenicsConstrPtr;
		ConstrPtr AromaticsConstrPtr;
		vector<ConstrPtr> MoreLumpingConstrPtr;

	public:
		LumpingStrategy(bool);
		LumpingStrategy();
		void setParameters(int,int, int, int, int, int);
		bool shoudLump();
		int getchainParameter();
		int getringParameter();
		int getParaffinParameter();
		int getOlefinParameter();
		int getAromaticsParameter();
		int getNaphthenicsParameter();
		vector<int> getMoreLumpingParameter();
		void setFunctionalLumpingConstraints(ConstrPtr);
		void setParaffinConstraints(ConstrPtr);
		void setOlefinConstraints(ConstrPtr);
		void setNaphthenicsConstraints(ConstrPtr);
		void setAromaticsConstraints(ConstrPtr);
		void setMoreLumpingConstraints(ConstrPtr,int);//the integer specifies what kind of representative is desired -- see above
		ConstrPtr getParaffinConstraints();
		ConstrPtr getOlefinConstraints();
		ConstrPtr getNaphthenicsConstraints();
		ConstrPtr getAromaticsConstraints(); 
		ConstrPtr getFunctionalLumpingConstraints();
		vector<ConstrPtr> getMoreLumpingConstraints();
		bool isThereFunctionalConstraints();
		bool isThereParaffinConstraints();
		bool isThereOlefinConstraints();
		bool isThereNaphthenicsConstraints();
		bool isThereAromaticsConstraints();
		bool isThereMoreLumpingConstraints();

};

class LumpedReaction
{
	protected:
		vector<int> reactantsLumpSet;
		vector<int> productsLumpSet;
		int ruleIndex;
	public:
		LumpedReaction(vector<int>, vector<int>, int);
		vector<int> getReactantLumps() const; 
		vector<int> getProductLumps() const;
		int getRule() const;
		
};

struct LumpedReactionCompare{
  bool operator() (LumpedReaction a, LumpedReaction b) const
  {
	  if (a.getReactantLumps()<b.getReactantLumps()) return true;
	  else if (a.getReactantLumps()>b.getReactantLumps()) return false;
	  else
	  {
		  if (a.getProductLumps()< b.getProductLumps()) return true;
		  else if (a.getProductLumps()>b.getProductLumps()) return false;
		  else return (a.getRule()<b.getRule());
	  }
	  
  }

};

enum SiteType {  Homogeneous, Heterogeneous};
enum ThermoType {EnthalpyType, EntropyType, CpType, FreeEnergyType, logPType};
class ThermoValues
{
	protected:
		map<int, double> EnthalpyFormation; // the ints refer to indices in AvailableTemps
		map<int, double> EntropyFormation;
		map<int, double> FreeEnergyFormation;
		map<int, double> CpMolecule;
		map<int, double> logP;//incomplete
		vector<double> AvailableTempsForEnthalpy;
		vector<double> AvailableTempsForEntropy;
		vector<double> AvailableTempsForFreeEnergy;
		vector<double> AvailableTempsForCp;
		vector<double> AvailableTempsForlogP;
	public:
		double getThermo (ThermoType, double);
		void setThermo(ThermoType, double, double);
		int isTempAvailable(ThermoType,double); //returns the index of the vector AvailableTemps; -1 if not available!
};

class MoreAdditionalLumpingInfo //Note this class keeps track of all additional lumping stuff the language allows for (that is in addition to PONA)
{
	public:
		int Size;
		int HydCount;
		int DoubleBonds;
		int TripleBonds;
		int AromaticBonds;
		int NumberOfRings;
		int MaxRingSize;
		int WhichMoreLumpingDcl;//keeps the index of which of the possibly many declarations of more lumping
};

struct MoreLumpingCompare{
  bool operator() (MoreAdditionalLumpingInfo a, MoreAdditionalLumpingInfo b) const
  {
	  if (a.WhichMoreLumpingDcl < b.WhichMoreLumpingDcl) return true;
	  else if (a.WhichMoreLumpingDcl > b.WhichMoreLumpingDcl) return false;
	  else
	  {
		  if (a.Size<b.Size)return true;
		  else if (a.Size>b.Size) return false;
		  else
		  {
			 if (a.HydCount<b.HydCount) return true;
			 else if (a.HydCount>b.HydCount) return false;
			 else
			 {
				 if (a.DoubleBonds <b.DoubleBonds) return true;
				 else if (a.DoubleBonds > b.DoubleBonds) return false;
				 else
				 {
					 if (a.TripleBonds < b.TripleBonds) return true;
					 else if (a.TripleBonds > b.TripleBonds) return false;
					 else
					 {
						 if (a.AromaticBonds < b.AromaticBonds) return true;
						 else if (a.AromaticBonds > b.AromaticBonds) return false;
						 else 
						 {
							 if (a.NumberOfRings < b.NumberOfRings) return true;
							 else if (a.NumberOfRings > b.NumberOfRings) return false;
							 else return(a.MaxRingSize < b.MaxRingSize); 
						 }
					 }
				}
			 }
		  }
	  }
  }
						 

					 
					 

};

typedef bool (*KineticParamPtr)(vector<Molecule>&, vector<Molecule>&, double&, double&, double&, bool&, bool&, bool&, double&,double&,bool&,double);


typedef pair<unsigned int, unsigned int> UnsignedIntPair;
typedef pair<pair<int,int>,int> PONALumpingInfo; //the inner pair holds size and hydrogen count, while the outer pair's second value is PONAElectronic value

class rxn_net_gen
{
	protected:
		multimap<int, Molecule*> unprocessedmol; //molecules yet to be processed
		vector<Molecule*> processedmol; //molecules already processed
		set<string*, classcomp> InitialReactants; //the initial reactants input by the user
		list<string>* inputStrings; //list of strings corresponding to the input specified by the user (note these are not in canonical SMILES necessarily)
		vector<pair<string, SiteType> > CompositeSites; //all composite atoms that have been specified as sites
		vector<string> CompositeAtoms; // all composite atoms 
		multimap<int, int> ReactionsMap; //a map of rules to all reactions (referred to by their reaction number NOT index, so starts from one) belonging to that rule.
		vector<generated_rxn> AllReactions; // an array of all generated reactions obtained through network generation
		vector<generated_rxn> ReconstructedReactions; // an array of all reactions generated through the process of reconstruction (populated only when reconstructed from the original network - kicked in ONLY when lumping is OFF)
		vector<generated_rxn> LumpedReconstructedNetwork; //final array of all lumped and reconstructed network.
		vector<LumpInfo> MolLumps; //An array of different lumps
		map<LumpedReaction, int, LumpedReactionCompare> LumpedReactionMap; //A map of each lumped reaction and it's frequency: consider moving it in to LumpedReaction() function
		multimap<int,int>LumpedReconstructedReactionsMap; //similar to reactions map but for LumpedReconstructed reactions
		
		map<PONALumpingInfo,int> ParaffinSizeLumpMap;//map with key as pair of Size-Hydrogen count pair and PONAElectronicValue; value = lump;
		map<PONALumpingInfo,int> OlefinSizeLumpMap;
		map<PONALumpingInfo,int> AlkylAromaticsBranchLumpMap;
		map<PONALumpingInfo,int> AlkenylAromaticsBranchLumpMap;
		map<PONALumpingInfo,int> NapthenicsSizeLumpMap;
		map<MoreAdditionalLumpingInfo,int,MoreLumpingCompare> MoreAdditionalLumpingSizeLumpMap; //for all additional specifications to do more lumping --in addition to PONA
		map<int,int> AdditionalLumpMap; //a map of additional lumps that maps the MolLumps lumps to other MolLumps lumps
		set<int> ToUndergoAdditionalLumping; //stores a list of MolLumps indices that need to undergo additional lumping


		map<string*,int, classcomp> AllMolecules;
		map<string*,int, classcomp> AllLumpedMolecules;//this keeps the final map of strings representative of the lump and it's rank - much like AllMolecules, except this is for the final set of lumps! 
		
		set<string*,classcomp> Intermediates;
		multimap<UnsignedIntPair, int > LumpHashMap;// map of hash value and the corresponding lump - given by the index of the lump in MolLumps
		map<string*,int, classcomp> MolLumpMap; // map of molecule and lump - given by the index in MolLumps

		multimap<string*,int> MolReactantMap;//says which reactions have the molecule as a reactant
		multimap<string*,int> MolProductMap;//says which reactions have the molecule as a product

		multimap<string*,int> MolLumpReactantMap;//says which lumpedreconstructed reaction has molecule as a reactant
		multimap<string*,int> MolLumpProductMap; // says which lumpedreconstructed reaction has a molecule as a product
		
		map<int,int> AtomtypeIndex;
		vector<Reactiontype> Rtlist;
		map<int, int> AtomtypingMap;//first int stores the value of the atomtype and the second stores its index corresponding to a number in the prime array. 
		Patternmatch checkmatch(vector<Molecule>&,int&, int);
		Patternmatch check_reactant_pattern(Molecule&, Substructure&);//checks for presence of reactant pattern instance
		ConstrPtr GlobalConstraints;//global constraints
		LumpingStrategy LumpStrat;
		bool shouldCalcThermo;
		map<string*, ThermoValues, classcomp> AllMolThermo;
		vector<int> RxnsWithInterchangeableReactants;//reactions where the reactants can both swap places but will lead to the same reaction. 
		
		double Temperature;
		bool isunique(string);
		bool globalconstraintcheck(Molecule&);
		bool check_combined_match(vector<Molecule>&,int&);
		bool check_product_constraints(Molecule&,int);
		//pair<unsigned int, unsigned int> CreateHashValue(const Molecule&);
		int GetAtomValue(string);
		void LumpMolecule(const Molecule &, string*);
		void LumpParaffins();
		void LumpOlefins();
		void LumpHydrocarbonAromatics();
		void LumpReactions();
		void LumpNaphthenics();
		void SetPONALumpsMap();
		LumpedReaction GenerateLumpedReaction(generated_rxn);
		generated_rxn GetRxnFromLumpedReaction(const LumpedReaction&);
		void PrintLumps();
		void printFinalLumpSpeciesAndRxns();

		string getRxnStringFromLumpedRxn(LumpedReaction);
		int findLumpofMolecule(string*);
		bool isSiteIntermediate(string);//this checks if a species (that is not an initial reactant) is an intermediate
		bool containsSiteAtom(string); // this checks if a species includes a composite site. -- a superset of isIntermediate
		double calculateThermoOfRxn(ThermoType, generated_rxn&, double);
		void add_unique_molecules_reactions(Reaction&, int, bool);//add the new molecule and reactions into the container of molecules and reactions. the int refers to the Reactiontype index. the bool says if we need to keep track of reactions where the popped molecule can participate as either reactant.
		void UpdateMolAsReactantsProductsMap(generated_rxn&);
		void GetHighestRankInfoForRxn(generated_rxn&, int&, int&, bool&);
		void GetProductParentsInfoForRxn(generated_rxn&, map<int,int>&, int, bool, int);

		bool setInitialReactants();
		bool MolValencyErrorStatements(pair<string, int>&);
		void DoPONALumping();
		void DoMoreLumping();
		void printOutputInfo();
		void calculateThermoValues();
		void GenerateMonoMolecularRxns(int, vector<Molecule>&, vector<Patternmatch>&);
		void GenerateBimolecularRxns(int, vector<Molecule>&, vector<Patternmatch>&, vector<vector<pair<Molecule*,Patternmatch> > >&, int&);
		void ReconstructForOriginalNetwork(int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, string);
		void ReconstructForLumpedNetwork (int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, string);
		generated_rxn GetOneRxnFromTwo(generated_rxn&, generated_rxn&);
		void SetAllLumpedMoleculesInfo(generated_rxn& r1);
		set<int> RulesRemovedByReconstruction;
		bool calcMolThermo(string*,ThermoType,double,double&);
		vector<KineticParamPtr>* kinetics;

		int getDeltaNGasPhase(generated_rxn&);
		int firstGasSpecies(int,int);//first argument is reaction index, second is 0,1 for reactant or product

		/* -------for getting network from files ------*/

		bool ReadReactionFile(const char*, map<int,string*>&);
		bool ReadSpeciesFile(const char*, map<int,string*>&);
		bool TraverseNetworkToGetRanks();
		void LumpAllMolecules();

		bool InterpretRxnsAndUpdate(vector<string>&, int, map<int,string*>&,int);

		/*-------------------------*/

		
		/*------- GAMS --------*/
		bool DoSimultaneousRxns;
		vector< set<int> > SimultaneousRxns;
		map<int, vector<int> > RxnsToBeSummedUpMap;
		set<int> SimultRxnsGenerated;
		set< set<string> > reactantpairsForBiMol;//pairs of reactants for the case of bimolecular reactions
		int FindRxnWithReactantsAndRule(const set<string>&, int);
		set<int> FindReactionsWithReactantsAndRule(const set<string>&, int);
		void findRxnsOfSimultRules(vector<pair<int,int> >&);
		void findTransformRxns(vector<pair<int,int> >&);
		void IdentifySimultaneousBimolRxns();
		map<string*, string, classcomp> SMILESGAMSSpeciesMap;
		multimap<string, string> DenticityInfoMap;
		int CalculateNetStoichDifferenceWithDenticity(generated_rxn&, string);

		/*-------- end ----------*/

	public:
		rxn_net_gen(list<string>&, vector<Reactiontype>&, ConstrPtr, LumpingStrategy&, vector<string>&, vector<string>&, bool);
		//Incomplete -- the above constructor takes in a list of strings (input SMILES strings), vector of all the reaction rules, a constraint pointer for global constraints, Lumping strategy and a list of composite atoms-sites
		
		rxn_net_gen(){DoSimultaneousRxns = false; Temperature =298.0;}
		void AddInitialReactants(list<string>&);//add inital reactants
		void AddReactionRules(vector<Reactiontype>&);//add reaction rule vector
		void AddGlobalConstraints(ConstrPtr);//add global constraints pointer
		void AddLumpingStrategy(LumpingStrategy&); // add lumping strategy info
		void AddCompositeAtoms(vector<string>&);// add composite atoms
		void AddCompositeSites(vector<pair<string, SiteType> >&);// add composite atoms that are specifically sites
		void SetCalcThermo(bool);//set to true/ false if thermo has to be calculated
		void GenerateNetwork();//begin generating the network;
		void SetAllLumpedRxns();//sets lumped reactions into a new vector of reactions -- NOTE that this SHOULD be called after network generation!
		void GetNetworkFromFile(const char*, const char*);
		rxn_net_gen(rxn_net_gen*, LumpingStrategy&); //TODO
		void generateStoichMatrix();
		void ReconstructReactions(int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, string);
		void print_rxnlist();
		generated_rxn getReaction(int);
		int findSpeciesRank(string);
		void GenerateInfoForAthena();
		void SetTemp(double);
		void StoreRxnsAndSpecies(const char*, const char*, bool);
		void setKinetics(vector<KineticParamPtr>&);

		pair<unsigned int, unsigned int> CreateHashValue(const Molecule&);
		


		
		//for GAMS
		void CalculateParametersFileForGAMS();
		void GatherSimultaneousAndJointReactionsForGAMS(vector<pair<int,int> >&, vector<pair<int,int> >& );
		void setSimultaneousRxns(bool value){DoSimultaneousRxns = value;}
		void PrepareSpeciesInfoForGAMS();
		void CheckDensityInfo(multimap<string,string>&);


		//---------------------------------------------------------


		
		friend class Pathways;	
		friend class Mechanisms;
		friend class MoleculeQuery;
		friend class ReactionQuery;
		friend class KineticModel;
		friend class KineticsInfo;
		friend class CHEMKinFiles;
		friend class GAMSFiles;

};




class RxnInfo
{
	protected:
		generated_rxn* rxn;
		double Enthalpy;
		double Entropy;
		double FreeEnergy;
		double ActEnergy;
	public:
		RxnInfo (generated_rxn*, double, double, double, double);
		double ActE();
};

typedef bool (*RxnConstraintPtr)(RxnInfo&);


class RxnPathway 
{
	protected:
		//vector<int> * rxn_indices;
		vector<generated_rxn> * rxn_index;
		vector<double>* RxnEnthalpy;
		vector<double>* RxnEntropy;
		vector<double>* RxnFreeEnergy;
		vector<double>* RxnActE;

	public:
		//RxnPathway(vector<generated_rxn>*, vector<double>*, vector<double>*, vector<double>*, vector<double>*);
		RxnPathway(vector<generated_rxn>*);
		void setRxnActE(vector<double>*);
		int numRule(int);  // index of rule to number of occurences
		int numRuleWithParticipants(int, ConstrPtr, int);  // index of rule, molecule constraint, and an integer indicator to indicate constraints is on reactants (0) or products (1) to number of occurences
		int numRuleWithParticipants(int, CombinedConstrPtr, int);  // index of rule, combined molecule constraint, and an integer indicator as above to number of occurences
		int numMolecule(ConstrPtr,int); // how often this molecule appears in this pathway as reactant or product.
		int numRuleIntramolecular(int);//how many intramolecular reactions of a given type
		int numIntramolecular();//how many intramolecular reactions in total!
		int numSatisfying(int, RxnConstraintPtr); //number of instances of reactions of a given rule satifying some reaction constraint parameters
			
};



typedef bool (*PathwayConstrPtr)(RxnPathway&);

class PathwayConstraints
{
	protected:
		int MaxLength;
		int MaxCost;
		int MinLength;
		int MinCost;
		int MaxRelativeLength;//difference over and above the shortest path
		//vector<PCRule> RuleConstr;
		//vector<PCRuleMol> RulePlusMolConstr;
		//vector<PCRuleConstr> RulePlusFeatureConstr;
		PathwayConstrPtr PathConstraints; 
				
	public:
		PathwayConstraints();
		void AddPathwayConstrPtr(PathwayConstrPtr);
		//void AddRuleConstr(PCRule);
		//void AddRulePlusMolConstr(PCRuleMol);
		//void AddRulePlusFeatuerConstr(PCRuleConstr);
		void SetMaxLength(int);
		void SetMaxRelativeLength(int);
		int getMaxRelativeLength();
		void SetMaxCost(int);
		int getMaxCost();
		void SetMinLengthAndCost(int,int);//min length, min cost
		int getMinCost();
		int getMinLength();
		int getMaxLength();
		bool IsSatisfied(RxnPathway&);
};



class Pathways
{
	protected:
		const char* filename;
		rxn_net_gen* Network;
		vector<double> ReactionRates;
		multimap<string*, vector<int> > FindAllPathways(string*, int);//returns the pathway reaction indices for input molecule string pointer and max path length
		string* CloserMolecule(int, string*);
		PathwayConstraints Constr;
		bool HasCycles(vector<string*>); 
		bool isReverseAlreadyPresent(vector<int>, int);
		bool similarRulePresent(multimap<int,string*>*, int, string*);
		bool DistinctNature; //set to true to find pathways that are distinct in terms of the number of different elemetary steps
		vector<KineticParamPtr>* KineticFunctions;
		bool calculateActEnergy(vector<int>&, vector<double>&);//setting kinetics for each reaction
		bool canCalculateActE;
		void DominantReaction(vector<int>&, vector<double>&);
		void DominantReactionClass(vector<int>&, vector<double>&);
		vector<int> MajorDominantRxns(vector<int>&, vector<double>&);
		multimap<string*, vector<int> > FindAllDominantPathways(string*, int);
		void PopulateDominantProductMap(string*);
		multimap<string*, int> MolDominantProductMap; //says which reactions have the dominant rates of formation for a product
	public:
		Pathways(rxn_net_gen*,PathwayConstraints&, const char*);
		void SetConstraints(PathwayConstraints);
		void GeneratePathways();//pathways to products that are final
		//void GeneratePathways(vector<string>);//pathways to products specified in the argument
		void QueryPlusGenerate(ConstrPtr);//pathways to products that satisfy constraints
		int CalculateCost(vector<int>);
		void setDistinctNature(bool);
		void AddKineticFunctions(vector<KineticParamPtr>&); // add kinetics info
		pair<int,map<int,int> > CalculatePathwayValue(vector<int>);
		void SetRateVector(vector<double>);
		void getDominantPrevStep(string);

		
};

class MoleculeQuery
{
	protected:
		const char* filename;
		rxn_net_gen* Network;
		vector<string*> PassedMolecules;
		void GenerateParameterFileForGAMS();
	public:
		MoleculeQuery(rxn_net_gen*, ConstrPtr, const char*, bool, bool);
		
};

class ReactionQuery
{
	protected:
		const char* filename;
		rxn_net_gen* Network;
		bool checkRxn(generated_rxn, PathwayConstrPtr&);
	public:
		ReactionQuery (rxn_net_gen*, PathwayConstrPtr, const char*,int);
};

class partialMechanism
{
	protected:
		vector < pair<generated_rxn,int> >* rxn_index;	
		map<string*,int> stoichiometry;
		vector<double>* RxnActE;

	public:
		partialMechanism(vector < pair<generated_rxn,int> >*);
		bool isContained(partialMechanism&);
		bool isNetZero();
		string printformula();
		string printSMILES();
		map<string*, int> getStoichiometry();
		int getFreq(string*);
		bool isEqual(partialMechanism*);//is equal in terms of overall stoichiometry!
		int numRule(int);  // index of rule to get number of occurences
		int numRuleWithParticipants(int, ConstrPtr, int);  // index of rule, molecule constraint, and an integer indicator to indicate constraints is on reactants (0) or products (1) to get number of occurences
		int numRuleWithParticipants(int, CombinedConstrPtr, int);  // index of rule, combined molecule constraint, and an integer indicator as above to get number of occurences
		int numMolecule(ConstrPtr); // how often this molecule appears in this mechanism.
		int numMoleculeOverall(ConstrPtr);//similar to numMolecule, applied to the overall rxn, the second argument says reactant (0) or product (1)
		set<string> getProducts(); //gets all the products of this partial mechanism 
		void MechinfoForAthena();
		int numSatisfying(int, RxnConstraintPtr); //number of instances of reactions of a given rule satifying some reaction constraint parameters
		void setRxnActE(vector<double>*);

		friend class Mechanisms;
};




class OverallMechanism
{
	protected:
		vector<pair<string*, int> > MoleculeDirectMechPair;/*stores a pair of which molecule's direct mechanism is part of the overall mechanism and the integer index corresponding to the direct mechanism stored in 
															MoleculeMechanismsMap of class Mechanisms*/
		vector<int> Frequency;//stores how many of each direct mechanism is to be counted in the overall mechanism	
	public:
		OverallMechanism();
		void insertDirectMech(string*,int,int);
		pair<string*,int> getDirectMechanism(int);
		int NumberOfDirectMechs();
		int getFrequency(int);
		map<string*, pair<int,int> > generateStringPairMap();
		void removeLastDirectMech();
		void printMechDetails();
		
};

typedef bool (*DirectMechConstrPtr)(partialMechanism&);

typedef bool (*CompleteMechConstrPtr)(partialMechanism&);

class DirectMechConstraints
{
	protected:
		int MaxLength;
		int MaxCost;
		int MinLength;
		int MinCost;
		DirectMechConstrPtr DMConstraints;
	public:
		DirectMechConstraints();
		void AddDMConstraints(DirectMechConstrPtr);
		void SetMaxLength(int);
		void SetMaxCost(int);
		void SetMinLengthAndCost(int,int);
		bool isSatisfied(partialMechanism&);
		int getMaxLength();
		int getMaxCost();
		int getMinLength();
		int getMinCost();

};


class CompleteMechConstraints
{
	protected:
		int MaxLength;
		int MaxCycleCount;
		int MaxCost;
		int MinLength;
		int MinCycleCount;
		int MinCost;
		CompleteMechConstrPtr MechConstraints;
	public:
		CompleteMechConstraints();
		void AddMechConstraints(CompleteMechConstrPtr);
		void SetMaxLength(int);
		void SetMaxCost(int);
		void SetMinLengthAndCost(int, int);
		void SetMaxCycleCount(int);
		void SetMinCycleCount(int);
		bool isSatisfied(partialMechanism&);
		int getMaxLength();
		int getMaxCost();
		int getMinLength();
		int getMinCost();
		int getMaxCycleCount();
		int getMinCycleCount();
		
};


class OverallRxn
{
	protected:
		map<string*,int> SpeciesStoichio;
	public:
		OverallRxn(partialMechanism*);
		int getFreq(string*);
};

class Mechanisms
{
	protected:
		rxn_net_gen* Network;
		const char* filename;
		map<string*, vector< map< int, int> > > MoleculeMechanismsMap;
		vector< map<int,int> > FindDirectMechanisms(string*, string*,bool);
		vector<OverallMechanism> AllMechanisms;
		set< map<string*, pair<int,int> > > UniqueMechs;
		bool DistinctDirectMechs;
		DirectMechConstraints DirectConstr;
		CompleteMechConstraints CompleteConstr;
		map<int,int> getMechanismMap(vector <pair<int, int> >&);
		bool willFormInternalCycles(generated_rxn, vector< pair<generated_rxn, int> >&);
		vector<KineticParamPtr>* KineticFunctions;
		bool calculateActEnergy(vector<pair<int,int> >&, vector<double>&);//setting kinetics for each reaction
		bool canCalculateActE;
		map<int,double> RxnActEnergyMap;//stores the activation barrier of each reaction -note that this reaction index starts from 1!!!. 

	public:
		Mechanisms (rxn_net_gen*, const char*);
		void GenerateMechanisms(ConstrPtr );
		void GenerateDirectMechanismsOlny(ConstrPtr);
		bool isComplete(partialMechanism*);
		void SetDistinctDirectMech(bool);
		bool ContainsInitialReactantsAsProducts(partialMechanism*);
		bool FoundOverallMechanism(partialMechanism*);
		int GenerateOverallMechanisms(string*);
		int getTotalCost(vector< pair<generated_rxn, int> > *);
		map<string*,int> getIntermediates(partialMechanism*);
		bool isReverseAlreadyPresent(vector<pair<generated_rxn,int> > *, generated_rxn *);
		void SetDirectMechConstraints(DirectMechConstraints);
		void SetCompleteMechConstraints(CompleteMechConstraints);
		void AddKineticFunctions(vector<KineticParamPtr>&); // add kinetics info
};

class Surfqsar
{
	protected:
		double KHZeroCalc();
		Molecule* mol;
		int NumEOUnits;
		double AICSecond();
		double RelativeOxygensCount();
		double KierShapeThird();
		double ABInfoContentZero();
		double StructInfoContentOne();
		
	public:
		Surfqsar(Molecule*, int);
		double calculateCMC();
		double CMCerror();
		double calculateCP();
		double CPerror();
		double calculateHLB();
		double calculateSurfTens();
};


class GATempDiscretePts
{
	protected:
		static vector<double> TempPts;
    public:
		static vector<double> getTemperaturePoints();
		static void InsertTempPoints(double);
};

class GAdata
{
	protected:
		double Enthalpy;
		double Entropy;
		double logP;
		map<int, double> CpMap;
		vector<double> Cp;//assumed to be start from 298
		vector<double> A;
		vector<double> B;
		//double getCpAtEndPoints(double);
		void CalculateAandBForEachRange();
		double getTempCorrection(ThermoType, double) const;
		double getCpAtArbitraryTemp(int) const;

	public:
		GAdata(vector<double>);
		GAdata(){}
		double getEnthalpy();
		double getEntropy();
		double getCpAt(int) const;
		void setEnthalpy(double);
		void setEntropy(double);
		void setCp(int, double);//int is for Temp and double is the value
		void setlogP(double);
		void PrepareGA();// calculates A and B for each range.
		bool getdata(ThermoType, double, double&) const;
		
};

typedef map<string, GAdata> GAMap; //the actual map that stores the groups and its additivity value
typedef map < pair<int,int>, GAMap > HashGAMap;// this maps HashPair with its appropriate GAMAP --> the HashPair keeps track of the num of double and triple bonds of the neighbors and if the neighboring atoms have special electronic features  (positive/ negative/radical). 

enum GADataType { AdditivityType, CorrectionsType};

typedef bool (*ThermoCorrelationPtr)(const Molecule &, double&, double&, double&, double&);

class ThermoGA
{
	protected:
		static map < string, HashGAMap* > AtomtypeGAMap;//this is the big map that keeps track of all groups starting with a given atomtype.
		static bool IncrementAdded(const Molecule&, int, const GAMap&, map<int, pair<int,int> >&, set<int>&, ThermoType, double&, double);
		static multimap <string, pair<ConstrPtr, GAdata> > Corrections;
		static multimap <string, pair<ConstrPtr, double> > BECorrections;
		static vector<pair<ConstrPtr, GAdata> > MolecularCorrections;
		static map <string, HashGAMap* > logPGAMap; //this keeps track of all logP groups
		static double getCorrectionsValue(const Molecule&, ThermoType, double);
		static bool AHashPairExists(const Molecule&, string, int, pair<int,int>&);// constructs reasonable hashpairs and checks if they exist in HashGAMap.
		static bool calculateThermoProperty(const Molecule&, double&, ThermoType, map < string, HashGAMap* >&, double, bool);
		static bool CalculateIncrement(const Molecule&, string, int, map<int, pair<int,int> >&, set<int>&, ThermoType, map < string, HashGAMap* >&, double&, double);
		static double calculateTempCorrections(double, ThermoType);
		static double getSymmetryCorrection(const Molecule&, ThermoType);
		static void AddAdditivityData(string&, GAdata&,map < string, HashGAMap* >&);
		
		

	public:
		static void AddGA(string, GAdata);
		static void AddCorrections(string, ConstrPtr, GAdata);
		static void AddBECorrections(string, ConstrPtr, double);
		static void AddMolecularCorrections(ConstrPtr,GAdata);
		static bool calculateDeltaH(const Molecule&, bool, double, double&);//calculate deltaH formation of a molecule, true or false as to whether or not the molecule is site intermediate, at a given temperature, and store it in the fourth argument passed by reference
		static bool calculateDeltaS(const Molecule&, bool, double, double&);//calculate deltaS formation of a molecule, true or false as to whether or not the molecule is site intermediate, at a given temperature, and store it in the fourth argument passed by reference
		static bool calculateDeltaG(const Molecule&, bool, double, double&);
		static double calculateBE (const Molecule&);
		static bool calculateCp(const Molecule&, bool, double, double&);
		static bool ReadInputsFromFile(const char*,GADataType);//reads from a file the GA and corrections input.
		static void AddlogPGA(string, GAdata);
		static bool calculateLogP(const Molecule&, double&);//calculate logP of a molecule
		static ThermoCorrelationPtr CorrelationPtr;
		static bool HasCorrelationsForThermo;
		
};


class Automorphs
{
	protected:
		const Molecule* mol;
		int NumAuto;
		int Symmetry;
		int RootClass;
		int RootAtom;
		
		map<int,vector<int> > ClassRankAtomMap;
		map<int,set<int> > AtomPredecessorMap;
		vector<int> DistanceFromRoot;
		vector< set<int> >  AllPermutationClasses;
		map<int, vector<set<int> > > NeighboringClassSetsForDepthOne;// for each atom 'i' at depth 1 from root, this map stores the sets of symmetric classes of atoms that are in the atom i's immediate neighbors
		void MapClassRanksAndAtoms();
		void SetRootClassAndAtom();
		void BFSFromRoot();
		void FindAllPermutationClasses();
		int GetSymmetryNumberForAllLeaves();
		int SymmetryCorrections();
		int GetSymmetryNumberFromEqSet(set<int>, int);
		int GetExtDbBondCorrection();
		void traverseDoubleBonds(int, vector<int>*);
		void CalculateAutomorphsAndSymmetryNumbers();
		

	public:
		Automorphs(const Molecule&);
		int NumberOfAutomorphs();
		int SymmetryNumber();
};

class KineticsEstimates 
{
	public:
		KineticsEstimates();
		int kParam;
		double kRelative;
		int EaParam;
		double EaRelative;
		double kRefTemp;
		bool isReverse;
		double TempExp;
};

typedef KineticsEstimates (*AltKineticParamPtr)(vector<Molecule>&, vector<Molecule>&);



class KineticsEvaluationFromParam//this calculates kinetics based on a set of parameters
{
	public:
		KineticsEvaluationFromParam(KineticsEstimates&, vector<double>*);
		KineticsEstimates estimates;
		vector<double>* Param;//the list of parameters is set as a pointer because then IDA can change the parameter value to get sensitivity stuff
		//defining a bunch of base parameter references to get reference kinetics and a bunch of relative values to get the actual rate
		double getKineticsValue(double);//provide temperature 
};


class KineticsInfo
{
	protected:
		vector<KineticParamPtr>* KineticFunctions;
		vector<AltKineticParamPtr>* AltKineticFunctions;
		generated_rxn* rxn;
		double Rg;
		double Temp;
		double RxnEnthalpy;// NOTE - in J/mol
		double RxnEntropy;
		double del_n;
		double siteDensity;
		double calculateBEofMolecules(vector<Molecule>&);
		int firstGasPhaseReactant;
		int firstGasPhaseProduct;
		vector<double>* parameters;
		bool calcCollisionFreq;
		
	public:
		KineticsInfo(vector<KineticParamPtr>&, generated_rxn&, double, double, double, double, double, double,int,int);
		//The above constructor takes in all the kinetics functions, the particular reaction, gas constant, temperature, dHrxn in J/mol, dS in J/mol/K, and delN = the change in moles of gas phase reactants, site density in molecules/cm^2, and an integer that indicates the first gas phase species in reactants, and first in products (only used when sticking frequencies are used) 
		void setCollisionFreq(bool);
		
		KineticsInfo(vector<AltKineticParamPtr>&, vector<double>*, generated_rxn&, double, double, double, double, double);
		/*The above constructor takes in all the kinetics functions in the alternative form, the list of parameters, 
		the particular reaction, gas constant, temperature, dHrxn in J/mol, dS in J/mol/K, and delN = the change in moles of gas phase reactants*/
		bool getKineticParameters(double&, double&, double&, double&,bool&,bool&, bool&, double&, double&, bool&);
		/*the order is Preexponential factor, activation barrier, temperature index, kinetic value, boolean that outputs calcK, boolean that outputs usesBEP, boolean that outputs usesLFER, alpha, beta, boolean that outputs stick*/
		bool getKineticParameters(double&, int&, double&, int&, double&,double&, double&, bool&);
		/*this gives actual kinetics value, index of the parameter list for reference k, relative multiplier for reference k, index of Ea from paramter list, relative adduct to this base Ea in kJ/mol, reference temperature for k, and boolean flag for calculating K for reverse step*/
		/*IMPORTANT NOTE: the alternative kinetics uses Ea in kJ/mol while the original form of kinetics has Ea in J/mol!!*/
};





class UserData
{
	public:	
		UserData(vector<generated_rxn>*, vector<double>*, map<int,double>*, double*, 
			double*, double*, int*, map<int,string*>*, multimap<int,pair<int,int> >*, 
			vector<set<int> >*, double, double, set<int>*, set<int>*, map<int, pair<double,int> >*,
			multimap<string, pair<int,int> >*, multimap<int, string>*, multimap<int,int>*,int*,int*,
			double*,map<int,int>*, vector<KineticsEvaluationFromParam>*, vector<double>);
		vector<generated_rxn>* Reactions;
		double* Temp;
		double* Volume;
		double* Pressure;
		int* NumOfEquations;
		map<int,string*>* SpeciesIndex;
		multimap<int,pair<int,int> >* StoichInfo;
		vector<set<int> >* RateInfo;
		vector<double>* Parameters;
		vector<KineticsEvaluationFromParam>* kinetics;
		double AbsTolerance;
		map<int,double>* RevKinSet;
		double Rg;
		set<int>* surfaceSpecies;
		set<int> * AllSites;
		map<int, pair<double,int> >* InitialMolecFlows;
		multimap<string, pair<int,int> >* SiteOccupants;
		multimap<int, string>* SiteBalanceInfo;
		multimap<int,int>* ReverseRxns;
		int* FirstFlowSpecies;
		int* FirstFoundProductFlowSpecies;
		double* InletMassFlowRate;
		map<int,int>* SpeciesMolWeights;
		vector<double> PreconditionerInverse;
		
		


};



class KineticModel
{
	protected:
		rxn_net_gen* Network;
		char* filename;
		vector<double>* PreExpFactor;
		vector<double>* ActEnergy;
		vector<double>* N;
		vector<KineticParamPtr>* kinetics; // kinetic functions
		vector<AltKineticParamPtr>* AlternativeKinetics; // alternative kinetics functions -- these are set as relative to the base kinetics
		bool useAlternativeKinetics; //boolean flag for alternative kinetics
		vector<double> AlternativeKinParams; //the parameters for alternative kinetics
		map<int,pair<double,int> > InitialMolecFlows;//int - molecule, double - molar flow rate (if int is zero), or site density (if int is one)
		map<string, pair<double,int> > InitialFlowSpecified;//the initial flow rates specified by the user
		multimap<string, string> SiteBalanceInput; //site balance info specified by the user - key is the site while values are strings that signify that molecules containing that string must be included in the site balance
		double Temp;
		double Volume; 
		double InletPressure;
		double Rg;
		int NumOfEquations;
		int NumOfParameters;
		bool ProceedToSolve;
		bool doSensitivity;
		map<int,string*> SpeciesIndex;//maps species integer index and its string pointer
		multimap<int, pair<int,int> > StoichInfo; //key value is the row, mapped value is a pair of the columns (rxns) and actual stoichiometric coefficient!
		//Note that the reactions of a particular molecule are not necessarily arranged in any order.
		vector<set<int> > RateInfo;//each element of the vector corresponds to rate expression of a reaction - the set<int> contains the indices of the species in the network
		vector<double> KineticsValueInfo;//each element of the vector corresponds to kinetic value! -> RateInfo[i]*kineticsInfo[i] is the actual rate of that reaction!
		vector<KineticsEvaluationFromParam> ActualKineticsUsed; //the actual kinetics values used to calculate the residual. Relatved to KineticsValueInfo 
		void GenerateStoichAndRateVectorInfo(map<string,pair<double,int> >&, multimap<string,string>&);//populates the stoich and rate vector datastructures! 
		void SetRateInfo(int);//updates the RateInfo vector for each molecule -updates the sets of those reactions the molecule is a reactant of
		void SetSpeciesIndex();//populates the map SpeciesIndex
		bool SetKinetics();//updates kineticsInfo vector with calculated rate values! 
		int getStoichCoeff(int, int);//get the stoichiometric coeff of the first argument in the reaction corresponding to the second argument
		void SetStoichInfoFromMap(multimap<string*,int>*, int);//used for setting up the StoichInfo for a particular molecule (int) from MolReactantMap and MolProductMap of the Network
		void SetInitialMolecFlows();//sets the InitialMolecFlows		
		static int check_flag(void *, char *, int );
		int getNumEquations();
		ConstrPtr RequiredOutputs; //a constraint pointer to specify outputs
		bool OutputConstraintsSpecified; // a boolean that checks if the constraint pointer has been specified
		set<string> specifiedOutputs; // a set of string that lists user specified outputs.
		vector<int> IndicesOfOutputs;
		bool isRequiredOutput(string);
		void printOutputs(double, N_Vector, N_Vector);
		void PrintSensOutput(double, N_Vector*);
		void PrintDORC(double,N_Vector*, N_Vector);
		void PrintFinalStats(void *);
		map<int,double> RxnsWithRevKin;//keeps track of those reactions that have kinetics defined wrt a reverse, and the K value
		map<double, vector<double> > ReqOutValues;
		set<int> surfaceSpecies; //all surface species
		set<int> AllSites;//keeps track of the index of all the main sites
		multimap<string, pair<int,int> > SiteOccupants; //for each site specification A, gives which species occupies one of the sites of A and how many sites it occupies
		bool PopulateSiteOccupancyMap(string, int, string); //given the species (and it's internal index) and site, update the map with the number of site atoms the species has (i.e. the number of sites the species occupies) - return true if it populates.
		multimap<int, string> SiteBalanceInfo;// Key - input surface site (given in conc units), value - heterogeneous site composite atoms. This provides the info of which site occupants have to be accounted for in the site balance of the key species.
		void findReverseReactions(int);
		multimap<int,int> ReverseRxns;
		string GetMFofSpecies(int);
		bool UseLumpedNetwork;
		int FirstDcldFlowSpecies; //this keeps track of the index of the first initial reactant declared as non-surface flow species -- the flow rate of this species is calculated from mass balance
		int FirstFoundProductFlowSpecies; // this keeps track of the index of the first product that is gas phase
		double InletMassFlowRate;
		map<int,int> MolWtOfSpecies;//key is species index, value is the weight

		map<string*,int,classcomp>* MoleculeSetPtr;
		multimap<string*,int>* MolReactantMapPtr;
		multimap<string*,int>* MolProductMapPtr;
		vector<generated_rxn>* ReactionsPtr;
		vector<double> InitialConditionsFromDyCSTR;
		vector<double> InitialConditionsDerivativeFromDyCSTR;
		double tfinal;
		void PrepareKineticModel();
		map<string, double> OutputValues;
		map<string, vector<double> > SensitivityOutputs; 
		void storeOutputs(N_Vector, N_Vector, N_Vector*);
		bool shouldCalculateRates;
		vector<vector<double> > Rates;
		

		int Solve(int, bool);//argument provides the type -- 0 for normal, 1 for DynamicCSTR. Bool says if it is for setting initial conditions or not
	
	public:
		//KineticModel(rxn_net_gen*, char*, vector<KineticParamPtr>&, double, double, map<string,pair<double,int> >&, multimap<string, string>, ConstrPtr, double, bool);
		KineticModel();
		int SolveModel();
		void calculateCpAt(int);
		void setNetwork(rxn_net_gen*);
		void printTofile(char*);
		void setKineticFns(vector<KineticParamPtr>&);
		void setVolume(double);
		void setPressure(double);
		void setTemperature(double);
		void setInitialFlow(map<string,pair<double,int> >&);
		void setSiteBalance(multimap<string, string>);
		void setOutputConstraint(ConstrPtr);
		void setOutputSpecies(set<string>);
		void setRg(double);
		void setForLumpedNetwork(bool);
		void doSensitivityAnalysis(); //to be specified before solveModel is called
		void setAlternativeKinetics(vector<AltKineticParamPtr>&, vector<double>);
		map<string, double> getModelOutputs();//gives the output flow rates of species
		map<string, vector<double> > getModelSensitivityOutputs();//gives the sensitivity of each species wrt parameters --only for alternative kinetics case; for the other case this will be empty!
		vector<vector<double> >  getRates();
		void calculateRates(); //set this before solve model

		
};

class CHEMKinFiles
{
	protected:
		rxn_net_gen* Network;
		set<int> TempsForCpCalc;
		vector<KineticParamPtr>* kinetics;
		double siteDensity;
		double Rg;
		map<string,double> AllSites;
		multimap<string,string> DenticityInfoMap;
		set<string> SitesRequiringDenticity;

		
		int CreateCHEMKINfileAndWriteSpecies();
		vector<string> CHEMKINRxnStrings;
		vector<string> CHEMKINPreExpInfo;
		vector<string> CHEMKINIndexInfo;
		vector<string> CHEMKINActEInfo;
		set<int> RxnsWithStickingCoeff;

		string generateCHEMKINRxnStrings(int);
		void generateCHEMKINkineticsStrings(double, double, double, bool);
		int WriteRxnInfoInCHEMKINfile();
		map<string*,string, classcomp> CHEMKINSpeciesName;
		map<string*,string, classcomp> CHEMKINBulkName;
		set<string> AllElements;
		void FindAllElements();
		bool isGasPhaseRxn(int);
		bool SetKinetics();
		multimap<string, pair<string*,int> > SiteOccupancyMap;//this just checks for composite heterogeneous atoms
		multimap<string, pair<string*,int> > HeterogeneousSitesOccupancyMap;//this checks for the actually defined heterogeneous sites -note the small but important difference - a heterogeneous site will require a heterogeneous composite site atom but can include more. for e.g. Zeo could be a heterogeneous composite site atom, but the actual free heterogeneous site is [{Zeo}H]
		bool PopulateSiteOccupancyMap(string*, string);
		void PopulateHeterogeneousSiteOccupanyMap(string*);
		set<string*> surfaceSpecies;
		int calculateNetSiteCount(string,int);

	public:
		CHEMKinFiles(rxn_net_gen*, set<int>&,vector<KineticParamPtr>&, double, map<string,double>&,multimap<string,string>&);
		void generateFiles();
};

class GAMSFiles
{
	protected:
		rxn_net_gen* Network;
		char* filename;
		vector<KineticParamPtr>* kineticsFns;
		set<string> AllSites;
		multimap<string,string> DenticityInfoMap;
		set<string> SitesRequiringDenticity;
		set<string*> surfaceSpecies;
		bool PopulateSiteOccupancyMap(string*, string);
		multimap<string, pair<string*,int> > SiteOccupancyMap;
		map<string, double> SiteDensity; //assuming right now that only one type of site exists -- need to allow for more
	public:
		GAMSFiles(rxn_net_gen*, vector<KineticParamPtr>&, set<string>&, map<string, double>&, multimap<string,string>&, char*);
		bool generateFiles();
};

class ThermoGroupsIdentification
{
	
	public:
		ThermoGroupsIdentification(vector<pair<string, string> >&, vector<pair<string, string> >&, char*,char*);
};




