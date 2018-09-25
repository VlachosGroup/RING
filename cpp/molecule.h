#ifndef MOLECULE_H
#define MOLECULE_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <sstream>
#include <vector>
#include <utility>
//#include <list>
#include <set>
//#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>

#include "atomcontainer.h"

class Molecule: public Atomcontainer//class Molecule derived from Atomcontainer
{
	protected:
		std::string smilesstring;//contains the SMILES like description of pattern
		void dfsvisit();//function for DFS traversing
		void find_all_rings();//function to find the Set of All Rings (SAR)
		void find_aromatic_rings();//find which of the SAR are aromatic rings 
		void find_allylic_atoms();//finds which of the carbon atom is allylic
		std::vector <int> aromatic_rings;//std::vector containing the indices of Allrings that are aromatic
		std::string MolecularFormula;//stores the molecular formula - in the order Carbon, Hydrogen,Oxygen, Nitrogen, Sulphur, and Phosphorous. Following this will be compositeatomtypes. 
		void update_aromaticity_details();//function that makes appropriate updation upon finding/ noting aromatic atoms
        std::set<int> allylic_atoms;//std::vector containing the indices of allylic atoms;
		void EvaluateMF();//evaluates the Molecular Formula of the molecule
		std::string SetWithoutHydrogenCount(std::string);//removes from the std::string, any Hydrogen and number following it within square brackets 
		std::string SetWithoutSquareBrackets(std::string);//removes from the std::string any Square Brackets
		std::string SetWithHydrogenCount(std::string, int);//adds into the std::string the Hydrogen count before "]" if any
		std::string SetWithSquareBrackets(std::string);//adds into the std::string Square Brackets. This assumes ring identifiers are not yet appended to the std::string
		void UpdateSquareBraces();//updates the atomsymbol to remove unwanted square braces and add in necessary ones
		void ModEqRankUsingPrev(std::vector<int>&);//modifies the ranks of atoms if two atoms have equal ranks in one iteration of canonical SMILES algorithm but had unequal ranks previously
		std::vector<int> RankOrder(std::vector<int>&, std::vector<int>&) const;
		void findRanksOfAllAtoms();
		
	public:
		Molecule(std::string, int);//constructor
		Molecule(Atomcontainer&);
		std::string moleculestring() const;
		std::string unique_smiles();//generate unique smiles
		std::string unique_smiles(int);
		void print_smiles();//print smilesstd::string function
		void set_smiles(std::string);// set the smilesstd::string
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
		int totalupElectrons() const;//calculates the total unstd::paired electron count
		int MaxRingsize() const;//Calculates the largest Ring size
		int MaxRingSizeOfAtom(int) const;//calculates the largest ring in which the atom resides
		int SmallestRingSize() const;//calculates the size of the smallest ring		
		int AdjacentBondInRingWithOrder(Path&, int, int ) const;//counts the number of adjacent bonds of a bond type (third arg) present in a ring (first arg) for an atom (second arg)
		int totalDoubleBonds() const;//calculates the total number of double bonds
		int totalTripleBonds() const;//calculates the total number of triple bonds
		int totalAromaticBonds() const; //calculates the total number of aromatic bonds
		bool HasBridgeHead() const;//checks if the molecule has a bridge head
		IntPair Hydrogen_counter(int) const;//return a std::pair- first being total H count, second being D count
		std::string GetMF() const;//get molecular formula
		std::string getMFwithoutAtoms(std::set<std::string>) const;
		std::pair<std::string,int> checkValencyError() const;//throws a std::pair <std::string,int> for error in valency!  
		std::pair<int,int> calcBranchandDistanceValue() const;//calculates a value used for lumping based on distances between leaves etc.
		bool isIntermediate() const;//checks if a molecule is an intermediate --if molecule has a charge, unstd::paired electron or has non-bonded interactions
		bool HasCompositeAtoms() const;//checks if the molecule has a composite atom
		int NumberOfAromaticRings() const;
		int NumberOfRings() const;
		int totalAtomsIncludingH() const;
		int totalAtomsOfType(std::string) const;//counts number of atoms of a given atomtype
		int MolecularWeight() const;
		bool isAtomChiral(int) const;
		int NumberChiralAtoms() const;
		int OptIsomers() const; //gives the number of optical isomers - min 1. 
		int NumberOfLeaves() const;
		
		friend class Patternmatch;//Patternmatch is a friend class
};

typedef bool (*ConstrPtr)(const Molecule &);
typedef bool (*CombinedConstrPtr)(const Molecule &, const Molecule &);

#endif
