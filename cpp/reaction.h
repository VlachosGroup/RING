#ifndef REACTION_H
#define REACTION_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <utility>
#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;

#include "atomcontainer.h"
#include "molecule.h"
#include "patternmatch.h"
#include "generated_rxn.h"

class rxn_net_get;

class Reactiontype
{
	protected:
		std::vector <Substructure> reactant_pattern;//reactant pattern list
		bool IntraMolecularRxnAlso;
		bool IntraMolecularRxnOnly;
        std::vector <bond_operations> bondchanges;// a vector of bond_operations to store connectivity changes
        std::map<int,std::string> mod_atomtype;//a map indicating how an atomtype of an atom with the given label changes
        std::vector <Substructure> struct_constraints;
        std::vector<int> FragmentCopyIndex;//-1 for no duplicate, 0 for the first reactant, 1 for the second reactant
        std::vector< std::map<int,int> > FragmentCopyLabels;
        std::list<ConstrPtr> RxnConstraints;
		ConstrPtr ProductConstraints;
		CombinedConstrPtr CombinedConstraint;
		void add_bondchanges(bond_operations);
		int Cost;
        std::string RuleName;//name assigned to the rule
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
		void add_mod_atomtype(int,std::string);
		int get_molecularity();
		void AllowIntraMolecularRxnOnly();
		void AllowIntraMolecularRxnAlso();
		bool isIntraMolecularAlso();
		bool isIntraMolecularOnly();
		void AddReactantCopy(int, std::map<int,int>);//specifies which reactant is copied and what the atom maps are
		void disconnect_bond(int,int);
		void connect_bond(int,int);
		void increaseBO(int,int, int);
		void decreaseBO(int,int, int);
		void setCost(int);
		int getCost();
		void setRuleName(std::string);
        std::string getRuleName();
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


class Reaction
{
	protected:
        std::vector<Molecule> M;
		Reactiontype Rt;
        std::vector<generated_rxn> gen_rxns;
        std::vector <Patternmatch> reactpattlist;
		bool GenerateIntraMolecularRxnsAlso;
		bool GenerateIntraMolecularRxnsOnly;
		bool OnlyOneReactantUsed;
		void generate_reactions();
		void findmatch(int);		
		void addproducts(Atomcontainer&, std::vector<Molecule>&, std::vector<std::pair<std::string, int> >&, int);
		bool perform_changes(Atomcontainer&, std::vector<std::pair<std::string, int> >&);//performs the atomtype modifications and bond order changes. Returns true if changes can be performed, false if it runs into a problem
		void GenerateRxnWithFragmentCopies(Atomcontainer, int, int);

	public:
		Reaction (std::vector<Molecule> &, Reactiontype &, std::vector<Patternmatch>& );
		generated_rxn get_generated_rxns(int);
		int number_rxns_generated();
}; 

#include "rng.h"
#endif
