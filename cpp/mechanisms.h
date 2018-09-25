#ifndef MECHANISM_H
#define MECHANISM_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <vector>
#include <utility>
//#include <list>
#include <set>
#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;
//
//#include <stdio.h>
//#include <stdlib.h>
//
//#include <idas/idas.h>
//#include <idas/idas_dense.h>
//#include <nstd::vector/nstd::vector_serial.h>
//#include <sundials/sundials_math.h>
//#include <sundials/sundials_types.h>
//
//#include "common.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "atomcontainer.h"
//#include "automorphs.h"
//#include "substructure.h"
//#include "patternmatch.h"
#include "molecule.h"
#include "pathways.h"
//#include "reaction.h"
//#include "lumping.h"
//#include "groupadditivity.h"
#include "generated_rxn.h"
#include "rng.h"

typedef bool (*RxnConstraintPtr)(RxnInfo&);

class partialMechanism
{
	protected:
		std::vector < std::pair<generated_rxn,int> >* rxn_index;	
        std::map<std::string*,int> stoichiometry;
		std::vector<double>* RxnActE;

	public:
		partialMechanism(std::vector < std::pair<generated_rxn,int> >*);
		bool isContained(partialMechanism&);
		bool isNetZero();
        std::string printformula();
        std::string printSMILES();
        std::map<std::string*, int> getStoichiometry();
		int getFreq(std::string*);
		bool isEqual(partialMechanism*);//is equal in terms of overall stoichiometry!
		int numRule(int);  // index of rule to get number of occurences
		int numRuleWithParticipants(int, ConstrPtr, int);  // index of rule, molecule constraint, and an integer indicator to indicate constraints is on reactants (0) or products (1) to get number of occurences
		int numRuleWithParticipants(int, CombinedConstrPtr, int);  // index of rule, combined molecule constraint, and an integer indicator as above to get number of occurences
		int numMolecule(ConstrPtr); // how often this molecule appears in this mechanism.
		int numMoleculeOverall(ConstrPtr);//similar to numMolecule, applied to the overall rxn, the second argument says reactant (0) or product (1)
        std::set<std::string> getProducts(); //gets all the products of this partial mechanism 
		void MechinfoForAthena();
		int numSatisfying(int, RxnConstraintPtr); //number of instances of reactions of a given rule satifying some reaction constraint parameters
		void setRxnActE(std::vector<double>*);

		friend class Mechanisms;
};




class OverallMechanism
{
	protected:
		std::vector<std::pair<std::string*, int> > MoleculeDirectMechPair;/*stores a std::pair of which molecule's direct mechanism is part of the overall mechanism and the integer index corresponding to the direct mechanism stored in 
															MoleculeMechanismsMap of class Mechanisms*/
		std::vector<int> Frequency;//stores how many of each direct mechanism is to be counted in the overall mechanism	
	public:
		OverallMechanism();
		void insertDirectMech(std::string*,int,int);
		std::pair<std::string*,int> getDirectMechanism(int);
		int NumberOfDirectMechs();
		int getFrequency(int);
        std::map<std::string*, std::pair<int,int> > generateStringPairMap();
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
        std::map<std::string*,int> SpeciesStoichio;
	public:
		OverallRxn(partialMechanism*);
		int getFreq(std::string*);
};

class Mechanisms
{
	protected:
		rxn_net_gen* Network;
		const char* filename;
        std::map<std::string*, std::vector< std::map< int, int> > > MoleculeMechanismsMap;
		std::vector< std::map<int,int> > FindDirectMechanisms(std::string*, std::string*,bool);
		std::vector<OverallMechanism> AllMechanisms;
        std::set< std::map<std::string*, std::pair<int,int> > > UniqueMechs;
		bool DistinctDirectMechs;
		DirectMechConstraints DirectConstr;
		CompleteMechConstraints CompleteConstr;
        std::map<int,int> getMechanismMap(std::vector <std::pair<int, int> >&);
		bool willFormInternalCycles(generated_rxn, std::vector< std::pair<generated_rxn, int> >&);
		std::vector<KineticParamPtr>* KineticFunctions;
		bool calculateActEnergy(std::vector<std::pair<int,int> >&, std::vector<double>&);//setting kinetics for each reaction
		bool canCalculateActE;
        std::map<int,double> RxnActEnergyMap;//stores the activation barrier of each reaction -note that this reaction index starts from 1!!!. 

	public:
		Mechanisms (rxn_net_gen*, const char*);
		void GenerateMechanisms(ConstrPtr );
		void GenerateDirectMechanismsOlny(ConstrPtr);
		bool isComplete(partialMechanism*);
		void SetDistinctDirectMech(bool);
		bool ContainsInitialReactantsAsProducts(partialMechanism*);
		bool FoundOverallMechanism(partialMechanism*);
		int GenerateOverallMechanisms(std::string*);
		int getTotalCost(std::vector< std::pair<generated_rxn, int> > *);
        std::map<std::string*,int> getIntermediates(partialMechanism*);
		bool isReverseAlreadyPresent(std::vector<std::pair<generated_rxn,int> > *, generated_rxn *);
		void SetDirectMechConstraints(DirectMechConstraints);
		void SetCompleteMechConstraints(CompleteMechConstraints);
		void AddKineticFunctions(std::vector<KineticParamPtr>&); // add kinetics info
};


#endif
