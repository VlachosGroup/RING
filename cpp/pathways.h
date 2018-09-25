#ifndef PATHWAYS_H
#define PATHWAYS_H

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "kinetics.h"
#include "molecule.h"
#include "generated_rxn.h"
#include "additionalfunc.h"
#include "rng.h"


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
        std::vector<generated_rxn> * rxn_index;
        std::vector<double>* RxnEnthalpy;
        std::vector<double>* RxnEntropy;
        std::vector<double>* RxnFreeEnergy;
        std::vector<double>* RxnActE;

	public:
		//RxnPathway(std::vector<generated_rxn>*, std::vector<double>*, std::vector<double>*, std::vector<double>*, std::vector<double>*);
		RxnPathway(std::vector<generated_rxn>*);
		void setRxnActE(std::vector<double>*);
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
		//std::vector<PCRule> RuleConstr;
		//std::vector<PCRuleMol> RulePlusMolConstr;
		//std::vector<PCRuleConstr> RulePlusFeatureConstr;
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
        std::vector<double> ReactionRates;
        std::multimap<std::string*, std::vector<int> > FindAllPathways(std::string*, int);//returns the pathway reaction indices for input molecule string pointer and max path length
        std::string* CloserMolecule(int, std::string*);
		PathwayConstraints Constr;
		bool HasCycles(std::vector<std::string*>); 
		bool isReverseAlreadyPresent(std::vector<int>, int);
		bool similarRulePresent(std::multimap<int,std::string*>*, int, std::string*);
		bool DistinctNature; //set to true to find pathways that are distinct in terms of the number of different elemetary steps
        std::vector<KineticParamPtr>* KineticFunctions;
		bool calculateActEnergy(std::vector<int>&, std::vector<double>&);//setting kinetics for each reaction
		bool canCalculateActE;
		void DominantReaction(std::vector<int>&, std::vector<double>&);
		void DominantReactionClass(std::vector<int>&, std::vector<double>&);
        std::vector<int> MajorDominantRxns(std::vector<int>&, std::vector<double>&);
        std::multimap<std::string*, std::vector<int> > FindAllDominantPathways(std::string*, int);
		void PopulateDominantProductMap(std::string*);
        std::multimap<std::string*, int> MolDominantProductMap; //says which reactions have the dominant rates of formation for a product
	public:
		Pathways(rxn_net_gen*,PathwayConstraints&, const char*);
		void SetConstraints(PathwayConstraints);
		void GeneratePathways();//pathways to products that are final
		//void GeneratePathways(std::vector<std::string>);//pathways to products specified in the argument
		void QueryPlusGenerate(ConstrPtr);//pathways to products that satisfy constraints
		int CalculateCost(std::vector<int>);
		void setDistinctNature(bool);
		void AddKineticFunctions(std::vector<KineticParamPtr>&); // add kinetics info
        std::pair<int,std::map<int,int> > CalculatePathwayValue(std::vector<int>);
		void SetRateVector(std::vector<double>);
		void getDominantPrevStep(std::string);

		
};

class MoleculeQuery
{
	protected:
		const char* filename;
		rxn_net_gen* Network;
        std::vector<std::string*> PassedMolecules;
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



#endif
