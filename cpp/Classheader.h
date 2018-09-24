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

#include "common.h"
#include "clonable.h"
#include "element.h"
#include "atom.h"
#include "atomcontainer.h"
#include "automorphs.h"
#include "substructure.h"
#include "patternmatch.h"
#include "molecule.h"
#include "reaction.h"
#include "lumping.h"
#include "groupadditivity.h"
#include "generated_rxn.h"
#include "rng.h"

enum AstNodeType
{
	Undefined, OperatorAND, OperatorOR,OperatorNOT,Boolean
};

enum TokenType
{
	Error, AND, OR, NOT, EndOfText, Openbraces, Closedbraces, Size, Charge, IsAromatic, IsCyclic, MinRingSize, MaxRingSize, Fragment
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



typedef bool (*KineticParamPtr)(vector<Molecule>&, vector<Molecule>&, double&, double&, double&, bool&, bool&, bool&, double&,double&,bool&,double);


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




