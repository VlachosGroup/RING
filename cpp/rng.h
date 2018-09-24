#ifndef RXN_NET_GEN_H
#define RXN_NET_GEN_H

#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <utility>

//#include <Python.h>
//#include "..\python\pyreactiontype.h"

#include "molecule.h"
#include "reaction.h"
#include "lumping.h"
#include "groupadditivity.h"
#include "generated_rxn.h"

class rxn_net_gen
{
	protected:
		std::multimap<int, Molecule*> unprocessedmol; //molecules yet to be processed
		std::vector<Molecule*> processedmol; //molecules already processed
        std::set<std::string*, classcomp> InitialReactants; //the initial reactants input by the user
        std::list<std::string> inputStrings; //list of strings corresponding to the input specified by the user (note these are not in canonical SMILES necessarily)
		std::vector<std::pair<std::string, SiteType> > CompositeSites; //all composite atoms that have been specified as sites
		std::vector<std::string> CompositeAtoms; // all composite atoms 
		std::multimap<int, int> ReactionsMap; //a map of rules to all reactions (referred to by their reaction number NOT index, so starts from one) belonging to that rule.
		std::vector<generated_rxn> AllReactions; // an array of all generated reactions obtained through network generation
		std::vector<generated_rxn> ReconstructedReactions; // an array of all reactions generated through the process of reconstruction (populated only when reconstructed from the original network - kicked in ONLY when lumping is OFF)
		std::vector<generated_rxn> LumpedReconstructedNetwork; //final array of all lumped and reconstructed network.
		std::vector<LumpInfo> MolLumps; //An array of different lumps
        std::map<LumpedReaction, int, LumpedReactionCompare> LumpedReactionMap; //A map of each lumped reaction and it's frequency: consider moving it in to LumpedReaction() function
		std::multimap<int,int>LumpedReconstructedReactionsMap; //similar to reactions map but for LumpedReconstructed reactions
		
        std::map<PONALumpingInfo,int> ParaffinSizeLumpMap;//map with key as pair of Size-Hydrogen count pair and PONAElectronicValue; value = lump;
        std::map<PONALumpingInfo,int> OlefinSizeLumpMap;
        std::map<PONALumpingInfo,int> AlkylAromaticsBranchLumpMap;
        std::map<PONALumpingInfo,int> AlkenylAromaticsBranchLumpMap;
        std::map<PONALumpingInfo,int> NapthenicsSizeLumpMap;
        std::map<MoreAdditionalLumpingInfo,int,MoreLumpingCompare> MoreAdditionalLumpingSizeLumpMap; //for all additional specifications to do more lumping --in addition to PONA
        std::map<int,int> AdditionalLumpMap; //a map of additional lumps that maps the MolLumps lumps to other MolLumps lumps
        std::set<int> ToUndergoAdditionalLumping; //stores a list of MolLumps indices that need to undergo additional lumping


        std::map<std::string*,int, classcomp> AllMolecules;
        std::map<std::string*,int, classcomp> AllLumpedMolecules;//this keeps the final map of strings representative of the lump and it's rank - much like AllMolecules, except this is for the final set of lumps! 
		
        std::set<std::string*,classcomp> Intermediates;
		std::multimap<UnsignedIntPair, int > LumpHashMap;// map of hash value and the corresponding lump - given by the index of the lump in MolLumps
        std::map<std::string*,int, classcomp> MolLumpMap; // map of molecule and lump - given by the index in MolLumps

		std::multimap<std::string*,int> MolReactantMap;//says which reactions have the molecule as a reactant
		std::multimap<std::string*,int> MolProductMap;//says which reactions have the molecule as a product

		std::multimap<std::string*,int> MolLumpReactantMap;//says which lumpedreconstructed reaction has molecule as a reactant
		std::multimap<std::string*,int> MolLumpProductMap; // says which lumpedreconstructed reaction has a molecule as a product
		
        std::map<int,int> AtomtypeIndex;
		std::vector<Reactiontype > Rtlist;
        std::map<int, int> AtomtypingMap;//first int stores the value of the atomtype and the second stores its index corresponding to a number in the prime array. 
		Patternmatch checkmatch(std::vector<Molecule>&,int&, int);
		Patternmatch check_reactant_pattern(Molecule&, Substructure&);//checks for presence of reactant pattern instance
		ConstrPtr GlobalConstraints;//global constraints
		LumpingStrategy LumpStrat;
		bool shouldCalcThermo;
        std::map<std::string*, ThermoValues, classcomp> AllMolThermo;
		std::vector<int> RxnsWithInterchangeableReactants;//reactions where the reactants can both swap places but will lead to the same reaction. 
		
		double Temperature;
		bool isunique(std::string);
		bool globalconstraintcheck(Molecule&);
		bool check_combined_match(std::vector<Molecule>&,int&);
		bool check_product_constraints(Molecule&,int);
		//pair<unsigned int, unsigned int> CreateHashValue(const Molecule&);
		int GetAtomValue(std::string);
		void LumpMolecule(const Molecule &, std::string*);
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

        std::string getRxnStringFromLumpedRxn(LumpedReaction);
		int findLumpofMolecule(std::string*);
		bool isSiteIntermediate(std::string);//this checks if a species (that is not an initial reactant) is an intermediate
		bool containsSiteAtom(std::string); // this checks if a species includes a composite site. -- a superset of isIntermediate
		double calculateThermoOfRxn(ThermoType, generated_rxn&, double);
		void add_unique_molecules_reactions(Reaction&, int, bool);//add the new molecule and reactions into the container of molecules and reactions. the int refers to the Reactiontype index. the bool says if we need to keep track of reactions where the popped molecule can participate as either reactant.
		void UpdateMolAsReactantsProductsMap(generated_rxn&);
		void GetHighestRankInfoForRxn(generated_rxn&, int&, int&, bool&);
		void GetProductParentsInfoForRxn(generated_rxn&, std::map<int,int>&, int, bool, int);

		bool setInitialReactants();
		bool MolValencyErrorStatements(std::pair<std::string, int>&);
		void DoPONALumping();
		void DoMoreLumping();
		void printOutputInfo();
		void calculateThermoValues();
		void GenerateMonoMolecularRxns(int, std::vector<Molecule>&, std::vector<Patternmatch>&);
		void GenerateBimolecularRxns(int, std::vector<Molecule>&, std::vector<Patternmatch>&, std::vector<std::vector<std::pair<Molecule*,Patternmatch> > >&, int&);
		void ReconstructForOriginalNetwork(int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, std::string);
		void ReconstructForLumpedNetwork (int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, std::string);
		generated_rxn GetOneRxnFromTwo(generated_rxn&, generated_rxn&);
		void SetAllLumpedMoleculesInfo(generated_rxn& r1);
        std::set<int> RulesRemovedByReconstruction;
		bool calcMolThermo(std::string*,ThermoType,double,double&);
		//std::vector<KineticParamPtr>* kinetics;

		int getDeltaNGasPhase(generated_rxn&);
		int firstGasSpecies(int,int);//first argument is reaction index, second is 0,1 for reactant or product

		/* -------for getting network from files ------*/

		bool ReadReactionFile(const char*, std::map<int,std::string*>&);
		bool ReadSpeciesFile(const char*, std::map<int,std::string*>&);
		bool TraverseNetworkToGetRanks();
		void LumpAllMolecules();

		bool InterpretRxnsAndUpdate(std::vector<std::string>&, int, std::map<int,std::string*>&,int);

		/*-------------------------*/

		
		/*------- GAMS --------*/
		bool DoSimultaneousRxns;
		std::vector< std::set<int> > SimultaneousRxns;
        std::map<int, std::vector<int> > RxnsToBeSummedUpMap;
        std::set<int> SimultRxnsGenerated;
        std::set< std::set<std::string> > reactantpairsForBiMol;//pairs of reactants for the case of bimolecular reactions
		int FindRxnWithReactantsAndRule(const std::set<std::string>&, int);
        std::set<int> FindReactionsWithReactantsAndRule(const std::set<std::string>&, int);
		void findRxnsOfSimultRules(std::vector<std::pair<int,int> >&);
		void findTransformRxns(std::vector<std::pair<int,int> >&);
		void IdentifySimultaneousBimolRxns();
        std::map<std::string*, std::string, classcomp> SMILESGAMSSpeciesMap;
		std::multimap<std::string, std::string> DenticityInfoMap;
		int CalculateNetStoichDifferenceWithDenticity(generated_rxn&, std::string);

		/*-------- end ----------*/

	public:
		rxn_net_gen(std::list<std::string>&, std::vector<Reactiontype >&, ConstrPtr, LumpingStrategy&, std::vector<std::string>&, std::vector<std::string >&, bool);
		//Incomplete -- the above constructor takes in a list of strings (input SMILES strings), std::vector of all the reaction rules, a constraint pointer for global constraints, Lumping strategy and a list of composite atoms-sites
		
		rxn_net_gen(){DoSimultaneousRxns = false; Temperature =298.0;}
		void AddInitialReactants(std::list<std::string>&);//add inital reactants
		void AddReactionRules(std::vector<Reactiontype>&);//add reaction rule std::vector
		void AddGlobalConstraints(ConstrPtr);//add global constraints pointer
		void AddLumpingStrategy(LumpingStrategy&); // add lumping strategy info
		void AddCompositeAtoms(std::vector<std::string>&);// add composite atoms
		void AddCompositeSites(std::vector<std::pair<std::string, SiteType> >&);// add composite atoms that are specifically sites
		void SetCalcThermo(bool);//set to true/ false if thermo has to be calculated
		void GenerateNetwork();//begin generating the network;
		void SetAllLumpedRxns();//sets lumped reactions into a new std::vector of reactions -- NOTE that this SHOULD be called after network generation!
		void GetNetworkFromFile(const char*, const char*);
		rxn_net_gen(rxn_net_gen*, LumpingStrategy&); //TODO
		void generateStoichMatrix();
		void ReconstructReactions(int,int,ConstrPtr, ConstrPtr, CombinedConstrPtr, std::string);
		void print_rxnlist();
		generated_rxn getReaction(int);
		int findSpeciesRank(std::string);
		void GenerateInfoForAthena();
		void SetTemp(double);
		void StoreRxnsAndSpecies(const char*, const char*, bool);
		//void setKinetics(std::vector<KineticParamPtr>&);

        std::pair<unsigned int, unsigned int> CreateHashValue(const Molecule&);
		


		
		//for GAMS
		void CalculateParametersFileForGAMS();
		void GatherSimultaneousAndJointReactionsForGAMS(std::vector<std::pair<int,int> >&, std::vector<std::pair<int,int> >& );
		void setSimultaneousRxns(bool value){DoSimultaneousRxns = value;}
		void PrepareSpeciesInfoForGAMS();
		void CheckDensityInfo(std::multimap<std::string,std::string>&);


		//---------------------------------------------------------

		//friend class Pathways;	
		//friend class Mechanisms;
		//friend class MoleculeQuery;
		//friend class ReactionQuery;
		//friend class KineticModel;
		//friend class KineticsInfo;
		//friend class CHEMKinFiles;
		//friend class GAMSFiles;

};

#endif
