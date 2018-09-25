#ifndef KINETICS_H
#define KINETICS_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>

#include "molecule.h"
#include "generated_rxn.h"
#include "rng.h"

#include <nvector/nvector_serial.h>


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


typedef KineticsEstimates (*AltKineticParamPtr)(std::vector<Molecule>&, std::vector<Molecule>&);



typedef bool (*KineticParamPtr)(std::vector<Molecule>&, std::vector<Molecule>&, double&, double&, double&, bool&, bool&, bool&, double&,double&,bool&,double);

class KineticsInfo
{
	protected:
		std::vector<KineticParamPtr>* KineticFunctions;
		std::vector<AltKineticParamPtr>* AltKineticFunctions;
		generated_rxn* rxn;
		double Rg;
		double Temp;
		double RxnEnthalpy;// NOTE - in J/mol
		double RxnEntropy;
		double del_n;
		double siteDensity;
		double calculateBEofMolecules(std::vector<Molecule>&);
		int firstGasPhaseReactant;
		int firstGasPhaseProduct;
		std::vector<double>* parameters;
		bool calcCollisionFreq;
		
	public:
		KineticsInfo(std::vector<KineticParamPtr>&, generated_rxn&, double, double, double, double, double, double,int,int);
		//The above constructor takes in all the kinetics functions, the particular reaction, gas constant, temperature, dHrxn in J/mol, dS in J/mol/K, and delN = the change in moles of gas phase reactants, site density in molecules/cm^2, and an integer that indicates the first gas phase species in reactants, and first in products (only used when sticking frequencies are used) 
		void setCollisionFreq(bool);
		
		KineticsInfo(std::vector<AltKineticParamPtr>&, std::vector<double>*, generated_rxn&, double, double, double, double, double);
		/*The above constructor takes in all the kinetics functions in the alternative form, the list of parameters, 
		the particular reaction, gas constant, temperature, dHrxn in J/mol, dS in J/mol/K, and delN = the change in moles of gas phase reactants*/
		bool getKineticParameters(double&, double&, double&, double&,bool&,bool&, bool&, double&, double&, bool&);
		/*the order is Preexponential factor, activation barrier, temperature index, kinetic value, boolean that outputs calcK, boolean that outputs usesBEP, boolean that outputs usesLFER, alpha, beta, boolean that outputs stick*/
		bool getKineticParameters(double&, int&, double&, int&, double&,double&, double&, bool&);
		/*this gives actual kinetics value, index of the parameter list for reference k, relative multiplier for reference k, index of Ea from paramter list, relative adduct to this base Ea in kJ/mol, reference temperature for k, and boolean flag for calculating K for reverse step*/
		/*IMPORTANT NOTE: the alternative kinetics uses Ea in kJ/mol while the original form of kinetics has Ea in J/mol!!*/
};


class KineticsEvaluationFromParam//this calculates kinetics based on a set of parameters
{
	public:
		KineticsEvaluationFromParam(KineticsEstimates&, std::vector<double>*);
		KineticsEstimates estimates;
        std::vector<double>* Param;//the list of parameters is set as a pointer because then IDA can change the parameter value to get sensitivity stuff
		//defining a bunch of base parameter references to get reference kinetics and a bunch of relative values to get the actual rate
		double getKineticsValue(double);//provide temperature 
};


class UserData
{
	public:	
		UserData(std::vector<generated_rxn>*, std::vector<double>*, std::map<int,double>*, double*, 
			double*, double*, int*, std::map<int,std::string*>*, std::multimap<int,std::pair<int,int> >*, 
			std::vector<std::set<int> >*, double, double, std::set<int>*, std::set<int>*, std::map<int, std::pair<double,int> >*,
			std::multimap<std::string, std::pair<int,int> >*, std::multimap<int, std::string>*, std::multimap<int,int>*,int*,int*,
			double*,std::map<int,int>*, std::vector<KineticsEvaluationFromParam>*, std::vector<double>);
        std::vector<generated_rxn>* Reactions;
		double* Temp;
		double* Volume;
		double* Pressure;
		int* NumOfEquations;
        std::map<int,std::string*>* SpeciesIndex;
        std::multimap<int,std::pair<int,int> >* StoichInfo;
        std::vector<std::set<int> >* RateInfo;
        std::vector<double>* Parameters;
        std::vector<KineticsEvaluationFromParam>* kinetics;
		double AbsTolerance;
        std::map<int,double>* RevKinSet;
		double Rg;
        std::set<int>* surfaceSpecies;
        std::set<int> * AllSites;
        std::map<int, std::pair<double,int> >* InitialMolecFlows;
        std::multimap<std::string, std::pair<int,int> >* SiteOccupants;
        std::multimap<int, std::string>* SiteBalanceInfo;
        std::multimap<int,int>* ReverseRxns;
		int* FirstFlowSpecies;
		int* FirstFoundProductFlowSpecies;
		double* InletMassFlowRate;
        std::map<int,int>* SpeciesMolWeights;
        std::vector<double> PreconditionerInverse;
};




class KineticModel
{
	protected:
		rxn_net_gen* Network;
		char* filename;
        std::vector<double>* PreExpFactor;
        std::vector<double>* ActEnergy;
        std::vector<double>* N;
        std::vector<KineticParamPtr>* kinetics; // kinetic functions
        std::vector<AltKineticParamPtr>* AlternativeKinetics; // alternative kinetics functions -- these are set as relative to the base kinetics
		bool useAlternativeKinetics; //boolean flag for alternative kinetics
        std::vector<double> AlternativeKinParams; //the parameters for alternative kinetics
        std::map<int,std::pair<double,int> > InitialMolecFlows;//int - molecule, double - molar flow rate (if int is zero), or site density (if int is one)
        std::map<std::string, std::pair<double,int> > InitialFlowSpecified;//the initial flow rates specified by the user
        std::multimap<std::string, std::string> SiteBalanceInput; //site balance info specified by the user - key is the site while values are strings that signify that molecules containing that string must be included in the site balance
		double Temp;
		double Volume; 
		double InletPressure;
		double Rg;
		int NumOfEquations;
		int NumOfParameters;
		bool ProceedToSolve;
		bool doSensitivity;
        std::map<int,std::string*> SpeciesIndex;//maps species integer index and its string pointer
        std::multimap<int, std::pair<int,int> > StoichInfo; //key value is the row, mapped value is a pair of the columns (rxns) and actual stoichiometric coefficient!
		//Note that the reactions of a particular molecule are not necessarily arranged in any order.
        std::vector<std::set<int> > RateInfo;//each element of the vector corresponds to rate expression of a reaction - the set<int> contains the indices of the species in the network
        std::vector<double> KineticsValueInfo;//each element of the vector corresponds to kinetic value! -> RateInfo[i]*kineticsInfo[i] is the actual rate of that reaction!
        std::vector<KineticsEvaluationFromParam> ActualKineticsUsed; //the actual kinetics values used to calculate the residual. Relatved to KineticsValueInfo 
		void GenerateStoichAndRateVectorInfo(std::map<std::string,std::pair<double,int> >&, std::multimap<std::string,std::string>&);//populates the stoich and rate vector datastructures! 
		void SetRateInfo(int);//updates the RateInfo vector for each molecule -updates the sets of those reactions the molecule is a reactant of
		void SetSpeciesIndex();//populates the map SpeciesIndex
		bool SetKinetics();//updates kineticsInfo vector with calculated rate values! 
		int getStoichCoeff(int, int);//get the stoichiometric coeff of the first argument in the reaction corresponding to the second argument
		void SetStoichInfoFromMap(std::multimap<std::string*,int>*, int);//used for setting up the StoichInfo for a particular molecule (int) from MolReactantMap and MolProductMap of the Network
		void SetInitialMolecFlows();//sets the InitialMolecFlows		
		static int check_flag(void *, char *, int );
		int getNumEquations();
		ConstrPtr RequiredOutputs; //a constraint pointer to specify outputs
		bool OutputConstraintsSpecified; // a boolean that checks if the constraint pointer has been specified
        std::set<std::string> specifiedOutputs; // a set of string that lists user specified outputs.
        std::vector<int> IndicesOfOutputs;
		bool isRequiredOutput(std::string);
		void printOutputs(double, N_Vector, N_Vector);
		void PrintSensOutput(double, N_Vector*);
		void PrintDORC(double,N_Vector*, N_Vector);
		void PrintFinalStats(void *);
        std::map<int,double> RxnsWithRevKin;//keeps track of those reactions that have kinetics defined wrt a reverse, and the K value
        std::map<double, std::vector<double> > ReqOutValues;
        std::set<int> surfaceSpecies; //all surface species
        std::set<int> AllSites;//keeps track of the index of all the main sites
        std::multimap<std::string, std::pair<int,int> > SiteOccupants; //for each site specification A, gives which species occupies one of the sites of A and how many sites it occupies
		bool PopulateSiteOccupancyMap(std::string, int, std::string); //given the species (and it's internal index) and site, update the map with the number of site atoms the species has (i.e. the number of sites the species occupies) - return true if it populates.
        std::multimap<int, std::string> SiteBalanceInfo;// Key - input surface site (given in conc units), value - heterogeneous site composite atoms. This provides the info of which site occupants have to be accounted for in the site balance of the key species.
		void findReverseReactions(int);
        std::multimap<int,int> ReverseRxns;
        std::string GetMFofSpecies(int);
		bool UseLumpedNetwork;
		int FirstDcldFlowSpecies; //this keeps track of the index of the first initial reactant declared as non-surface flow species -- the flow rate of this species is calculated from mass balance
		int FirstFoundProductFlowSpecies; // this keeps track of the index of the first product that is gas phase
		double InletMassFlowRate;
        std::map<int,int> MolWtOfSpecies;//key is species index, value is the weight

        std::map<std::string*,int,classcomp>* MoleculeSetPtr;
        std::multimap<std::string*,int>* MolReactantMapPtr;
        std::multimap<std::string*,int>* MolProductMapPtr;
        std::vector<generated_rxn>* ReactionsPtr;
        std::vector<double> InitialConditionsFromDyCSTR;
        std::vector<double> InitialConditionsDerivativeFromDyCSTR;
		double tfinal;
		void PrepareKineticModel();
        std::map<std::string, double> OutputValues;
        std::map<std::string, std::vector<double> > SensitivityOutputs; 
		void storeOutputs(N_Vector, N_Vector, N_Vector*);
		bool shouldCalculateRates;
        std::vector<std::vector<double> > Rates;
		

		int Solve(int, bool);//argument provides the type -- 0 for normal, 1 for DynamicCSTR. Bool says if it is for setting initial conditions or not
	
	public:
		//KineticModel(rxn_net_gen*, char*, vector<KineticParamPtr>&, double, double, std::map<std::string,std::pair<double,int> >&, std::multimap<std::string, std::string>, ConstrPtr, double, bool);
		KineticModel();
		int SolveModel();
		void calculateCpAt(int);
		void setNetwork(rxn_net_gen*);
		void printTofile(char*);
		void setKineticFns(std::vector<KineticParamPtr>&);
		void setVolume(double);
		void setPressure(double);
		void setTemperature(double);
		void setInitialFlow(std::map<std::string,std::pair<double,int> >&);
		void setSiteBalance(std::multimap<std::string, std::string>);
		void setOutputConstraint(ConstrPtr);
		void setOutputSpecies(std::set<std::string>);
		void setRg(double);
		void setForLumpedNetwork(bool);
		void doSensitivityAnalysis(); //to be specified before solveModel is called
		void setAlternativeKinetics(std::vector<AltKineticParamPtr>&, std::vector<double>);
        std::map<std::string, double> getModelOutputs();//gives the output flow rates of species
        std::map<std::string, std::vector<double> > getModelSensitivityOutputs();//gives the sensitivity of each species wrt parameters --only for alternative kinetics case; for the other case this will be empty!
        std::vector<std::vector<double> >  getRates();
		void calculateRates(); //set this before solve model

		
};

#endif
