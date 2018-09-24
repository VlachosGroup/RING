#ifndef GROUPADDITIVITY_H
#define GROUPADDITIVITY_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
#include <sstream>
#include <vector>
#include <utility>
//#include <list>
#include <set>
#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;

//#include "common.h"
//#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "atomcontainer.h"
#include "molecule.h"
//#include "patternmatch.h"
//#include "substructure.h"

enum ThermoType {EnthalpyType, EntropyType, CpType, FreeEnergyType, logPType};

class GATempDiscretePts
{
	protected:
		static std::vector<double> TempPts;
    public:
		static std::vector<double> getTemperaturePoints();
		static void InsertTempPoints(double);
};

class GAdata
{
	protected:
		double Enthalpy;
		double Entropy;
		double logP;
        std::map<int, double> CpMap;
		std::vector<double> Cp;//assumed to be start from 298
		std::vector<double> A;
		std::vector<double> B;
		//double getCpAtEndPoints(double);
		void CalculateAandBForEachRange();
		double getTempCorrection(ThermoType, double) const;
		double getCpAtArbitraryTemp(int) const;

	public:
		GAdata(std::vector<double>);
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

typedef std::map<std::string, GAdata> GAMap; //the actual std::map that stores the groups and its additivity value
typedef std::map < std::pair<int,int>, GAMap > HashGAMap;// this std::maps HashPair with its appropriate GAMAP --> the HashPair keeps track of the num of double and triple bonds of the neighbors and if the neighboring atoms have special electronic features  (positive/ negative/radical). 

enum GADataType { AdditivityType, CorrectionsType};

typedef bool (*ThermoCorrelationPtr)(const Molecule &, double&, double&, double&, double&);

class ThermoGA
{
	protected:
		static std::map < std::string, HashGAMap* > AtomtypeGAMap;//this is the big std::map that keeps track of all groups starting with a given atomtype.
		static bool IncrementAdded(const Molecule&, int, const GAMap&, std::map<int, std::pair<int,int> >&, std::set<int>&, ThermoType, double&, double);
		static std::multimap <std::string, std::pair<ConstrPtr, GAdata> > Corrections;
		static std::multimap <std::string, std::pair<ConstrPtr, double> > BECorrections;
		static std::vector<std::pair<ConstrPtr, GAdata> > MolecularCorrections;
		static std::map <std::string, HashGAMap* > logPGAMap; //this keeps track of all logP groups
		static double getCorrectionsValue(const Molecule&, ThermoType, double);
		static bool AHashPairExists(const Molecule&, std::string, int, std::pair<int,int>&);// constructs reasonable hashstd::pairs and checks if they exist in HashGAMap.
		static bool calculateThermoProperty(const Molecule&, double&, ThermoType, std::map < std::string, HashGAMap* >&, double, bool);
		static bool CalculateIncrement(const Molecule&, std::string, int, std::map<int, std::pair<int,int> >&, std::set<int>&, ThermoType, std::map < std::string, HashGAMap* >&, double&, double);
		static double calculateTempCorrections(double, ThermoType);
		static double getSymmetryCorrection(const Molecule&, ThermoType);
		static void AddAdditivityData(std::string&, GAdata&, std::map < std::string, HashGAMap* >&);
		
		

	public:
		static void AddGA(std::string, GAdata);
		static void AddCorrections(std::string, ConstrPtr, GAdata);
		static void AddBECorrections(std::string, ConstrPtr, double);
		static void AddMolecularCorrections(ConstrPtr,GAdata);
		static bool calculateDeltaH(const Molecule&, bool, double, double&);//calculate deltaH formation of a molecule, true or false as to whether or not the molecule is site intermediate, at a given temperature, and store it in the fourth argument passed by reference
		static bool calculateDeltaS(const Molecule&, bool, double, double&);//calculate deltaS formation of a molecule, true or false as to whether or not the molecule is site intermediate, at a given temperature, and store it in the fourth argument passed by reference
		static bool calculateDeltaG(const Molecule&, bool, double, double&);
		static double calculateBE (const Molecule&);
		static bool calculateCp(const Molecule&, bool, double, double&);
		static bool ReadInputsFromFile(const char*,GADataType);//reads from a file the GA and corrections input.
		static void AddlogPGA(std::string, GAdata);
		static bool calculateLogP(const Molecule&, double&);//calculate logP of a molecule
		static ThermoCorrelationPtr CorrelationPtr;
		static bool HasCorrelationsForThermo;
		
};

class ThermoValues
{
	protected:
        std::map<int, double> EnthalpyFormation; // the ints refer to indices in AvailableTemps
        std::map<int, double> EntropyFormation;
        std::map<int, double> FreeEnergyFormation;
        std::map<int, double> CpMolecule;
        std::map<int, double> logP;//incomplete
		std::vector<double> AvailableTempsForEnthalpy;
		std::vector<double> AvailableTempsForEntropy;
		std::vector<double> AvailableTempsForFreeEnergy;
		std::vector<double> AvailableTempsForCp;
		std::vector<double> AvailableTempsForlogP;
	public:
		double getThermo (ThermoType, double);
		void setThermo(ThermoType, double, double);
		int isTempAvailable(ThermoType,double); //returns the index of the std::vector AvailableTemps; -1 if not available!
};

#endif
