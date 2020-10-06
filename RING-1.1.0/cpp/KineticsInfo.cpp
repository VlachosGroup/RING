#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <vector>
#include <cctype>
#include <queue>
#include <algorithm>
#include <set>
#include <cmath>
using namespace std;



#include "Classheader.h"
#include "AdditionalFunctions.h"
#include "StringRegistry.h"


KineticsInfo::KineticsInfo(vector<KineticParamPtr>& kinetics, generated_rxn& reaction, double GasConst, double T, double H, double S, double n, double d, int r, int p)
{
	KineticFunctions = &kinetics;
	rxn = &reaction;
	
	Rg = GasConst;
	Temp = T;
	RxnEnthalpy = H;
	RxnEntropy = S;
	del_n = n;
	siteDensity = d;
	firstGasPhaseReactant = r;
	firstGasPhaseProduct = p;
	calcCollisionFreq = true;

}

void KineticsInfo::setCollisionFreq(bool coll)
{
	calcCollisionFreq = coll;
}

KineticsInfo::KineticsInfo(vector<AltKineticParamPtr>& kinetics, vector<double>* params, generated_rxn& reaction, double GasConst, double T, double H, double S, double n)
{
	AltKineticFunctions = &kinetics;
	rxn = &reaction;
	Rg = GasConst;
	Temp = T;
	RxnEnthalpy = H;
	RxnEntropy = S;
	del_n = n;
	parameters = params;
	calcCollisionFreq = true;

}

bool KineticsInfo::getKineticParameters( double& PreExp, double& ActE, double& TempIndex, double& kinValue, bool& calcK, bool& usesBEP, bool& usesLFER, double& alpha, double& beta, bool& stick)
{
	//TODO:: The tasks involved in calculating parameters defined through reverse reaction have unnecessarily repetitive steps. This is done in purpose to be sure they are correct! 
	
	int rule = rxn->get_rule();
		
	vector<Molecule> reactants;
	vector<Molecule> products;



	for (int m_i = 0;m_i<rxn->number_pdcts();m_i++)
	{
		string* sptr = rxn->get_products(m_i);
		Molecule mol(*sptr, moleculesize(*sptr));
		mol.unique_smiles();
		products.push_back(mol);
		
	}
	
	for (int m_i = 0;m_i<rxn->number_reactants();m_i++)
	{
		string* sptr = rxn->get_reactants(m_i);
		Molecule mol(*sptr, moleculesize(*sptr));
		mol.unique_smiles();
		reactants.push_back(mol);
	}

	bool isBEP = false;
	
	if (!KineticFunctions->at(rule)(reactants,products,PreExp,TempIndex, ActE, calcK, isBEP, usesLFER, alpha,beta,stick,Temp))
	{
		cout<<"cannot set kinetics of reaction "<<endl;
		cout<<rxn->reactionstring()<<endl;
		return false;
	}

	

	//if (calcK || isBEP) RxnEnthalpy = Network->calculateThermoOfRxn(EnthalpyType,Network->AllReactions[i],Network->Temperature)*1000.0;

	if (isBEP || usesLFER)
	{
		if (calcK)
		{
			ActE = -alpha*RxnEnthalpy + beta;//we shouldn't do 1-alpha yet! because kinValue is corrected with eq const K
			if (usesLFER) 
			{
				
				ActE+= (alpha-1)*calculateBEofMolecules(products); //because I need to add the total BE of reactants of the reverse reaction
			
			}
		}
		else
		{
			ActE = alpha*RxnEnthalpy + beta;
			if (usesLFER) ActE+= (alpha-1)*calculateBEofMolecules(reactants);
		}
	}
	//Note up to this point, for those reactions that are calculated with reverse steps, we cacluate kinetics of that reverse step. 
		

	if (stick && calcCollisionFreq)
	{
		//PreExp at this point only contains sticking coefficient and any conversion from atm to volUnits/MolUnits (taken care of by translation) and any other correction factors input by the user (like site density for CHEMKIN formats)

		//PreExp is StickingFactor*Asite * 1/kb1 *sqrt(kb2/(2*pi*Temp)) * 1/sqrt(m)  * 1/10 
		//this is now in (volUnits/molUnits)^2/TimeUnits

		//Asite is in cm^2/site, kb1 is boltzmann constant in l.atm/K, kb2 is boltzmann constant in m^2Kg/s^2/K,
		// m is moleclar mass in Kg, and 1/10 is a conversion factor so that the product is in 1/s/atm
		double Asite = 1/siteDensity; //TODO - this needs to be input!
		double kb1 = 0.0821/6.023e23;
		double kb2 = 1.38e-23;

		//sticking coefficients require the gas phase reactant information. These sticking coefficients are to be applied only to cases where there is one gasphase reactant only
		double m=0.0;
		if (calcK)
			m = products[firstGasPhaseProduct].MolecularWeight();
		else
			m = reactants[firstGasPhaseReactant].MolecularWeight();
		
		m=m/(1000.0*6.023e23);

		

		//cout<<PreExp<<" "<<Asite<<" "<<(1/kb1)<<" "<<sqrt(kb2/(2*3.14*Temp))<<"  "<<1/sqrt(m)<<endl;

		PreExp = PreExp*Asite*(1/kb1)*sqrt(kb2/(2*3.14*Temp))*1/sqrt(m)*0.1;
		
	}
	double RT = RCONST(Temp*8.314);
	double EaUponRT = ActE/RT;
	kinValue = PreExp*pow(Temp, TempIndex)*exp(-EaUponRT);


	
	if (calcK)// now, if the specification is in terms of reverse step, we calculate the actual kinetic value
	{
		//double del_n = Network->getDeltaNGasPhase(Network->AllReactions[i]);

		//double dS = Network->calculateThermoOfRxn(EntropyType,(*ReactionsPtr)[i],Temp);
		double deltaG = RxnEnthalpy - RxnEntropy*Temp;

		double K = exp(-deltaG/(RCONST(Temp*8.314)));
		K=K*pow(Rg*Temp,-del_n);//NOTE: Here, R depends on what value of Rg we use. 
		kinValue= kinValue*K;//(for reaction A--> B, kf = kr*K(A->B) and kinvalue will now be kr!

		ActE+=RxnEnthalpy; // we need to add reaction enthalpy to this-- this is because we need to also output ActE
		//if (stick) cout<<del_n<<"   "<<dS<<"  "<<exp(dS/RCONST(8.314))<<"  "<<pow(Rg*Temp,del_n)<<endl;
		PreExp = PreExp*exp(RxnEntropy/RCONST(8.314))*pow(Rg*Temp,-del_n); //for reaction A--> B, Af = Ar*exp(dS(A->B)/R) and PreExp will now be Ar! an RT^(-delN) is also included to correct the units! 



	}

	usesBEP = isBEP;
	return true;

}

bool KineticsInfo::getKineticParameters(double& kinValue, int& kIndex, double& kRelative, int& EaIndex, double& EaRelative,double& RefTemp, double& TempExp, bool& calcK)
{
	//THERE is a lot of repetition from the previous function - need to make it more efficient. 
	int rule = rxn->get_rule();
		
	vector<Molecule> reactants;
	vector<Molecule> products;

	for (int m_i = 0;m_i<rxn->number_pdcts();m_i++)
	{
		string* sptr = rxn->get_products(m_i);
		Molecule mol(*sptr, moleculesize(*sptr));
		mol.unique_smiles();
		products.push_back(mol);
		
	}
	
	for (int m_i = 0;m_i<rxn->number_reactants();m_i++)
	{
		string* sptr = rxn->get_reactants(m_i);
		Molecule mol(*sptr, moleculesize(*sptr));
		mol.unique_smiles();
		reactants.push_back(mol);
	}

	KineticsEstimates kEstimates = AltKineticFunctions->at(rule)(reactants,products);

	kIndex = kEstimates.kParam;
	kRelative = kEstimates.kRelative;
	EaIndex = kEstimates.EaParam;
	EaRelative = kEstimates.EaRelative;
	RefTemp = kEstimates.kRefTemp;
	calcK = kEstimates.isReverse;
	TempExp = kEstimates.TempExp;

	//NOTE teh kinValue is not calculated here - it's the parameter and relative info that's required here


	
	if (calcK)
	{
	
		double deltaG = RxnEnthalpy - RxnEntropy*Temp;

		double K = exp(-deltaG/(RCONST(Temp*8.314)));
		K=K*pow(Rg*Temp,-del_n);//NOTE: Here, R depends on what value of Rg we use. 
		kinValue= kinValue*K;//(for reaction A--> B, kf = kr*K(A->B) and kinvalue will now be kr!
		
		//updating kRelative with K. NOT changing EaRelative because K captures this! 
		//This is one place where this function differs from the previous function!
		kRelative = kRelative*K;
	}



	return true;



}
		
double KineticsInfo::calculateBEofMolecules(std::vector<Molecule> & mols)
{
	double TotalBE = 0.0;
	for (int i =0;i<mols.size();i++)
	{
		TotalBE+=ThermoGA::calculateBE(mols[i]);
	}
	
	return TotalBE;
}
		


