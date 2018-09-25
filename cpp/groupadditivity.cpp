#include <iostream>
#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
#include <map>
#include <vector>
#include <utility>
//#include <cctype>
//#include <queue>
//#include <algorithm>
#include <set>

//#include <stdio.h>
//#include <stdlib.h>
#include <cmath>

//#include "common.h"
#include "additionalfunc.h"
//#include "stringreg.h"
//#include "clonable.h"
//#include "element.h"
//#include "atom.h"
//#include "singleatom.h"
//#include "compositeatom.h"
//#include "atomcontainer.h"
//#include "molecule.h"
#include "patternmatch.h"
#include "substructure.h"
#include "automorphs.h"
#include "groupadditivity.h"

using std::vector; using std::string; using std::map; using std::multimap;
using std::pair; using std::set; 
using std::log; using std::pow; //using std::multimap;
using std::cout; using std::endl; using std::ifstream; //using std::endl; 

vector<double> GATempDiscretePts::TempPts;

vector<double> GATempDiscretePts::getTemperaturePoints()
{
	return GATempDiscretePts::TempPts;
}

void GATempDiscretePts::InsertTempPoints(double Pt)
{
	GATempDiscretePts::TempPts.push_back(Pt);
}

GAdata::GAdata(vector<double> data)
{
	if (data.size()!=9)
	{
		cout<<"there are not 9 data points here, assuming the remaining to be zero!"<<endl;
		data.resize(9,0.0);

	}

	Enthalpy = data[0];
	Entropy = data[1];
	for (int i=2;i<data.size();i++)
		Cp.push_back(data[i]);

	CalculateAandBForEachRange();
}

double GAdata::getEnthalpy()
{
	return Enthalpy;
}

double GAdata::getEntropy()
{
	return Entropy;
}

void GAdata::setEnthalpy(double H)
{
	Enthalpy = H;
}

void GAdata::setEntropy(double S)
{
	Entropy = S;
}

void GAdata::setCp(int temp, double value)
{
	if (CpMap.count(temp)==0)
		CpMap.insert(pair<int,double>(temp,value));
}


void GAdata::setlogP(double value)
{
	logP = value;
}

double GAdata::getCpAt(int temp)const
{
	
	if (CpMap.count(temp)==0)
	{
		return getCpAtArbitraryTemp(temp);
	}
	else return CpMap.find(temp)->second;
		
}

double GAdata::getCpAtArbitraryTemp(int temp) const
{
	double Cplower = 0.0;
	int Tlower = 0;
	double CpHighest = 0.0;
	map<int,double>::const_iterator it;
	for ( it= CpMap.begin();it!=CpMap.end();it++)
	{
		if (it->first<=temp)
		{
			Cplower = it->second;
			Tlower = it->first;
		}
		CpHighest=it->second;
		if (it->first > temp)
			return Cplower + ((it->second-Cplower)/(it->first-Tlower))*(temp-Tlower);
	}
	if (CpMap.size()>0)
		return CpHighest;
	else
	{
		cout<<"no value available for Cp at temperature "<<temp<<"; returning zero!"<<endl;
		return 0.0;
	}
}
	


void GAdata::PrepareGA()
{
	CalculateAandBForEachRange();
}

/*
double GAdata::getCpAtEndPoints(double Temp)
{
	if (Temp==298. || Temp==300.)return Cp300;
	else if (Temp== 400.) return Cp400;
	else if (Temp== 500.) return Cp500;
	else if (Temp== 600.) return Cp600;
	else if (Temp== 800.) return Cp800;
	else if (Temp== 1000.) return Cp1000;
	else if (Temp== 1500.) return Cp1500;
	else 
	{
		cout<<"oops! that's not one of the end points of the range"<<endl;
		return 0.0;
	}
}*/

bool GAdata::getdata(ThermoType t, double Temp, double& dataVal) const
{
	try{
		if (t== EnthalpyType)
		{
			//cout << "Enthalpy1: " << Enthalpy << endl;
			//cout <<"dataVal1: " << dataVal <<endl;
			//cout <<"Temp: "<<Temp<<endl;
			//cout<< "Corr: "<<getTempCorrection(t, Temp)/1000.<<endl;
			dataVal = Enthalpy; //+ getTempCorrection(t, Temp)/1000.;
			//cout << "Enthalpy2: " << Enthalpy << endl;
			//cout <<"dataVal2: " << dataVal <<endl;
			return true;
		}
		else if (t==EntropyType)
		{
			cout << "Entropy1: " << Entropy << endl;
			dataVal = Entropy; //+ getTempCorrection(t,Temp);
			cout << "Entropy2: " << Entropy << endl;
			return true;
		}
		else if (t==CpType)
		{
			dataVal = getCpAt(Temp);//TODO - allow for double values for temperature. This still works because I never ask for cp at non integer values of temperature! 
			return true;
		}
		else if (t==logPType)
		{
			dataVal = logP;
			return true;
		}
		else 
		{
			cout<<"requesting a non-existing GA data type!"<<endl;
			return false;
		}
	}
	catch (int er)
	{
		cout<<"cannot get the corresponding thermo data"<<endl;
		return false;
	}
}


double GAdata::getTempCorrection(ThermoType type, double Temp) const
{
	double TempCorrection = 0.0;
	double ArbitrarySmallNum = 0.0001;
	vector<int>	TemperatureRanges;
	map<int,double>::const_iterator it;

	for (it =CpMap.begin();it!=CpMap.end();it++)
		TemperatureRanges.push_back(it->first);
	if (Temp<double(TemperatureRanges.front())-ArbitrarySmallNum)
	{
		cout<<"Temperature sought ("<<Temp<<") is less than what Cp has been specified for ("<<TemperatureRanges.front()<<")"<<endl;
		throw 1;
	}
	if (Temp>=double(TemperatureRanges.back())- ArbitrarySmallNum)
	{
		if (type == EnthalpyType)
			return getCpAt(TemperatureRanges.back())*(Temp-double(TemperatureRanges.back()));//Cp[6] refers to Cp at 1500 K
		else if (type == EntropyType) 
			return getCpAt(TemperatureRanges.back())*log(Temp/double(TemperatureRanges.back()));//if entropy return the corresponding correction term
		else 
		{
			cout<<"wrong thermo type!"<<endl;
			return 0.0;
		}
	}
	else 
	{
		for (int i=0;i<TemperatureRanges.size()-1;i++)
		{
			if (Temp>=double(TemperatureRanges[i])- ArbitrarySmallNum && Temp<double(TemperatureRanges[i+1])+ ArbitrarySmallNum)
			{
				if (Temp>double(TemperatureRanges[0]) + ArbitrarySmallNum)
				{
					if (type == EnthalpyType)
					{
						TempCorrection+= A[i]*(pow(Temp,2)-pow(double(TemperatureRanges[i]),2))/2 + B[i]*(Temp-TemperatureRanges[i]);
						TempCorrection+= getTempCorrection(EnthalpyType,double(TemperatureRanges[i]));
					}
					else if (type == EntropyType) 
					{
						TempCorrection+=A[i]*(Temp-TemperatureRanges[i]) + B[i]*log(Temp/TemperatureRanges[i]);
						TempCorrection+= getTempCorrection(EntropyType,double(TemperatureRanges[i]));
					}
					else
					{
						cout<<"wrong thermo type"<<endl;
						throw 1;
					}
				}//if not TempCorrection is zero -- this would break the recursive loop
				break;
			}
		}

		return TempCorrection;
	}
	
}

void GAdata::CalculateAandBForEachRange()
{
	vector<int> TemperatureRanges;
	map<int,double>::iterator it;
	for (it =CpMap.begin();it!=CpMap.end();it++)
		TemperatureRanges.push_back(it->first);
	for (int i=0;i<TemperatureRanges.size()-1;i++)
	{
		A.push_back((CpMap[TemperatureRanges[i+1]]-CpMap[TemperatureRanges[i]])/(TemperatureRanges[i+1]-TemperatureRanges[i]));
		B.push_back((CpMap[TemperatureRanges[i]]*TemperatureRanges[i+1]-CpMap[TemperatureRanges[i+1]]*TemperatureRanges[i])/(TemperatureRanges[i+1]-TemperatureRanges[i]));
	}
}


map < string, HashGAMap* > ThermoGA::AtomtypeGAMap;
multimap <string, pair<ConstrPtr, GAdata > > ThermoGA::Corrections;
map < string, HashGAMap* > ThermoGA::logPGAMap;
multimap <string, pair<ConstrPtr, double > > ThermoGA::BECorrections;
vector<pair<ConstrPtr, GAdata> > ThermoGA::MolecularCorrections;
bool DummyTrueFunction(const Molecule&)
{
	return true;
}

void ThermoGA::AddAdditivityData(string& Str, GAdata& f, map < string, HashGAMap* >& hashmap)
{
	Substructure S(Str, patternsize(Str));

	pair<int,int> HashPair;
	HashPair.first = S.getGroupHash(0);
	HashPair.second = S.getNNElecHash(0)*S.getNNDoubleTripleBondFactors(0);
	string AtomType = S.getatomtype(0);
	if (hashmap.count(AtomType)==0)//if AtomType not found here, it means no group with this first atomtype has been defined yet
	{
		HashGAMap* H = new HashGAMap(); //creating a new HashGAMap here
		GAMap g;
		g.insert(pair<string,GAdata>(Str,f));//setting up a GAMap to store the group and its value. 

		H->insert(pair<pair<int,int>, GAMap>(HashPair,g));

		hashmap.insert(pair<string, HashGAMap*>(AtomType, H));
	}
	else
	{
		if(hashmap[AtomType]->count(HashPair)==0)
		{
			GAMap g;
			g.insert(pair<string,GAdata>(Str,f));

			hashmap[AtomType]->insert(pair<pair<int,int>, GAMap>(HashPair,g));

		}
		else
		{
			(*hashmap[AtomType])[HashPair].insert(pair<string, GAdata>(Str,f));
		}
	}
}


void ThermoGA::AddGA(string Str, GAdata f)
{
	AddAdditivityData(Str, f,AtomtypeGAMap);
	
}

void ThermoGA::AddlogPGA(string Str, GAdata f)
{
	AddAdditivityData(Str, f, logPGAMap);
}

bool ThermoGA::calculateThermoProperty(const Molecule & mol, double& value, ThermoType type, map < string, HashGAMap* >& hashmap, double Temp, bool isIntermediate)
{
	value = 0.0;
	
	// if (HasCorrelationsForThermo)
	// {
		// double EnthalpyA = 0.0;
		// double EnthalpyB = 0.0;
		// double EntropyA = 0.0;
		// double EntropyB = 0.0;

		// if ((CorrelationPtr)(mol,EnthalpyA, EnthalpyB, EntropyA, EntropyB))
		// {

			// if (type ==EnthalpyType)
			// {
				// value = EnthalpyA*mol.totalAtomsOfType("C") + EnthalpyB;
				// return true;
			// }
			// if (type ==EntropyType)
			// {
				// value = EntropyA*mol.totalAtomsOfType("C") + EntropyB;
				// return true;
			// }

			// return true; //if not enthalpy or entropy type, default set things to zero FOR NOW
		// }
	
	// }


	set<int>AtomsToBeCovered;//set of atoms to be checked with atom-centered groups

	map<int, pair<int,int> > AtomHashPairMap;
	for (int i=0;i<mol.getsize();i++)
	{
		AtomsToBeCovered.insert(i);//setting it initially with all atoms
		pair<int,int> P;
		P.first = mol.getGroupHash(i);
		P.second =mol.getNNElecHash(i)*mol.getNNDoubleTripleBondFactors(i);
		//cout << "first: " << P.first << endl;
		//cout << "second :" << P.second << endl;
		AtomHashPairMap[i] = P; // calculating the atom-hash pair mapping for each atom upfront
	}

	//mol.print_adjacency_list();

	try{
		cout<<"In here1"<<endl;
		while (!AtomsToBeCovered.empty())
		{
			int FirstAtom = (*AtomsToBeCovered.begin());// the first atom of the set
		
			string S = mol.getatomtype(FirstAtom);// atomtype
			//AtomsToBeCovered.erase(atom);
			//cout << "S: " << S << endl;
			//system("pause");
	
			if (hashmap.count(S)>0)//if group additivity exists for the atomtypes
			{
				set<int> AtomsCoveredNow;
				cout<<"value 1: "<<value<<endl;
				if (CalculateIncrement(mol, S, FirstAtom, AtomHashPairMap, AtomsCoveredNow, type, hashmap, value, Temp))
				{
					set<int>::iterator it;
					for (it = AtomsCoveredNow.begin(); it!=AtomsCoveredNow.end();it++)
						AtomsToBeCovered.erase(*it);
				}
				else
				{
					cout<<"group corresponding to"<<mol.AtomCenteredGroupForGA(FirstAtom)<<" was not declared for "<<S<<endl;
					cout<<"Note: d or t just refers to double or triple bonded atoms, the other atom of the bond can be anything!"<<endl;
					cout<<"Exiting"<<endl;
					return false;
				}
			}
			else
			{
				cout<<"groups for "<<S<<" was not declared at all!"<<endl;
				cout<<"Exiting"<<endl;
				return false;
			}
			cout<<"value now is  "<<value<<endl;
			
		}
	

		cout<<"value is "<<value<<endl;
		if (type!=logPType)
		{
			value+=getCorrectionsValue(mol,type, Temp);
			cout<<" value after corrections is "<<value<<endl;
			
			if (!isIntermediate) //allowing symmetry correction to apply only when it's an intermediate!! 
			value-=getSymmetryCorrection(mol,type);
			cout<<"value after symmetry correction is "<<value<<endl;
			
		}
		system("pause");
	}
	catch (int er)
	{
		cout<<"Exiting"<<endl;
		return false;
	}
	return true;
}

bool ThermoGA::CalculateIncrement(const Molecule & mol, std::string S, int FirstAtom, map<int, pair<int,int> >& AtomHashPairMap, std::set<int> & AtomsAlreadyCovered, ThermoType valueType, map < string, HashGAMap* >& hashmap, double& value, double Temp)
{
	pair<int,int> P;

	P = AtomHashPairMap[FirstAtom];

	//P.first = mol.getGroupHash(FirstAtom);
	//P.second =mol.getNNElecHash(FirstAtom);

	int NNDoubleBonds = mol.getNNDoubleBondCount(FirstAtom);
	int NNTripleBonds = mol.getNNTripleBondCount(FirstAtom);

	try{
		while (NNDoubleBonds+NNTripleBonds >=0)
		{
			//int NNFactor = (int)(pow(11.0,NNDoubleBonds)*pow(13.0,NNTripleBonds));//11 and 13 because NNElecHash takes prime numbers up to 7
	
			//P.second=P.second*NNFactor;
	
			if ((hashmap[S]->count(P)>0))
			{
				if (IncrementAdded(mol, FirstAtom, (*AtomtypeGAMap[S])[P], AtomHashPairMap, AtomsAlreadyCovered, valueType, value, Temp))	
					return true;
			}
		
			//P.second=P.second/NNFactor;
			if (NNDoubleBonds>0)
			{
				NNDoubleBonds--;
				P.second = P.second/11;
			}
			else if (NNTripleBonds>0)
			{
				NNTripleBonds--;
				P.second = P.second/13;
			}
			else NNDoubleBonds--;//making it -1 essentially so the sum is also -1;
		
		}


	}catch (int er){
		throw;
	}
	cout<<"None of the groups match for atom "<<FirstAtom<<endl;
	cout<<P.first<<"  "<<P.second<<endl; 
	//throw 1;
	return false;
}


bool ThermoGA::IncrementAdded(const Molecule & mol, int atomIndex, const GAMap& g, map<int, pair<int,int> >& AtomHashPairMap, set<int>& AtomsCovered, ThermoType valueType, double& thermoValue, double Temp)
{
	GAMap::const_iterator it;

	for (it=g.begin();it!=g.end();it++)
	{
		Substructure S(it->first, patternsize(it->first));

		Patternmatch P(mol,S,0);
		if (P.GetDistinctMatches()>0)
		{
			set<int> atomsCoveredByFirstSubAtom = P.AtomsCoveredBySubstr(0); // the set of atoms matching the first atom of the pattern described by the substructure S
			if (atomsCoveredByFirstSubAtom.count(atomIndex)>0)//if the current atom is covered
			{
				//AtomsCovered.insert(atomsCoveredByFirstSubAtom.begin(),atomsCoveredByFirstSubAtom.end());
				set<int>::iterator set_it;

				//cout<<" atom hash pair "<<AtomHashPairMap[atomIndex].second<<endl;
				for (set_it = atomsCoveredByFirstSubAtom.begin(); set_it!=atomsCoveredByFirstSubAtom.end();set_it++)
				{
					if (AtomHashPairMap[*set_it] == AtomHashPairMap[atomIndex] )//checking if the hashpair of the matched atom is equal to that of the first atom (atomindex)!
					{
						AtomsCovered.insert(*set_it);//add the atoms into teh set of newly covered atom iff it has the same hash pair as the first atom!
						//cout<<*set_it<<endl;
					}
				}

				double dataVal = 0.0;
				if (it->second.getdata(valueType,Temp,dataVal))
				{
					cout << "thermovalue1: " << thermoValue << endl;
					cout<<"AtomsSize: "<<AtomsCovered.size()<<endl;
					cout<<"Dataval: "<<dataVal<<endl;
					thermoValue+= AtomsCovered.size()*dataVal;
					cout << "thermovalue2: " << thermoValue << endl;
					
					return true;
				}
				else throw 1;
			}
			
		}
	}
	
	system("pause");
	return false;
}

void ThermoGA::AddCorrections(string S, ConstrPtr Cp, GAdata value)
{
	Corrections.insert(pair<string,pair<ConstrPtr,GAdata> >(S,pair<ConstrPtr, GAdata> (Cp,value)));
}

void ThermoGA::AddBECorrections(string S, ConstrPtr Cp, double value)
{
	BECorrections.insert(pair<string,pair<ConstrPtr,double> >(S,pair<ConstrPtr, double> (Cp,value)));
}


double ThermoGA::getCorrectionsValue(const Molecule& m, ThermoType valueType, double Temp)
{
	double correction = 0.0;
	multimap<string,pair<ConstrPtr,GAdata> >::iterator it;

	for (it = Corrections.begin();it!=Corrections.end();it++)
	{
		if (it->second.first(m))//if the constraints are satisfied
		{
			Substructure S(it->first, patternsize(it->first));

			Patternmatch P(m,S,0);
			pair<int,int> DistinctMatches = P.GetDistinctAllAndRingMatches();
			double data = 0.0;
			if (DistinctMatches.first > 0)
			{
				if (it->second.second.getdata(valueType,Temp, data))
				{
					correction+=(DistinctMatches.first)*data;
					//cout<<"new intermediate value of correction is "<<correction<<endl;
					//cout<<"this was after correcting for "<<it->first<<endl;
					//if (DistinctMatches.first - DistinctMatches.second >0)cout<<it->first<<"  "<<data<<endl;
				}
				else throw 1;
			}
		}
	}

	for (int i =0;i<MolecularCorrections.size();i++)
	{
		if ((MolecularCorrections[i].first)(m))
		{
			double data = 0.0;
			if (MolecularCorrections[i].second.getdata(valueType,Temp,data))
				correction+=data;
		}
	}
	cout<<"correction is "<<correction<<endl;
	return correction;
}

double ThermoGA::calculateBE(const Molecule & m)
{
	double BEcorrection = 0.0;
	multimap<string,pair<ConstrPtr,double> >::iterator it;

	for (it = BECorrections.begin();it!=BECorrections.end();it++)
	{
		if (it->second.first(m))//if the constraints are satisfied
		{
			Substructure S(it->first, patternsize(it->first));

			Patternmatch P(m,S,0);
			pair<int,int> DistinctMatches = P.GetDistinctAllAndRingMatches();
			double data = 0.0;
			if (DistinctMatches.first > 0)
				BEcorrection+=(DistinctMatches.first)*it->second.second;
		}
	}
	//cout<<"correction is "<<correction<<endl;
	return BEcorrection;


}

bool ThermoGA::ReadInputsFromFile(const char * filename, GADataType type)
{
	string line;
	
	int linenumber = 0;

	ifstream file(filename);
	string thermo_strfile(filename);
	cout<<"reading thermochemistry values from "<<thermo_strfile<<endl;

	if (file.is_open())
	{
		while (file.good())
		{
			getline(file,line);
			linenumber++;
			vector<string> tokens;
			if(TokenizeIntoWords(line,0,tokens))
			{
				/*if (tokens.size()!=10)
				{
					cout<<"error reading "<<thermo_strfile<<endl;
					cout<<"unexpected number of words in  line "<<linenumber<<" -- expected 10, but found "<<tokens.size()<<endl;
					return false;
				}*/
				
				vector<double> dataValues;
				for (int i =1;i<tokens.size();i++)
				{
					double dataPoint = 0.0;
					if (tokens[i].compare("*")!=0)
						dataPoint = StringToDouble(tokens[i]);
					dataValues.push_back(dataPoint);
				}

						
				GAdata data(dataValues);
				if (type == AdditivityType)ThermoGA::AddGA(tokens[0],data);
				else if (type == CorrectionsType)ThermoGA::AddCorrections(tokens[0],&DummyTrueFunction, data);			
				else 
				{
					cout<<"wrong choice: "<<type<<" should be Additivity (0) or Corrections (1)"<<endl;
					return false;
				}
			}
			else 
			{
				cout<<"cannot tokenize the line "<<linenumber<<endl;
				return false;
			}
		}
		return true;
	}
	else
	{
		cout<<"file not open"<<endl;
		return false;
	}
	
}

bool ThermoGA::calculateDeltaH(const Molecule& mol, bool isIntermediate, double Temp, double& value)
{
	return calculateThermoProperty(mol,value,EnthalpyType, AtomtypeGAMap, Temp, isIntermediate);
}

bool ThermoGA::calculateDeltaS(const Molecule& mol, bool isIntermediate, double Temp, double& value)
{
	return calculateThermoProperty(mol,value,EntropyType, AtomtypeGAMap, Temp, isIntermediate);
}

bool ThermoGA::calculateDeltaG(const Molecule& mol, bool isIntermediate, double Temp, double& value)
{
	double deltaH = 0.0;
	double deltaS = 0.0;

	if (ThermoGA::calculateDeltaH(mol,isIntermediate,Temp,deltaH) && ThermoGA::calculateDeltaS(mol,isIntermediate,Temp,value))
	{
		value = deltaH - Temp*deltaS;
		return true;
	}
	return false;
}

bool ThermoGA::calculateCp(const Molecule& mol, bool isIntermediate, double Temp, double& value)
{
	return calculateThermoProperty(mol,value,CpType,AtomtypeGAMap, Temp, isIntermediate);
}

bool ThermoGA::calculateLogP(const Molecule & mol, double & value)
{
	return calculateThermoProperty(mol,value,logPType, logPGAMap, 298.0, true);
}


double ThermoGA::getSymmetryCorrection(const Molecule& mol, ThermoType type)
{
	cout<<mol.moleculestring()<<"  "<<Automorphs(mol).SymmetryNumber()<<"  "<<mol.NumberChiralAtoms()<<endl;
	cout<<"symmetry number/ opt isomers "<<Automorphs(mol).SymmetryNumber()<<"  "<<mol.OptIsomers()<<endl;
	if (type == EntropyType)
		return 8.314*log(Automorphs(mol).SymmetryNumber()/double(mol.OptIsomers()));
	else return 0.0;
}



double ThermoValues::getThermo(ThermoType type, double Temp)
{
	int index = isTempAvailable(type, Temp);
	if (index ==-1)
	{
		cout<<"invalid temp requested for thermo, returning zero!"<<endl;
		return 0.0;
	}

		
	if (type == EnthalpyType) return EnthalpyFormation[index];
	else if (type ==EntropyType)return EntropyFormation[index];
	else if (type ==CpType) return CpMolecule[index];
	else if (type ==FreeEnergyType) return FreeEnergyFormation[index];
	else if (type == logPType) return logP[index];
	else {
		cout<<"invalid type requested for thermo, returning zero!"<<endl;
		return 0.0;
	}
}

void ThermoValues::setThermo(ThermoType type, double Temp, double value)
{
	
	if (isTempAvailable(type,Temp) !=-1)
	{
		cout<<"the values for this temp already exists, ignoring this input"<<endl;
	}

	
	if (type == EnthalpyType)
	{
		AvailableTempsForEnthalpy.push_back(Temp);
		EnthalpyFormation[AvailableTempsForEnthalpy.size()-1] = value;
	}
	else if (type ==EntropyType)
	{
		AvailableTempsForEntropy.push_back(Temp);
		EntropyFormation[AvailableTempsForEntropy.size()-1] = value;
		
	}
	else if (type ==CpType)
	{
		
		AvailableTempsForCp.push_back(Temp);
		CpMolecule[AvailableTempsForCp.size()-1] = value;
	}
	else if (type ==FreeEnergyType)
	{
		AvailableTempsForFreeEnergy.push_back(Temp);
		FreeEnergyFormation[AvailableTempsForFreeEnergy.size()-1] = value;
		
	}
	else if (type == logPType)
	{
		AvailableTempsForlogP.push_back(Temp);
		logP[AvailableTempsForlogP.size()-1] = value;
		
	}
}

int ThermoValues::isTempAvailable(ThermoType type, double Temp)
{
	int Index = -1;
	vector<double> temp;

	if (type == EnthalpyType)
		temp = AvailableTempsForEnthalpy;
	else if (type == EntropyType)
		temp = AvailableTempsForEntropy;
	else if (type == FreeEnergyType)
		temp = AvailableTempsForFreeEnergy;
	else if (type == CpType)
		temp = AvailableTempsForCp;
	else if (type == logPType) 
		temp = AvailableTempsForlogP;
	else{
		cout<<"wrong type specified for checking temp availability"<<endl;
	}

	for (int i = 0; i< temp.size();i++)
	{
		if (abs(temp[i] - Temp) < 0.01)
		{
			Index = i;
		}
	}

	return Index;
}







	

		


			






	
