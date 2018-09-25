#include <string>
#include <vector>
#include <cmath>

#include "lumping.h"

using std::string; using std::vector; using std::pow;

LumpInfo::LumpInfo(int Molsize, int Hydrogens, bool Molshape, int brValue, int dist, string * molstr)
{
	size = Molsize;
	HydCount = Hydrogens;
	isCyclic = Molshape;
	BranchValue = brValue;
	ringleavesDistance = dist;
	MoleculeString = molstr;
	leaves = -1;
	PONAcharacteristic = -1;
	PONAElectronicValue = 0;
	DoubleBondCount = 0;
	TripleBondCount = 0;
	AromaticBondCount = 0;
	MaxRingSize = 0;
	RingCount = 0;
}

void LumpInfo::setMoleculeString(string* S)
{
	MoleculeString = S;
}

string LumpInfo::getMoleculeString()
{
	return *MoleculeString;
}

string* LumpInfo::getMolStringPtr()
{
	return MoleculeString;
}
int LumpInfo::getSize()
{
	return size;
}

int LumpInfo::getHydrogens()
{
	return HydCount;
}

int LumpInfo::getBranchValue()
{
	return BranchValue;
}
int LumpInfo::getringDistanceValue()
{
	return ringleavesDistance;
}

void LumpInfo::setBranchValue(int b)
{
	BranchValue  = b;
}

void LumpInfo::setringDistanceValue(int b)
{
	ringleavesDistance = b;
}

int LumpInfo::getleaves()
{
	return leaves;
}

void LumpInfo::setleaves(int i)
{
	leaves = i;
}

void LumpInfo::setRank(int i)
{
	rank = i;
}

int LumpInfo::getRank()
{
	return rank;
}

void LumpInfo::setPONAcharacteristic(int i)
{
	PONAcharacteristic = i;
}

int LumpInfo::getPONAcharacteristic()
{
	return PONAcharacteristic;
}

int LumpInfo::getPONAElectronicValue()
{
	return PONAElectronicValue;
}

void LumpInfo::setPONAElectronicValue(int v)
{
	PONAElectronicValue = v;
}

void LumpInfo::setBondCounts(int d, int t, int a)
{
	DoubleBondCount = d;
	TripleBondCount = t;
	AromaticBondCount = a;
}

int LumpInfo::getDoubleBondCount()
{
	return DoubleBondCount;
}

int LumpInfo::getTripleBondCount()
{
	return TripleBondCount;
}

int LumpInfo::getAromaticBondCount()
{
	return AromaticBondCount;
}

int LumpInfo::getMaxRingSize()
{
	return MaxRingSize;
}

int LumpInfo::getRingCount()
{
	return RingCount;
}

void LumpInfo::setMaxRingSize(int m)
{
	MaxRingSize = m;
}

void LumpInfo::setRingCount(int c)
{
	RingCount = c;
}

int LumpInfo::getScore()
{
	int MaxRingSizeValue = 0;
	if (MaxRingSize>0) MaxRingSizeValue = MaxRingSize - 2;//ring will always be greater than size 2.
	int AromaticBondCountValue = 0;
	if (AromaticBondCount>0) AromaticBondCountValue = AromaticBondCount - 4; //there will always be 5 or more aromatic bonds in an aromatic molecule

	return (int)(pow(3.0, MaxRingSizeValue) * pow(5.0,DoubleBondCount) * pow(7.0, RingCount) * pow (11.0, TripleBondCount) * pow(13.0,AromaticBondCountValue));
}


LumpingStrategy::LumpingStrategy(bool lump)
{
	toLump = lump;
	chainParameter = 0;
	ringParameter = 0;
	paraffinParameter = -1;
	olefinParameter = -1;
	aromaticsParameter =-1;

	HasSetFunctionalLumpingConstraints = false;
	HasSetParaffinConstraints = false;
	HasSetOlefinConstraints = false;
	HasSetNaphthenicsConstraints = false;
	HasSetAromaticsConstraints = false;
	HasSetMoreLumpingConstraints = false;

}

LumpingStrategy::LumpingStrategy()
{
	toLump = false;
	chainParameter = 0;
	ringParameter = 0;
	paraffinParameter = -1;
	olefinParameter = -1;
	aromaticsParameter = -1;

}

void LumpingStrategy::setParameters(int a, int b, int c, int d, int e, int f)
{
	chainParameter = a;
	ringParameter = b;
	paraffinParameter = c;
	olefinParameter = d;
	naphthenicsParameter = e;
	aromaticsParameter = f;

}

bool LumpingStrategy::shoudLump()
{
	return toLump;
}

int LumpingStrategy::getchainParameter()
{
	return chainParameter;
}

int LumpingStrategy::getringParameter()
{
	return ringParameter;
}

int LumpingStrategy::getParaffinParameter()
{
	return paraffinParameter;
}

int LumpingStrategy::getOlefinParameter()
{
	return olefinParameter;
}

int LumpingStrategy::getAromaticsParameter()
{
	return aromaticsParameter;
}

int LumpingStrategy::getNaphthenicsParameter()
{
	return naphthenicsParameter;
}

vector<int> LumpingStrategy::getMoreLumpingParameter()
{
	return MoreLumpingParameter;
}

void LumpingStrategy::setParaffinConstraints(ConstrPtr CP)
{
	ParaffinConstrPtr=CP;
	HasSetParaffinConstraints = true;
}

void LumpingStrategy::setOlefinConstraints(ConstrPtr CP)
{
	OlefinConstrPtr=CP;
	HasSetOlefinConstraints = true;
}

void LumpingStrategy::setNaphthenicsConstraints(ConstrPtr CP)
{
	NaphthenicsConstrPtr=CP;
	HasSetNaphthenicsConstraints  = true;
}

void LumpingStrategy::setAromaticsConstraints(ConstrPtr CP)
{
	AromaticsConstrPtr=CP;
	HasSetAromaticsConstraints = true;
}

void LumpingStrategy::setMoreLumpingConstraints(ConstrPtr CP, int i)
{
	MoreLumpingConstrPtr.push_back(CP);
	HasSetMoreLumpingConstraints = true;
	MoreLumpingParameter.push_back(i);
}

void LumpingStrategy::setFunctionalLumpingConstraints(ConstrPtr CP)
{
	HasSetFunctionalLumpingConstraints = true;
	FunctionalLumpingConstrPtr = CP;
}

ConstrPtr LumpingStrategy::getParaffinConstraints()
{
	return ParaffinConstrPtr;
}

ConstrPtr LumpingStrategy::getOlefinConstraints()
{
	return OlefinConstrPtr;
}

ConstrPtr LumpingStrategy::getNaphthenicsConstraints()
{
	return NaphthenicsConstrPtr;
}

ConstrPtr LumpingStrategy::getAromaticsConstraints()
{
	return AromaticsConstrPtr;
}

vector<ConstrPtr> LumpingStrategy::getMoreLumpingConstraints()
{
	return MoreLumpingConstrPtr;
}

ConstrPtr LumpingStrategy::getFunctionalLumpingConstraints()
{
	return FunctionalLumpingConstrPtr;
}

bool LumpingStrategy::isThereParaffinConstraints()
{
	return HasSetParaffinConstraints;
}

bool LumpingStrategy::isThereOlefinConstraints()
{
	return HasSetOlefinConstraints;
}

bool LumpingStrategy::isThereNaphthenicsConstraints()
{
	return HasSetNaphthenicsConstraints;
}

bool LumpingStrategy::isThereAromaticsConstraints()
{
	return HasSetAromaticsConstraints;
}

bool LumpingStrategy::isThereMoreLumpingConstraints()
{
	return HasSetMoreLumpingConstraints;
}

bool LumpingStrategy::isThereFunctionalConstraints()
{
	return HasSetFunctionalLumpingConstraints;
}

LumpedReaction::LumpedReaction(vector<int> r, vector<int> p, int ru)
{
	reactantsLumpSet = r;
	productsLumpSet = p;
	ruleIndex = ru;

}

vector<int> LumpedReaction::getReactantLumps() const
{
	return reactantsLumpSet;
}

vector<int> LumpedReaction::getProductLumps() const
{
	return productsLumpSet;
}

int LumpedReaction::getRule()const
{
	return ruleIndex;
}




