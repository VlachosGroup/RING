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
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Classheader.h"
#include "additionalfunc.h"


Surfqsar::Surfqsar(Molecule * M, int EOUnits)
{
	mol = M;
	NumEOUnits = EOUnits;
}

double Surfqsar::calculateCMC()
{
	//cout<<"KHZero is "<<KHZeroCalc()<<endl;
	//cout<<"AICSecond is "<<AICSecond()<<endl;
	return -1.80 -0.567*KHZeroCalc() +1.054*AICSecond() + 7.5*RelativeOxygensCount();


}


double Surfqsar::calculateSurfTens()
{
	
	string molstr = mol->moleculestring();
	string ethoxyl = "CCO";
	for (int i=0;i<NumEOUnits;i++)
	{
		molstr.insert(molstr.find("O")+1,ethoxyl);
	}
	double dH = 0.0;
	ThermoGA::calculateDeltaH(Molecule(molstr,moleculesize(molstr)),false, 298.0,dH);
		
	
	return 11.63 + 0.6750*(NumEOUnits+1) + 0.6857*KHZeroCalc() - 0.01261*dH - 0.1387*KHZeroCalc()*(NumEOUnits+1);
}
double Surfqsar::CMCerror()
{
	return -0.16 -0.009*KHZeroCalc() + 0.048*AICSecond()+ RelativeOxygensCount();
}
double Surfqsar::KHZeroCalc()
{
	double DvSum=0.0;
	for (int i = 0;i< mol->getsize();i++)
	{
		
		Atom* atomptr = mol->getatom(i);
		if (atomptr->get_atomtype_name()=="C")
		{
			double num = atomptr->get_valency()-mol->getHydrogens(i);
			double den = atomptr->get_atomic_number()-atomptr->get_valency()-1;
			DvSum=DvSum + 1/sqrt(num/den);
		}
	}
	return DvSum;
}

double Surfqsar::RelativeOxygensCount()
{
	//return double(wholeMol->totalAtomsOfType("O"))/double(wholeMol->totalAtomsIncludingH());

	return double(NumEOUnits+1)/(double(mol->totalAtomsIncludingH()+ 7*NumEOUnits));
}
double Surfqsar::AICSecond()
{
	map<int,int> ValueFreqMapForC;
	map<int,int> ValueFreqMapForH;
	int totalAtomCount =0;
	for (int i=0;i<mol->getsize();i++)
	{
		Atom* atomptr = mol->getatom(i);
		int ValueSecond = 1;//for Carbons
		int ValueFirst = 1;//for Hydrogens
		if (atomptr->get_atomtype_name()=="C")
		{
			totalAtomCount++;
			totalAtomCount+=mol->getHydrogens(i);
			for (int j=0;j<mol->getNN(i);j++)
			{
				int ImmediateNeighbor = mol->get_adjacency(i,j);
				for (int k=0;k<mol->getNN(ImmediateNeighbor);k++)
				{
					int SecondNeighbor = mol->get_adjacency(ImmediateNeighbor, k);
					//a hack for now!
					if (mol->getatom(SecondNeighbor)->get_atomtype_name()=="H")
						ValueSecond = ValueSecond*3;
					if (mol->getatom(SecondNeighbor)->get_atomtype_name()=="C")
						ValueSecond = ValueSecond*5;
					if (mol->getatom(SecondNeighbor)->get_atomtype_name()=="O")
						ValueSecond = ValueSecond*7;
				}
				ValueSecond=ValueSecond*int(pow(double(3),double(mol->getHydrogens(ImmediateNeighbor))));
				
				if (mol->getatom(ImmediateNeighbor)->get_atomtype_name()=="H")
						ValueFirst = ValueFirst*3;
				if (mol->getatom(ImmediateNeighbor)->get_atomtype_name()=="C")
						ValueFirst = ValueFirst*5;
				if (mol->getatom(ImmediateNeighbor)->get_atomtype_name()=="O")
						ValueFirst = ValueFirst*7;

				ValueFirst=ValueFirst*int(pow(double(3),double(mol->getHydrogens(i)-1)));

			}
			if (ValueFreqMapForC.count(ValueSecond)!=0)
				ValueFreqMapForC[ValueSecond]++;
			else
				ValueFreqMapForC[ValueSecond]=1;
		
			if (ValueFirst!=0)
			{
				if (ValueFreqMapForH.count(ValueFirst)!=0)
					ValueFreqMapForH[ValueFirst]+=mol->getHydrogens(i);
				else
					ValueFreqMapForH[ValueFirst]=mol->getHydrogens(i);
			}


		}
	}

	double AICValue=0.0;

	map<int,int>::iterator it;
	
	for (it =ValueFreqMapForC.begin();it!=ValueFreqMapForC.end();it++)
	{
		
		
		double fraction =(double(it->second)/double(totalAtomCount));
		//cout<<fraction<<endl;	
		AICValue = AICValue - fraction*(log10(fraction)/log10(2.0));
	}
	//cout<<"  "<<AICValue<<endl;
	for (it =ValueFreqMapForH.begin();it!=ValueFreqMapForH.end();it++)
	{
		
		double fraction =(double(it->second)/double(totalAtomCount));
		//cout<<fraction<<endl;
		AICValue = AICValue - fraction*(log10(fraction)/log10(2.0));
	}
	//cout<<"   "<<AICValue<<endl;
	return AICValue;
	

}

double Surfqsar::KierShapeThird()
{
	double K = 0.0;

	double P = 0.0;
	double alpha = 0.0;
	/*for (int i=0;i<mol->getsize();i++)
	{
		if (mol->getatom(i)->get_atomtype_name()=="C")
		{
			for (int j=0;j<mol->getNN(i);j++)
			{
				P=P+ double(mol->getNN(i)-1)*double(mol->getNN(mol->get_adjacency(i,j))-1);
			}
			

			//if (mol->getatom(i)->get_atomtype_name()=="C")
			//	alpha = alpha+1 ;
			//if (mol->getatom(i)->get_atomtype_name()=="O")
			//	alpha = alpha - 0.2;//ratio of Oxygen'atomic radius to Carbon's
		}
		
		
	}*/

	int startingAtO = 0;
	for (int i=0;i<mol->getsize();i++)
	{
		for (int j=0;j<mol->getNN(i);j++)
		{
			int atomIJ = mol->get_adjacency(i,j);
			
			for (int k=0;k<mol->getNN(atomIJ);k++)
			{
				int atomIJK = mol->get_adjacency(atomIJ,k);
				if (atomIJK!=i)
				{
					P=P+mol->getNN(atomIJK)-1;
				}
			}
			
		}
	}


		


 
	
	P=P/2;
	//cout<<P<<"  "<<alpha<<endl;
	int N = (mol->getsize()-1);
	if (N%2!=0)
	{
		//cout<<"here"<<endl;
		K = (N + alpha -1)*pow((double(N) + alpha - 3.0),2)/pow(P+2.0+alpha,2);
		//K = (N -1)*pow(double(N-3),2)/pow(P,2);
	}
	else
	{
		//cout<<N+alpha-3<<"  "<<pow((N + alpha - 2),2)<<"  "<<pow(P+alpha,2)<<endl;
		K = (N + alpha -3)*pow((N + alpha - 2),2)/pow(P+2.0+alpha,2);
		//K = double((N+alpha-3)*(N-2)*(N-2))/pow(P+alpha,2);
	}

	//cout<<"K is "<<K<<endl;
	
	return K;
}

double Surfqsar::calculateCP()
{
	double cp = 0.0;
	
	cp = -264 + 86.1*log(double(NumEOUnits)) + 8.02*KierShapeThird() + 1284*ABInfoContentZero() - 14.26*StructInfoContentOne();

	return cp;
}

double Surfqsar::ABInfoContentZero()
{
	double ABIC = 0.0;
	int Hydrogens = 0;
	int Carbons=0;
	int bonds = 0;
	for (int i=0;i<mol->getsize();i++)
	{
		if (mol->getatom(i)->get_atomtype_name()=="C")
		{
			Hydrogens+=mol->getHydrogens(i);
			Carbons+=1;
		}
		bonds+=mol->getNN(i);
	}
	bonds=bonds/2;
	bonds--;//removing one C-O bond
	bonds+=Hydrogens;//accounting for C-H bond
	double carbonFraction = double(Carbons)/double(Hydrogens+Carbons);
	double hydrogenFraction = 1.0 - carbonFraction;
	//cout<<Carbons<<"  "<<Hydrogens<<endl;
	double IC = -carbonFraction*(log10(carbonFraction)/log10(2.0)) -hydrogenFraction*(log10(hydrogenFraction)/log10(2.0));
	//cout<<bonds<<"  "<<IC<<endl;
	ABIC = IC/(log10(double(bonds))/(log10(2.0)));
	//cout<<"ABIC is"<<ABIC<<endl;
	return ABIC;
}

double Surfqsar::StructInfoContentOne()
{
	map<int,int> ValueFreqMapForC;
	int Hydrogens = 0;
	for (int i=0;i<mol->getsize();i++)
	{
		int AtomValue = 1;
		if (mol->getatom(i)->get_atomtype_name()=="C")
		{
			Hydrogens+=mol->getHydrogens(i);
			for (int j=0;j<mol->getNN(i);j++)
			{
				int Neighbor = mol->get_adjacency(i,j);
				if (mol->getatom(Neighbor)->get_atomtype_name()=="C")
					AtomValue=AtomValue*3;
				if (mol->getatom(Neighbor)->get_atomtype_name()=="O")
					AtomValue=AtomValue*5;
				AtomValue = AtomValue * int(pow(2,double(mol->getHydrogens(i))));
			}

			if (ValueFreqMapForC.count(AtomValue)!=0)
				ValueFreqMapForC[AtomValue]++;
			else
				ValueFreqMapForC[AtomValue]=1;
		}

	}
	double ICOne =0.0;

	map<int,int>::iterator it;
	
	for (it =ValueFreqMapForC.begin();it!=ValueFreqMapForC.end();it++)
	{
		
		//cout<<it->second<<"  "<<Hydrogens+mol->getsize()-1<<endl;
		double fraction =(double(it->second)/double(Hydrogens+mol->getsize()-1));
		
		ICOne = ICOne - fraction*(log10(fraction)/log10(2.0));
	}
	double HydFraction = double(Hydrogens)/double(Hydrogens+ mol->getsize()-1);
	ICOne = ICOne - HydFraction*(log10(HydFraction)/log10(2.0));
	double SIC = (Hydrogens+ mol->getsize()-1)*ICOne/(log10(double(Hydrogens+mol->getsize()-1))/log10(2.0));
	//cout<<"SIC is "<<SIC<<endl;

	return SIC;
}

double Surfqsar::calculateHLB()
{
   return (3.674 + 10.215*pow(double(NumEOUnits+1), double(0.5)) + 0.362*KHZeroCalc()-1.119*double(NumEOUnits+1) - 4.825*pow(KHZeroCalc(),double(0.5)));
}



				
					





	






		


					






