     
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

#include "rng.h"
#include "Classheader.h"
#include "additionalfunc.h"
#include "stringreg.h"
#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <idas/idas_spgmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#define Ith(v,i)    NV_Ith_S(v,i)  

vector<double> getrateVector(N_Vector yy, void *user_data )
{
	realtype *yval;
	yval = NV_DATA_S(yy);
	UserData* data;
	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);
	vector<double> rateVector;
	vector<double> yvaltemp;

	ofstream kinfile("kineticsRecheck.txt");

	

	/*for (int i=0;i<NumOfEquations;i++)
	{
		if (yval[i]>=0.0)
			yvaltemp.push_back(yval[i]);
		else
			yvaltemp.push_back(0.0);
	}*/

	for (int i=0;i<data->RateInfo->size();i++)
	{
		//realtype rateValue = data->KineticsValueInfo->at(i);
		realtype rateValue = data->kinetics->at(i).getKineticsValue(*data->Temp);
		//cout<<"kinetics of "<<i<<" is "<<data->kinetics->at(i).getKineticsValue(*data->Temp)<<endl;
		kinfile<<rateValue<<endl;
		/*If the reaction kinetics is for a forward reaction, then multiply with a frequency factor, else not*/
		if (data->RevKinSet->count(i)==0)
			rateValue = rateValue*(data->Reactions->at(i).get_frequency());
		
		//cout<<"rate value is "<<rateValue<<endl;
		/*
		calculate the rate of each reaction
		1. for each reaction, use RateInfo vector to find out the species that go into forming the rate expression
		2. figure out the stoichiometric coeff of these species and appropriately calculate the product of species concentration raised to the stoichiometric index
		3. Note that the yy vector contains molar flow rate values and needs to be converted into concentration by diviing each value by the vol flow rate (the last entry of yy)
		*/
		set<int>::iterator it2;
			
		for (it2=data->RateInfo->at(i).begin();it2!=data->RateInfo->at(i).end();it2++)
		{
			int ReactantStoichCoeff =(*data->Reactions)[i].getReactantStoich((*data->SpeciesIndex)[*it2]);//coeff of the species corresponding to it2 in reactants side of reaction rxn
			//if (abs(ReactantStoichCoeff)!=1)cout<<"coeff is not +/-1 but "<<ReactantStoichCoeff<<endl;
			realtype ConcValue = yval[*it2]; //Concentration -> note that yval is flow rate for non heterogeneous surface intermediates! 
			if (data->surfaceSpecies->count(*it2)==0)//if not surface species
				ConcValue = ConcValue/(yval[NumOfEquations-1]);
			if (abs(ReactantStoichCoeff)!=1 )
				rateValue*=pow(ConcValue,abs(ReactantStoichCoeff));
			else rateValue*=ConcValue;

			//cout<<"rateValue now is "<<rateValue<<endl;
				
			//NOTE:yval[NumOfEquations-1] is the vol flow rate - so I am calculating the ratio of molar flow rate and vol flow rate to get the concentration
		}


		rateVector.push_back(rateValue);
	}

	

	/*for (int i =0;i<rateVector.size();i++)
	{
		cout<<rateVector[i]<<endl;
	}*/
	

	
	//kinfile.close();

	

	return rateVector;
}

vector<double> NumJacColumnCalc(N_Vector yy, void *user_data, int j)
{
	realtype *yval, *yval1, *yval2;
	yval = NV_DATA_S(yy);
	UserData* data;
	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);

	vector<double> JacColumnJ, JacColumnJ2;

	N_Vector yy1, yy2;

	yy1 = NULL;
	yy2 = NULL;

	yy1 = N_VNew_Serial(NumOfEquations);
	
	
	yy2 = N_VNew_Serial(NumOfEquations);
	
	yval1 = NV_DATA_S(yy1);
	yval2 = NV_DATA_S(yy2);

	for (int i = 0; i <NumOfEquations;i++)
	{
		yval1[i]=yval[i];
		yval2[i]=yval[i];
	}

	double deltachange = 1.0e-8;
	//if (deltachange<=1.0e-20)
		//deltachange = 1.0e-18;
	yval2[j]+=deltachange;


	JacColumnJ = getrateVector(yy2, data);
	JacColumnJ2 = getrateVector(yy1,data);

	for (int i =0;i<JacColumnJ.size();i++)
	{
		JacColumnJ[i] = (JacColumnJ[i]-JacColumnJ2[i])/deltachange;
	}


	N_VDestroy_Serial(yy1);
	N_VDestroy_Serial(yy2);

	return JacColumnJ;
}


vector<double> DrUponDyColumnCalc(N_Vector yy, void *user_data, int j, vector<double>& rateVector)
{
	realtype *yval;
	yval = NV_DATA_S(yy);
	UserData* data;
	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);

	vector<double> DrUponDyJ; //jth column of DrUponDy

	for (int i=0;i<data->RateInfo->size();i++)
	{
		if (j<NumOfEquations-1 && data->RateInfo->at(i).count(j)==0)
			DrUponDyJ.push_back(0.0);
		else
		{
			int totalIndexForNu  = 0;
			set<int>::iterator it2;

			for (it2=data->RateInfo->at(i).begin();it2!=data->RateInfo->at(i).end();it2++)
			{			
				int ReactantStoichCoeff =(*data->Reactions)[i].getReactantStoich((*data->SpeciesIndex)[*it2]);//coeff of the species corresponding to it2 in reactant side of the reaction rxn
				totalIndexForNu+=abs(ReactantStoichCoeff);
				if (*it2 ==j)
				{
					if (yval[i]==0.0)cout<<"oops zero value"<<endl;
					double drdy = RCONST(abs(ReactantStoichCoeff))*rateVector[i]/yval[i];
					
					
					DrUponDyJ.push_back(drdy);
					break;
				}
			}
			if (j==NumOfEquations-1)
			{
				DrUponDyJ.push_back(rateVector[i]*(-totalIndexForNu)/yval[NumOfEquations-1]);
			}
		}
	}
	return DrUponDyJ;
}		

vector<double> JacobianColumn(N_Vector yy, void *user_data, int j)
{
	realtype *yval;
	yval = NV_DATA_S(yy);
	UserData* data;
	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);

	vector<double> JacobianVector;

	for (int i=0;i<data->RateInfo->size();i++)
	{
		if (j<NumOfEquations-1 && data->RateInfo->at(i).count(j)==0)
			JacobianVector.push_back(0.0);
		else
		{
			//realtype jacobianValue = (data->KineticsValueInfo->at(i));
			realtype jacobianValue = data->kinetics->at(i).getKineticsValue(*data->Temp);
			if (data->RevKinSet->count(i)==0)
				jacobianValue = jacobianValue*(data->Reactions->at(i).get_frequency());
			
				set<int>::iterator it2;
		
			int ToatlIndexForNu = 0; //counts the number of conc terms in jacobianValue -> that gives the number of Nu (vol flow rate) exists in that term. This is needed for calculating jacobian entry wrt Nu.
			
			for (it2=data->RateInfo->at(i).begin();it2!=data->RateInfo->at(i).end();it2++)
			{
				int ReactantStoichCoeff =(*data->Reactions)[i].getReactantStoich((*data->SpeciesIndex)[*it2]);//coeff of the species corresponding to it2 in reactant side of the reaction rxn
				if (data->surfaceSpecies->count(*it2)==0)
					ToatlIndexForNu+=abs(ReactantStoichCoeff);
				realtype concValue;
				if (*it2!=j)//taking the product of all the terms that are not the variable j
				{
					concValue = yval[*it2];
					if (data->surfaceSpecies->count(*it2)==0)
						concValue = concValue/yval[NumOfEquations-1];
											
					if (abs(ReactantStoichCoeff)!=1 )
						jacobianValue*=pow(concValue,abs(ReactantStoichCoeff));
					else jacobianValue*=concValue;
					//NOTE:yval[NumOfEquations-1] is the vol flow rate - so I am calculating the ratio of molar flow rate and vol flow rate to get the concentration
				}
				//accounting for >1 stoich coefficient of variable j here
				else
				{
					//cout<<"yval[NumOfEquations] is "<<yval[NumOfEquations]<<endl;
					concValue = yval[*it2];
					if (data->surfaceSpecies->count(*it2)==0)
						concValue = concValue/yval[NumOfEquations-1];
					if (abs(ReactantStoichCoeff)>1)
						jacobianValue*=RCONST(abs(ReactantStoichCoeff))*pow(concValue,abs(ReactantStoichCoeff-1));
					if (data->surfaceSpecies->count(*it2)==0)
						jacobianValue = jacobianValue/yval[NumOfEquations-1]; //we need to ensure yval[NumOfEquations-1] is always raised to inverse power of TotalIndexForNu for flows

				}
						
			}
			if (j==NumOfEquations -1)//if jth variable is the last one -> Nu (vol flow rate)
			{
				jacobianValue*=(-ToatlIndexForNu)/yval[NumOfEquations-1];
				
			}
			JacobianVector.push_back(jacobianValue);
		}
	}

	return JacobianVector;
}	


vector<double> getNetRateVector(vector<double>& rateVector, void *user_data)
{
	vector<double> NetRate;
	UserData* data;

	data = (UserData*)user_data;

	int NumOfEquations = *(data->NumOfEquations);
	NetRate.resize(NumOfEquations-1,0.0);

	for (int i=0;i<NumOfEquations-2;i++)//note that stoichiometry is required for the first numOfEqns -1 only, the last one being hte equation for vol flow rate
	{
		
		multimap<int,pair<int,int> >::iterator stoich_it;
		pair<multimap<int,pair<int,int> >::iterator, multimap<int,pair<int,int> >::iterator > ret;
		
		ret = data->StoichInfo->equal_range(i);
			
						 
		for (stoich_it=ret.first;stoich_it!=ret.second;stoich_it++)
		{
			int rxn = stoich_it->second.first;
			realtype addTerm = rateVector.at(rxn);
			int ithCoeff=(*data->Reactions)[rxn].NetOccurence((*data->SpeciesIndex)[i]);
			//if (*(*data->SpeciesIndex)[i]=="O(C)C" && abs(addTerm)>0.00000000000001)cout<<i<<"  "<<*(*data->SpeciesIndex)[i]<<" "<<rxn<<"  "<<addTerm<<" "<<ithCoeff<<endl;
				
			NetRate[i]+=ithCoeff*addTerm;//because it is dy/dt - S.(nu) where S is stoichiometric matrix and nu is the rate vector! 
		
		}
		
	}

	return NetRate;
}

double getSiteBalanceDynamicResidual(N_Vector yp, int i, void *user_data)
{
	realtype *ypval;
	ypval = NV_DATA_S(yp); 

	UserData* data;

	data = (UserData*)user_data;
	double SiteBalanceDyn = 0.0;

	SiteBalanceDyn = ypval[i];

	multimap<string, pair<int, int> >::iterator site_it;
	pair<multimap<string, pair<int, int> >::iterator, multimap<string, pair<int, int> >::iterator > pair_it;

	multimap<int,string>::iterator SiteBalanceCompAtoms;

	pair<multimap<int,string>::iterator, multimap<int,string>::iterator> SiteBalancePairIt;

	SiteBalancePairIt = data->SiteBalanceInfo->equal_range(i);

	for (SiteBalanceCompAtoms = SiteBalancePairIt.first; SiteBalanceCompAtoms!= SiteBalancePairIt.second; SiteBalanceCompAtoms++)
	{
		
		pair_it = data->SiteOccupants->equal_range(SiteBalanceCompAtoms->second);

		for (site_it=pair_it.first; site_it!=pair_it.second;site_it++)
		{
			if (i!=site_it->second.first)//add only if not site itself!
				SiteBalanceDyn+=ypval[site_it->second.first]*site_it->second.second;
				

		}
	}
	
	return SiteBalanceDyn;
}

double getSiteBalanceResidual(N_Vector yy, int i,void *user_data)
{
	realtype *yval;
	yval = NV_DATA_S(yy); 

	UserData* data;

	data = (UserData*)user_data;
	double SiteBalance = 0.0;

	SiteBalance = yval[i];

	multimap<string, pair<int, int> >::iterator site_it;
	pair<multimap<string, pair<int, int> >::iterator, multimap<string, pair<int, int> >::iterator > pair_it;

	multimap<int,string>::iterator SiteBalanceCompAtoms;

	pair<multimap<int,string>::iterator, multimap<int,string>::iterator> SiteBalancePairIt;

	SiteBalancePairIt = data->SiteBalanceInfo->equal_range(i);

	for (SiteBalanceCompAtoms = SiteBalancePairIt.first; SiteBalanceCompAtoms!= SiteBalancePairIt.second; SiteBalanceCompAtoms++)
	{
		
		pair_it = data->SiteOccupants->equal_range(SiteBalanceCompAtoms->second);

		for (site_it=pair_it.first; site_it!=pair_it.second;site_it++)
		{
			if (i!=site_it->second.first)//add only if not site itself!
				SiteBalance+=yval[site_it->second.first]*site_it->second.second;
				

		}
	}
	
	SiteBalance-=(*data->InitialMolecFlows)[i].first;//subtracting the initial site density!
	
	return SiteBalance;
}


/*
int DynamicCSTRresidual(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data)
{
	realtype *yval, *ypval, *rval;
	yval = NV_DATA_S(yy); 
	ypval = NV_DATA_S(yp); 
	rval = NV_DATA_S(resval);
	UserData* data;

	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);
	rval[NumOfEquations-1]=RCONST(0.0);
	double Volume = *(data->Volume);
	double VolFlow = yval[NumOfEquations-1];

	double AbsTol = data->AbsTolerance;
	if ((VolFlow < RCONST(-AbsTol*1e5)) || (yval[NumOfEquations-2] < RCONST(-AbsTol*1e5)) )
		return (1);

	

	for (int i =0;i<NumOfEquations;i++)
		rval[i]=0.0;
	

	vector<double> rateVector = getrateVector(yy, data);
	vector<double> NetSpeciesRate = getNetRateVector(rateVector,data);


	for (int i =0;i<NumOfEquations-2;i++)
	{
		if (data->AllSites->count(i)>0)//if site, do only site balance
		{
			rval[i]=getSiteBalanceDynamicResidual(yp,i,data);
		}
		else if (data->surfaceSpecies->count(i)==0)
		{

			rval[i]+=ypval[i]*(Volume/VolFlow)-yval[i]*Volume/(VolFlow*VolFlow)*ypval[NumOfEquations-1] - NetSpeciesRate[i]*Volume;
				
			rval[i]+=yval[i];

			
			
			if(data->InitialMolecFlows->count(i)!=0)
			{
				rval[i]-=(*data->InitialMolecFlows)[i].first;
				rval[NumOfEquations-2]-=(*data->InitialMolecFlows)[i].first;
			}
			rval[NumOfEquations-2]+=yval[i];
			rval[NumOfEquations-2]-=NetSpeciesRate[i]*Volume;
			rval[NumOfEquations-1]+=yval[i];

		}
		else
			rval[i]=ypval[i]-NetSpeciesRate[i];
	}

	double totalConc = yval[NumOfEquations-2]/(data->Rg*(*data->Temp));	

	rval[NumOfEquations-2]+= ypval[NumOfEquations-2]*Volume/data->Rg*(*data->Temp);

	rval[NumOfEquations-1]=-rval[NumOfEquations-1];
	rval[NumOfEquations-1]+=yval[NumOfEquations-1]*totalConc;

	ofstream outputfile;
	outputfile.open("residualCSTR.txt");
	
	for (int i=0;i<NumOfEquations-2;i++)
	{
		outputfile<<*((*data->SpeciesIndex)[i])<<" "<<rval[i]<<"  "<<yval[i]<<"   "<<ypval[i]<<endl;
		
	}
	outputfile<<"Pressure "<<rval[NumOfEquations-2]<<" "<<yval[NumOfEquations-2]<<"  "<<ypval[NumOfEquations-2]<<endl;
	outputfile<<"Nu "<<rval[NumOfEquations-1]<<" "<<yval[NumOfEquations-1]<<"  "<<ypval[NumOfEquations-1]<<endl;

	outputfile<<"At "<<tres<<endl;
	outputfile.close();

	cout<<"returning"<<endl;

	return (0);
}
*/


int DynamicCSTRresidual(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data)
{
	realtype *yval, *ypval, *rval;
	yval = NV_DATA_S(yy); 
	ypval = NV_DATA_S(yp); 
	rval = NV_DATA_S(resval);
	UserData* data;

	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);
	rval[NumOfEquations-1]=RCONST(0.0);
	double Volume = *(data->Volume);
	double VolFlow = yval[NumOfEquations-1];

	double AbsTol = data->AbsTolerance;
	//if ((VolFlow < RCONST(-AbsTol*1e5)) || (yval[NumOfEquations-2] < RCONST(-AbsTol*1e5)) )
	//	return (1);

	

	N_Vector FlowRates;//for calculating rate
	realtype *flow;

	FlowRates = NULL;
	flow = NULL;
	

	FlowRates = N_VNew_Serial(NumOfEquations);
	
	flow = NV_DATA_S(FlowRates);

	for (int i =0;i<NumOfEquations-1;i++)
	{
		//if site use only conc, else use flow rates! 
		if (data->surfaceSpecies->count(i)==0) flow[i]=yval[i]*yval[NumOfEquations-1];
		else flow[i]=yval[i];
		rval[i]=RCONST(0.0);
	}
	flow[NumOfEquations-1]=yval[NumOfEquations-1];
	
	rval[NumOfEquations-1]=0.0;
	

	vector<double> rateVector = getrateVector(FlowRates, data);
	vector<double> NetSpeciesRate = getNetRateVector(rateVector,data);


	for (int i =0;i<NumOfEquations-1;i++)
	{
		if (data->AllSites->count(i)>0)//if site, do only site balance
		{
			rval[i]=getSiteBalanceDynamicResidual(yp,i,data);
		}
		else if (data->surfaceSpecies->count(i)==0)
		{
			rval[i]+= ypval[i]*Volume - Volume*NetSpeciesRate[i];
			rval[i]+= yval[i]*yval[NumOfEquations-1];

			if(data->InitialMolecFlows->count(i)!=0)
				rval[i]-=(*data->InitialMolecFlows)[i].first;

			//rval[NumOfEquations-2]-=(data->Rg*(*data->Temp))*ypval[i];
			rval[NumOfEquations-1]+=yval[i];

		}
		else
			rval[i]=ypval[i]-NetSpeciesRate[i];
	}

	rval[NumOfEquations-1]-=yval[NumOfEquations-2]/(data->Rg*(*data->Temp));	


	ofstream outputfile;
	outputfile.open("residualCSTR.txt");
	
	for (int i=0;i<NumOfEquations-1;i++)
	{
		outputfile<<*((*data->SpeciesIndex)[i])<<" "<<rval[i]<<"  "<<yval[i]<<"   "<<ypval[i]<<endl;
		
	}
	outputfile<<"Nu "<<rval[NumOfEquations-1]<<" "<<yval[NumOfEquations-1]<<"  "<<ypval[NumOfEquations-1]<<endl;

	outputfile<<"At "<<tres<<endl;
	outputfile.close();

	//cout<<"returning"<<endl;

	return (0);
}






/*
int DynamicCSTRresidual(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data)
{

	realtype *yval, *ypval, *rval;
	yval = NV_DATA_S(yy); 
	ypval = NV_DATA_S(yp); 
	rval = NV_DATA_S(resval);
	UserData* data;

	data = (UserData*)user_data;
	int NumOfEquations = *(data->NumOfEquations);
	rval[NumOfEquations-1]=RCONST(0.0);
	double Volume = *(data->Volume);
	double VolFlow = yval[NumOfEquations-1];
	double RTOverP = (data->Rg*(*data->Temp))/(*data->Pressure);

	double AbsTol = data->AbsTolerance;
	if (VolFlow < RCONST(-AbsTol*1000))
		return (1);
	
	cout<<"entered CSTR residuals"<<endl;
	

	
	N_Vector FlowRates;//for calculating rate
	realtype *flow;

	FlowRates = NULL;
	flow = NULL;
	
	FlowRates = N_VNew_Serial(NumOfEquations);
	
	flow = NV_DATA_S(FlowRates);

	for (int i =0;i<NumOfEquations-1;i++)
	{
		//if site use only conc, else use flow rates! 
		if (data->surfaceSpecies->count(i)==0) flow[i]=yval[i]*yval[NumOfEquations-1];
		else flow[i]=yval[i];
		rval[i]=RCONST(0.0);
	}
	flow[NumOfEquations-1]=yval[NumOfEquations-1];
	rval[NumOfEquations-1]=0.0;

	vector<double> rateVector = getrateVector(FlowRates, data);
	vector<double> NetSpeciesRate = getNetRateVector(rateVector,data);

	//for (int i =0;i<NumOfEquations;i++)
	//	rval[i]=0.0;
	
	cout<<"first gas phase product is "<<*(*data->SpeciesIndex)[*data->FirstFoundProductFlowSpecies]<<endl;
	
	
	
	

	//rval[*data->FirstFoundProductFlowSpecies]+=ypval[*data->FirstFoundProductFlowSpecies]-ypval[*data->FirstFoundProductFlowSpecies]/RTOverP;
	
	rval[NumOfEquations-1]+=NetSpeciesRate[*data->FirstFoundProductFlowSpecies]*Volume - yval[NumOfEquations-1]/RTOverP;
	cout<<"first it is "<<rval[NumOfEquations-1]<<"  "<<yval[NumOfEquations-1]/RTOverP<<endl;
	rval[*data->FirstFoundProductFlowSpecies]+=yval[*data->FirstFoundProductFlowSpecies];

	

	for (int i = 0;i<NumOfEquations -1; i++)
	{
		if (data->AllSites->count(i)>0)//if site, do only site balance
		{
			//rval[i]=getSiteBalanceResidual(yy,i,data);
			rval[i]=getSiteBalanceDynamicResidual(yp,i,data);
		}
		else
		{
			if (data->surfaceSpecies->count(i)==0)
			{
						

				//rval[i]+=ypval[i]*(Volume/VolFlow)-yval[i]*Volume/(VolFlow*VolFlow)*ypval[NumOfEquations-1] - NetSpeciesRate[i]*Volume;
				
				//if (*(*data->SpeciesIndex)[i]=="C=C")
				//	cout<<ypval[i]*(Volume/VolFlow)<<"  "<<yval[i]*Volume/(VolFlow*VolFlow)*ypval[NumOfEquations-1]<<"  "<<NetSpeciesRate[i]*Volume<<rval[i]<<"  "<<yval[i]<<endl;
				//rval[i]+=yval[i];

				//if (*(*data->SpeciesIndex)[i]=="C=C")
				//	cout<<rval[i]<<endl;


				rval[i]+= ypval[i]*Volume - Volume*NetSpeciesRate[i];
				rval[i]+= yval[i]*yval[NumOfEquations-1];
				if(data->InitialMolecFlows->count(i)!=0)
				{
					rval[i]-=(*data->InitialMolecFlows)[i].first;
					//if (*(*data->SpeciesIndex)[i]=="C=C")
					//	cout<<rval[i]<<endl;
					rval[NumOfEquations-1]+=(*data->InitialMolecFlows)[i].first;
					if (*(*data->SpeciesIndex)[i]=="C=C")
						cout<<rval[NumOfEquations-1]<<endl;

				}
				rval[NumOfEquations-1]+= Volume*NetSpeciesRate[i];
				cout<<rval[NumOfEquations-1]<<endl;
				//rval[*data->FirstFoundProductFlowSpecies]+=ypval[i];
			
				//if (*(*data->SpeciesIndex)[i]=="O(C)C")cout<<i<<"  "<<*(*data->SpeciesIndex)[i]<<"  "<<rval[i]<<endl;
		
				//rval[NumOfEquations-1]+=ypval[i];
				rval[*data->FirstFoundProductFlowSpecies]+=yval[i];
				
			}
			else rval[i]=ypval[i]-NetSpeciesRate[i];
		}
	}

	rval[*data->FirstFoundProductFlowSpecies]-=(*data->Pressure)/(data->Rg*(*data->Temp));

	//rval[NumOfEquations-1]=(-rval[NumOfEquations-1]*data->Rg*(*data->Temp))/(*data->Pressure); //Gas constant is 8.314 Pa m^3/K.mol or 0.0821 l atm/K.mol
	//rval[NumOfEquations-1]+=yval[NumOfEquations-1];
 

	ofstream outputfile;
	outputfile.open("residualCSTR.txt");
	
	for (int i=0;i<NumOfEquations-1;i++)
	{
		outputfile<<*((*data->SpeciesIndex)[i])<<" "<<rval[i]<<"  "<<yval[i]<<endl;
		
	}
	outputfile<<"Nu "<<rval[NumOfEquations-1]<<" "<<yval[NumOfEquations-1]<<endl;

	outputfile.close();

	cout<<"returning"<<endl;

	return (0);

}

*/

/*
int ICforDerivativesDynCSTR(vector<double>& ypval, vector<double> Rval, void *user_data)
{
	UserData* data;
	data = (UserData*)user_data;
	double ypvalSumGasPhase = 1000.0;
	
	
	while (abs(ypval[NumOfEquations-1]-ypvalSumGasPhase)< 1.0e-3*ypval[NumOfEquations-1])
	{
		for (int i = 0;i<NumOfEquations-1;i++)
		{	
			if (data->surfaceSpecies->count(i)==0)
			{
				ypval[i] = -Rval[i];
				ypvalSumGasPhase+=ypval[i];
			}
		}
		ypval[NumOfEquations-1] = ypvalSumGasPhase;
	}
}*/

int Residual(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data)
{
	realtype *yval, *ypval, *rval;
	yval = NV_DATA_S(yy); 
	ypval = NV_DATA_S(yp); 
	rval = NV_DATA_S(resval);
	UserData* data;

	data = (UserData*)user_data;
	

	int NumOfEquations = *(data->NumOfEquations);
	rval[NumOfEquations-1]=RCONST(0.0);
	
	double AbsTol = data->AbsTolerance;
	vector<realtype> yvaltemp;

	

	/*for (int i=0;i<NumOfEquations;i++)
	{

		if (yval[i]<RCONST(-AbsTol*1.0e5))
		{
			return (1);
		}
	}*/
	

	vector<double> rateVector = getrateVector(yy,data);
	
	//for (int i =0;i<rateVector.size();i++)
		//cout<<i<<"  "<<rateVector[i]<<endl;
		
	
	for (int i=0;i<NumOfEquations-1;i++)//note that stoichiometry is required for the first numOfEqns -2 only, the last one being hte equation for vol flow rate
	{
		
		if (data->AllSites->count(i)>0)//if site, do only site balance
		{
			rval[i] = yval[i];

			multimap<string, pair<int, int> >::iterator site_it;
			pair<multimap<string, pair<int, int> >::iterator, multimap<string, pair<int, int> >::iterator > pair_it;

			multimap<int,string>::iterator SiteBalanceCompAtoms;

			pair<multimap<int,string>::iterator, multimap<int,string>::iterator> SiteBalancePairIt;

			SiteBalancePairIt = data->SiteBalanceInfo->equal_range(i);

			for (SiteBalanceCompAtoms = SiteBalancePairIt.first; SiteBalanceCompAtoms!= SiteBalancePairIt.second; SiteBalanceCompAtoms++)
			{
				
				pair_it = data->SiteOccupants->equal_range(SiteBalanceCompAtoms->second);

				for (site_it=pair_it.first; site_it!=pair_it.second;site_it++)
				{
					if (i!=site_it->second.first)//add only if not site itself!
						rval[i]+=yval[site_it->second.first]*site_it->second.second;
						

				}
			}
			
			rval[i]-=(*data->InitialMolecFlows)[i].first;//subtracting the initial site density!
			
		}
		/*else if (*data->FirstFoundProductFlowSpecies==i)//if first found product flow species
		{
			//cout<<"this is the flow species"<<*data->FirstFlowSpecies<<endl;
			//sum over all flow rates times their mol weight and equate it to inlet mass flow rate
			rval[i]=0.0;
			
			//cout<<"species "<<i<<*(*data->SpeciesIndex)[i]<<" is the first product species found"<<endl;
			for (int j = 0;j<NumOfEquations-1;j++)
			{
				if (data->surfaceSpecies->count(j)==0)//not a site
					rval[i]+=ypval[j] * (*data->SpeciesMolWeights)[j];
				//if ((*data->SpeciesMolWeights)[j]>100)cout<<(*data->SpeciesMolWeights)[j]<<endl;
			}
			//rval[i]-=*(data->InletMassFlowRate);
			rval[NumOfEquations-1]+=yval[i];//adding up the flow rates! 

			
		}*/
		else //if not site or first product flow species 
		{
			if (data->surfaceSpecies->count(i)==0)//derivative added only in the case of flows (non surface)
				rval[i]=ypval[i];//setting the yp value (ith derivative only present in the ith residue)
			else rval[i]=RCONST(0.0);
		
			multimap<int,pair<int,int> >::iterator stoich_it;
			pair<multimap<int,pair<int,int> >::iterator, multimap<int,pair<int,int> >::iterator > ret;
			
			ret = data->StoichInfo->equal_range(i);
			
			//if (*(*data->SpeciesIndex)[i]=="N#N")cout<<"for nitrogen "<<rval[i]<<"  "<<ypval[i]<<endl;
						 
			for (stoich_it=ret.first;stoich_it!=ret.second;stoich_it++)
			{
				int rxn = stoich_it->second.first;
				realtype addTerm = rateVector.at(rxn);
				int ithCoeff=(*data->Reactions)[rxn].NetOccurence((*data->SpeciesIndex)[i]);
				//if (*(*data->SpeciesIndex)[i]=="C" && abs(addTerm)>0.00000000000001)cout<<i<<"  "<<*(*data->SpeciesIndex)[i]<<" "<<rxn<<"  "<<addTerm<<" "<<ithCoeff<<"  "<<rval[i]<<endl;
				rval[i]-=ithCoeff*addTerm;//because it is dy/dt - S.(nu) where S is stoichiometric matrix and nu is the rate vector! 
				//if (*(*data->SpeciesIndex)[i]=="O(C)C" && abs(addTerm)>0.00000000000001)cout<<rval[i]<<endl;
			}
			if (data->surfaceSpecies->count(i)==0)
				rval[NumOfEquations-1]+=yval[i];//adding up the molar flow rates first, remainining part of the equation done outside the loop
		}
	}
	
	//cout<<rval[NumOfEquations-1]<<endl;
	rval[NumOfEquations-1]=(-rval[NumOfEquations-1]*data->Rg*(*data->Temp))/(*data->Pressure); //Gas constant is 8.314 Pa m^3/K.mol or 0.0821 l atm/K.mol
	rval[NumOfEquations-1]+=yval[NumOfEquations-1];//doing the volumetric flow rate balance here as mu (m^3/unit.time) = Ftotal*R*T/P
 
	//cout<<"Residual returning"<<endl;

	/*ofstream outputfile;
	outputfile.open("residual.txt");
	
	for (int i=0;i<NumOfEquations-1;i++)
	{
		outputfile<<*((*data->SpeciesIndex)[i])<<" "<<rval[i]<<"  "<<yval[i]<<"  "<<ypval[i]<<endl;
		//if (*(*data->SpeciesIndex)[i]=="O(C)C")
			//cout<<rval[i]<<"  "<<yval[i]<<"  "<<ypval[i]<<endl;
		
	}
	outputfile<<"Nu "<<rval[NumOfEquations-1]<<" "<<yval[NumOfEquations-1]<<endl;

	outputfile<<"At "<<tres<<endl;

	


	outputfile.close();*/
	
	return (0);

}


int Jacobian(long int Neq, realtype tt,  realtype cj, 
		   N_Vector yy, N_Vector yp, N_Vector resvec,
		   DlsMat JJ, void *user_data,
		   N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{

	/*constructing the jacobian is similar to that of the residual - only difference being
		that the jacobian elements of the reaction rate vector is calculated on-the-fly from RateInfo.*/

	realtype *yval;
	yval = NV_DATA_S(yy);
	UserData* data;
	data = (UserData*)user_data;

	int NumOfEquations = *data->NumOfEquations;
	
	//cout<<"called Jacobian"<<endl;
	//vector<double> rateVector = getrateVector(yy, data);

	map<int, vector<double> > jacobianVectors;

	for (int j=0;j<NumOfEquations;j++)
	{
		jacobianVectors[j]=JacobianColumn(yy, data, j);
	}

	
	for (int i=0;i<NumOfEquations-1;i++)
	{
		if (data->surfaceSpecies->count(i)==0)
			(DENSE_COL(JJ,i))[i] = cj;
		//DENSE_ELEM(JJ,i,i) = cj;

		if (data->AllSites->count(i)>0)
		{
			string site = *(*data->SpeciesIndex)[i]; //string corresponding to the site
			(DENSE_COL(JJ,i))[i] = RCONST(1.0);
			multimap<string, pair<int, int> >::iterator site_it;
			pair<multimap<string, pair<int, int> >::iterator, multimap<string, pair<int, int> >::iterator > pair_it;

			multimap<int,string>::iterator SiteBalanceCompAtoms;

			pair<multimap<int,string>::iterator, multimap<int,string>::iterator> SiteBalancePairIt;

			SiteBalancePairIt = data->SiteBalanceInfo->equal_range(i);

			for (SiteBalanceCompAtoms = SiteBalancePairIt.first; SiteBalanceCompAtoms != SiteBalancePairIt.second; SiteBalanceCompAtoms++)
			{
				
				pair_it = data->SiteOccupants->equal_range(SiteBalanceCompAtoms->second);

				for (site_it=pair_it.first; site_it!=pair_it.second;site_it++)
				(DENSE_COL(JJ,site_it->second.first))[i] = RCONST(site_it->second.second);
			}
			
		}
		else
		{

			for (int j=0;j<NumOfEquations;j++)
			{
				multimap<int,pair<int,int> >::iterator it;
				pair<multimap<int,pair<int,int> >::iterator, multimap<int,pair<int,int> >::iterator > ret;
			
				ret = data->StoichInfo->equal_range(i);

				vector<double> jacobianVectorForJ=jacobianVectors[j];

				
				//vector<double> jacobianVectorForJ=NumJacColumnCalc(yy, data, j);

				//vector<double> analytical = JacobianColumn(yy, data, j);
				
				/*for (int k = 0;k<analytical.size();k++)
				{
					double diff = (analytical[k]-jacobianVectorForJ[k]);

					if (diff > 1.0e-5 || diff < -1.0e-5)
					{
						cout<<"huge difference "<<diff<<"  "<<j<<"  "<<k<<endl;
						cout<<analytical[k]<<"  "<<jacobianVectorForJ[k]<<endl;
						//THIS is WRONGcout<<*(*data->SpeciesIndex)[j]<<"  "<<data->Reactions->at(k).reactionstring()<<endl;
					}
				}*/
				
				//vector<double> jacobianVectorForJ = DrUponDyColumnCalc(yy, data,j,rateVector);
				
				for (it=ret.first;it!=ret.second;it++)
				{

					//check if the rate vector of that reaction even has the jth variable
					int rxn = it->second.first;
					realtype addTerm= jacobianVectorForJ.at(rxn);
					int ithCoeff=(*data->Reactions)[rxn].NetOccurence((*data->SpeciesIndex)[i]);
					
					//if (ithCoeff*addTerm > 1.0e5)
					//	cout<<"adding a large value "<<ithCoeff*addTerm<<endl;
					
					(DENSE_COL(JJ,j))[i]-=ithCoeff*addTerm;//note that Jacobian J = I.cj - S(dnu/dy)
					
					//DENSE_ELEM(JJ,i,j) -=ithCoeff*addTerm;
					 
					
				}
				
				
			}
		}
	}

	//filling in the last row of the jacobian ->entries for hte algebraic equaitons of Nu! 
	
	double RTOverP = -(data->Rg*(*data->Temp))/(*data->Pressure);
	for (int j=0;j<NumOfEquations-1;j++)
	{
		if (data->surfaceSpecies->count(j)==0)
			(DENSE_COL(JJ,j))[NumOfEquations-1]=RTOverP;
		else (DENSE_COL(JJ,j))[NumOfEquations-1]= RCONST(0.0);
	}

		//DENSE_ELEM(JJ,NumOfEquations-1,j) = RTOverP;

	(DENSE_COL(JJ,NumOfEquations-1))[NumOfEquations-1] = RCONST(1.0);
	//DENSE_ELEM(JJ,NumOfEquations-1,NumOfEquations-1) = RCONST(1.0);
	
	/*ofstream outputfile;
	outputfile.open("jacobian.txt");
	int ZerosCounter = 0;
	for (int i=0;i<NumOfEquations;i++)
	{
		if (i<NumOfEquations -1)
			outputfile<<*((*data->SpeciesIndex)[i])<<" ";
		else outputfile<<"nu  ";
		outputfile<<"--> "<<yval[i]<<" ";
		for (int j=0;j<NumOfEquations;j++)
		{
			outputfile<<" "<<DENSE_COL(JJ,j)[i];
			if (DENSE_COL(JJ,j)[i]==0.0) ZerosCounter++;

		}
		outputfile<<endl;
	}
	outputfile<<cj<<endl;

	outputfile<<"% zeros "<<float(ZerosCounter)/float(NumOfEquations*NumOfEquations)<<endl;
	cout<<"called jacobian with "<<cj<<endl;
			

	
	outputfile.close();*/
	return 0;

		

}

int resS(int Ns, realtype t, 
                N_Vector yy, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

	UserData* data;
	data = (UserData*)user_data;

	

	

	vector< vector<double> > jacobianTransp;

	int NumOfEquations = *data->NumOfEquations;
	//cout<<"Number of equations is "<<NumOfEquations<<endl;

	DlsMat Jac;
	N_Vector dummy1, dummy2, dummy3, dummyresvec;

	dummyresvec = NULL;

	Jac = NULL;
	Jac = NewDenseMat(NumOfEquations, NumOfEquations);

	SetToZero(Jac);

	/*for (int i = 0;i<NumOfEquations;i++)
	{
		for (int j = 0;j<NumOfEquations;j++)
			DENSE_COL(Jacobian,j)[i]=RCONST(0.0);
	}*/

	Jacobian(NumOfEquations, t,  RCONST(0.0), yy, yp, dummyresvec,
         Jac, data, tmp1,tmp2,tmp3);

	
	//for (int i = 0;i<NumOfEquations;i++)
	//	jacobianTransp.push_back(JacobianColumn(yy, data, i));

	vector<double> rateVector = getrateVector(yy,data);
	
	for (int j = 0; j< Ns; j++)
	{
		for (int i = 0;i<NumOfEquations;i++)
		{
			Ith(resvalS[j],i)=0.0;
			if (data->surfaceSpecies->count(i)==0 && i!=NumOfEquations-1)
				Ith(resvalS[j],i)=Ith(ypS[j],i);
			for (int l = 0;l<NumOfEquations;l++)
				Ith(resvalS[j],i)+=DENSE_COL(Jac,l)[i]*Ith(yyS[j],l);
			
		}

		vector<double> rateVectorForJ;
		rateVectorForJ.resize(rateVector.size(),0.0);//it is nonzero only for j and and the reverse reactions of j -first set it to zero and then calcualte the nonzero ones
		set<int> reverseRxns;
		multimap<int,int>::iterator it;

		pair<multimap<int,int>::iterator, multimap<int,int>::iterator> pair_it;

		pair_it = data->ReverseRxns->equal_range(j);

		for (it = pair_it.first;it!=pair_it.second;it++)
			reverseRxns.insert(it->second);

		/*for (int i =0;i<rateVector.size();i++)
		{
			
			if(i==j)
				//rateVectorForJ[i]=rateVector[i]/data->KineticsValueInfo->at(i);
				rateVectorForJ[i]=rateVector[i]/data->kinetics->at(i).getKineticsValue()
			else if (reverseRxns.count(i)!=0)
			{
				
				if (data->RevKinSet->count(i)==1)//rxn j is the original reaction and i is in the reverse set
					rateVectorForJ[i] = rateVector[i]/data->KineticsValueInfo->at(i) * (*data->RevKinSet)[i];
				else if (data->RevKinSet->count(j)==1) //rxn j is in the reverse set
					rateVectorForJ[i]=rateVector[i]/(data->KineticsValueInfo->at(i) * (*data->RevKinSet)[j]);
				else cout<<"there appears to be a freak case of both reactions not being defined as a reverse in terms of kinetics"<<endl;
			}
			else rateVectorForJ[i] = 0.0;
			
			
		}*/

		//using stoich info add the term dfi/dkj

		for (int i=0;i<NumOfEquations-1;i++)
		{
			if (data->AllSites->count(i)==0)//not sites
			{
				multimap<int,pair<int,int> >::iterator stoich_it;
				pair<multimap<int,pair<int,int> >::iterator, multimap<int,pair<int,int> >::iterator > ret;
			
				ret = data->StoichInfo->equal_range(i);
				
				for (stoich_it=ret.first;stoich_it!=ret.second;stoich_it++)
				{
					int rxn = stoich_it->second.first;
					realtype addTerm = rateVectorForJ.at(rxn);
					int ithCoeff=stoich_it->second.second;
					Ith(resvalS[j],i)-=ithCoeff*addTerm;//because it is ds/dt + jacobian term - S.(d(rate)/dkj) where S is stoichiometric matrix 
					//if (i==0 && j==3) cout<<ithCoeff*addTerm<<"  "<<addTerm<<" "<<rxn<<endl;
				}
			}
		}

		
		//no need to set anything for the sensitivity parameter d(flowrate)/dkj because it is zero! 
				
	}

	ofstream residualfile("residualsSens.txt");

	for (int i = 0;i<NumOfEquations;i++)
	{
		if (i<NumOfEquations -1) residualfile<<*((*data->SpeciesIndex)[i]);
		else residualfile<<"nu";
		for (int j = 0; j<Ns;j++)
		{
			residualfile<<" "<<Ith(resvalS[j],i)<<" ("<<Ith(yyS[j],i)<<")";
		}
		residualfile<<endl;
		
	}

	residualfile.close();

	DestroyMat(Jac);
	//N_VDestroy_Serial(dummy1);
//	N_VDestroy_Serial(dummy2);
//	N_VDestroy_Serial(dummy3);
//	N_VDestroy_Serial(dummyresvec);
	return 0;
}



int PreCondSetup(realtype tt, 
               N_Vector yy, N_Vector yp, N_Vector rr, 
               realtype cj, void *user_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

{
	UserData* data;
	data = (UserData*)user_data;

	realtype *yval;
	yval = NV_DATA_S(yy);

	
	

	int NumOfEquations = *data->NumOfEquations;
	//cout<<NumOfEquations<<endl;

	vector<double> preCond;
	preCond.resize(NumOfEquations);
	
	
	map<int, vector<double> > jacobianVectors;

	for (int j=0;j<NumOfEquations;j++)
	{
		jacobianVectors[j]=JacobianColumn(yy, data, j);
		preCond[j] = 0.0;
	}

	
	for (int i=0;i<NumOfEquations-1;i++)
	{
		if (data->surfaceSpecies->count(i)==0)
			preCond[i] = cj;

		if (data->AllSites->count(i)>0)
			preCond[i] = RCONST(1.0);
		else
		{
			//similar to jacobian, but note that I need jacobian wrt to i only for equation i. 
			//note that the inner for loop will work fine even though it has been directly taken from Jacobian. for reactions where jacobian wrt to i is nonzero, 
			//by definition they have i as reactants. 
			
			multimap<int,pair<int,int> >::iterator it;
			pair<multimap<int,pair<int,int> >::iterator, multimap<int,pair<int,int> >::iterator > ret;
		
			ret = data->StoichInfo->equal_range(i);

			vector<double> jacobianVectorForJ=jacobianVectors[i];				
			
			for (it=ret.first;it!=ret.second;it++)
			{
				int rxn = it->second.first;
				realtype addTerm= jacobianVectorForJ.at(rxn);
				int ithCoeff=(*data->Reactions)[rxn].NetOccurence((*data->SpeciesIndex)[i]);				
				preCond[i]-=ithCoeff*addTerm;			
			}
		}
	}
	
	preCond[NumOfEquations-1] = RCONST(1.0);

	//cout<<preCond.size()<<"  "<<data->PreconditionerInverse.size()<<endl;

	for (int i = 0;i<NumOfEquations-1;i++)
		data->PreconditionerInverse[i] = 1/(preCond[i]+ 1.0e-20);
				
	
		
	return (0);
}


int PreCondSolve(realtype tt, 
               N_Vector yy, N_Vector yp, N_Vector rr, 
               N_Vector rvec, N_Vector zvec, 
               realtype cj, realtype delta, void *user_data, 
               N_Vector tmp)
{
	UserData* data;
	data = (UserData*)user_data;
	realtype *zval;
	zval = NV_DATA_S(zvec);
	realtype *rval;
	rval = NV_DATA_S(rvec);

	int NumOfEquations = *data->NumOfEquations;

	for (int i =0;i<NumOfEquations-1;i++)
		zval[i] = rval[i]*data->PreconditionerInverse[i];


  
  return(0);
}

/*KineticModel::KineticModel(rxn_net_gen * net, char * file, vector<KineticParamPtr>& kinParam, double vol, double Pres, map<string,pair<double,int> >& Init, multimap<string, string> SiteBalance, ConstrPtr Outputs, double R, bool lumpedNetwork)
{
	
	Network = net;
	filename = file;
	//PreExpFactor =&A;
	//ActEnergy = &Ea;
	//N = &tempindex;
	kinetics =&kinParam;
	UseLumpedNetwork = lumpedNetwork;

	Temp = net->Temperature;
	Volume = vol;
	InletPressure = Pres;
	Rg = R;
	
	
	
	
	
}*/

KineticModel::KineticModel()
{
	OutputConstraintsSpecified = false;
	useAlternativeKinetics = false;
	doSensitivity = false;
	shouldCalculateRates = false;
}

void KineticModel::setNetwork(rxn_net_gen* net)
{
	Network = net;
	Temp = net->Temperature;
}

void KineticModel::printTofile(char* file)
{
	filename = file;
}

void KineticModel::setKineticFns(vector<KineticParamPtr>& kinParam)
{
	kinetics =&kinParam;
}

void KineticModel::setVolume(double vol)
{
	Volume = vol;
}

void KineticModel::setPressure(double Pres)
{
	InletPressure = Pres;
}

void KineticModel::setTemperature(double T)
{
	Temp = T;
}

void KineticModel::setInitialFlow(map<string,pair<double,int> >& Init)
{
	InitialFlowSpecified = Init;
}

void KineticModel::setSiteBalance(multimap<string,string> SiteBalance)
{
	SiteBalanceInput = SiteBalance;
}

void KineticModel::setOutputConstraint(ConstrPtr Outputs)
{
	RequiredOutputs = Outputs;
	OutputConstraintsSpecified = true;
}

void KineticModel::setRg(double R)
{
	Rg = R;
}

void KineticModel::setForLumpedNetwork(bool lumpedNetwork)
{
	UseLumpedNetwork = lumpedNetwork;
}

void KineticModel::setOutputSpecies(set<string> outputMols)
{
	specifiedOutputs = outputMols;
}

void KineticModel::doSensitivityAnalysis()
{
	doSensitivity = true;
}
void KineticModel::setAlternativeKinetics(vector<AltKineticParamPtr>& kinFns, vector<double> params)
{
	AlternativeKinetics = &kinFns;
	useAlternativeKinetics = true;
	AlternativeKinParams = params;
}

map<string, double> KineticModel::getModelOutputs()
{
	return OutputValues;
}

map<string, vector<double> > KineticModel::getModelSensitivityOutputs()
{
	return SensitivityOutputs;
}

void KineticModel::PrepareKineticModel()
{
	cout<<"Building the kinetic model"<<endl;
	FirstDcldFlowSpecies = -1;//initializing - value set in function GenerateStoichAndRateVectorInfo
	FirstFoundProductFlowSpecies = -1;
	InletMassFlowRate = 0.0; //value set in function GenerateStoichAndRateVectorInfo

	if (UseLumpedNetwork)
	{
		//NumOfEquations = net->AllLumpedMolecules.size()+2;
		//NumOfParameters = net->LumpedReconstructedNetwork.size();
		MoleculeSetPtr = &Network->AllLumpedMolecules;
		MolReactantMapPtr = &Network->MolLumpReactantMap;
		MolProductMapPtr = &Network->MolLumpProductMap;
		ReactionsPtr = &Network->LumpedReconstructedNetwork;
	}
	else
	{
		//NumOfEquations = net->AllMolecules.size()+1;
		//NumOfParameters = net->AllReactions.size();
		MoleculeSetPtr = &Network->AllMolecules;
		MolReactantMapPtr = &Network->MolReactantMap;
		MolProductMapPtr = &Network->MolProductMap;
		ReactionsPtr = &Network->AllReactions;
	}
	
	NumOfEquations = MoleculeSetPtr->size()+1;
	if (useAlternativeKinetics)
		NumOfParameters = AlternativeKinParams.size();
	else NumOfParameters = ReactionsPtr->size();

	KineticsValueInfo.resize(ReactionsPtr->size(),0.0);

	RateInfo.resize(ReactionsPtr->size());
	GenerateStoichAndRateVectorInfo(InitialFlowSpecified, SiteBalanceInput);
	if (FirstFoundProductFlowSpecies==-1)
	{
		cout<<"no product bulk phase flow species found! This is impossible! Stopping..."<<endl;
		ProceedToSolve = false;
	}
	else ProceedToSolve = true;

	if (ProceedToSolve && !SetKinetics())
	{
		cout<<"kinetics of reactions not set completely"<<endl;
		ProceedToSolve = false;
	}
}


int KineticModel::SolveModel()
{
	PrepareKineticModel();
	return Solve(0,false);
}





int KineticModel::Solve(int Type, bool isSolvingIC)
{
	/*variable names have been kept generic to match IDA: For our purpose, 
	yy = the y vector which in our case is the molar flow rates vector. For adsorbed intermediates, this will also include concentrations!
	yp = the time derivative of the y vector
	t = the time which in our case translates to the volume of the reactor! */

	/*NOTE: Type here refers to whether PFR (0) or dynamic CSTR(1) */
	/*NOTE: isSolvingIC is mainly for letting the function know if this is for solving the initial conditions (using a dynamic CSTR) or not*/

	/*Note that NumOfEquations - 1 is always the volumetric flow rate*/ 
	if (!ProceedToSolve)
	{
		cout<<"stopping from proceeding to solve the model"<<endl;
		return 1;
	}

	
	void *mem;
	N_Vector yy, yp, avtol, constraints, id;
	realtype rtol, *yval, *ypval, *atval, *idval, sensrTol, *sensAtval;
	realtype t0, tout1, tout, tret, tout2;
	int retval, flag;
	N_Vector *yS, *ypS, *sensavtol;
	
	double AbsTolerance = 1.0e-4;
	double AbsTolerance2 = 1.0e-4;
	double sensAbsTol = 1.0e-12;

	mem = NULL;
	yy = yp = avtol = id = NULL;
	yval = ypval = atval = idval = sensAtval = NULL;

	bool NeedToCheckIC = false;

	double EffVolume = 0.0;

	EffVolume = Volume;
	if (isSolvingIC) EffVolume = Volume/1000;//only for initializing


	vector<double> PreConditionerInv;
	PreConditionerInv.resize(NumOfEquations,0.0);
	//setting user data with everything that's required;
	UserData data(ReactionsPtr,&KineticsValueInfo,&RxnsWithRevKin,
		&Temp, &EffVolume, &InletPressure, &NumOfEquations,
		&SpeciesIndex, &StoichInfo, &RateInfo,AbsTolerance, 
		Rg, &surfaceSpecies, &AllSites, &InitialMolecFlows, 
		&SiteOccupants, &SiteBalanceInfo, &ReverseRxns, 
		&FirstDcldFlowSpecies, &FirstFoundProductFlowSpecies, &InletMassFlowRate,&MolWtOfSpecies, &ActualKineticsUsed, PreConditionerInv);

	/*cout<<"initially"<<endl;
	for (int i =0;i<data.Reactions->size();i++)
	{
		cout<<data.kinetics->at(i).getKineticsValue(Network->Temperature)<<endl;
	}
	cout<<"now"<<endl;

	realtype* r = &KineticsValueInfo[0];
	r[2]=r[2]*1.05;

	for (int i =0;i<data.Reactions->size();i++)
	{
		cout<<data.kinetics->at(i).getKineticsValue(Network->Temperature)<<endl;
	}*/

	
	cout<<"solving the kinetic model"<<endl;
	/* Allocate N-vectors. */
	yy = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);
	yp = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);
	avtol = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

	//sensavtol = N_VNew_Serial(NumOfEquations);
	//if(check_flag((void *)sensavtol, "N_VNew_Serial", 0)) return(1);

	id = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);




	
	yval = NV_DATA_S(yy);
	ypval = NV_DATA_S(yp);
	atval = NV_DATA_S(avtol);
	idval = NV_DATA_S(id);
	//sensAtval = NV_DATA_S(sensavtol);

	map<int,pair<double,int> >::iterator it;

	//Initialize yval and idval
	//yval is how yy is manipulated. idval is how id is manipulated to tell IDA if something is a differential or algebraic variable or not

	for (int i=0;i<NumOfEquations-1;i++)
	{
		//if (surfaceSpecies.count(i)>0) yval[i]=RCONST(0.0001);
		yval[i]= RCONST(0.0);
		//for a PFR, if the species is a surface species (including site itself) then it is an algebraic variable
		if ((Type==0 && (surfaceSpecies.count(i)>0)) /*|| (Type==1 && i==FirstFoundProductFlowSpecies)*/)
		{
			idval[i] = RCONST(0.0);//if surface species (which also accounts for sites), then set id = 0.0. also set to zero for the last equation (for vol flow rate)
			if (Type==0)NeedToCheckIC = true;
			
		}
		else idval[i] = RCONST(1.0);

	}

	idval[NumOfEquations-1]=RCONST(0.0); // zero for PFR and for CSTR
	yval[NumOfEquations-1] = RCONST(0.0);

	
	double InitialyvalSum = 0.0;

	//calculating the total initial value (flow rate - to get hte initial estimate of vol flow rate)
	for (it=InitialMolecFlows.begin();it!=InitialMolecFlows.end();it++)
	{
		//cout<<"setting initial value "<<it->second.first<<endl;
		yval[it->first]=it->second.first;
		//cout<<yval[it->first]<<endl;
		if (it->second.second==0)// if zero, then molecular flows		
			InitialyvalSum+=yval[it->first];//adding all molecular flows
	
	}

	//calculating the volumetric flowrate
	yval[NumOfEquations-1]=InitialyvalSum*Rg*(Temp)/(InletPressure);//initializing the volumetric flow rate! 

	
	
	
	//Above holds for PFR. We need to modify this for CSTRs because we work with concentrations there
	//Also we need to modify the initial conditions 

	//For now, we assume that the initial conditions require that the reactor be filled up with an inert
	//BIG HACK-- this inert is assumed to be nitrogen. The code will probably explode if N2 is not an initial reactant! 
	if (Type==1)
	{
		for (it=InitialMolecFlows.begin();it!=InitialMolecFlows.end();it++)
		{
			if (it->second.second==0)
			{
				if (*SpeciesIndex[it->first]=="N#N")
					yval[it->first]=InitialyvalSum/yval[NumOfEquations-1];//converting into concentrations
				else
					yval[it->first]=0.0;
			}
		}

	}

	

	//Initialize ypval - first set it to ZERO then, run Residual to get the new ypval.
	
	N_Vector ypvalInit, rvalInit;
	realtype *ypvalInitPtr, *rvalInitPtr;

	ypvalInit = rvalInit = NULL;
	ypvalInitPtr = rvalInitPtr = NULL;
	
	ypvalInit = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)ypvalInit, "N_VNew_Serial", 0)) return(1);
	
	rvalInit = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)rvalInit, "N_VNew_Serial", 0)) return(1);

	constraints = N_VNew_Serial(NumOfEquations);
	if(check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);

	N_VConst(RCONST(1.0), constraints);/*Set constraints to all 1's for nonnegative solution values */

	ypvalInitPtr = NV_DATA_S(ypvalInit);
	rvalInitPtr = NV_DATA_S(rvalInit);

	for (int i=0;i<NumOfEquations;i++)
	{
		ypvalInitPtr[i] = RCONST(0.0);
		rvalInitPtr[i] = RCONST(0.0);
	}

	//running residual functions
	if (Type==0)Residual(RCONST(0.0),yy,ypvalInit,rvalInit,&data);
	else DynamicCSTRresidual(RCONST(0.0),yy,ypvalInit,rvalInit,&data);
	//cout<<"ypval is "<<endl;
	ofstream ratefile("rate.txt");

	double ypvalSum = 0.0;
	
	//setting ypval;
	for (int i=0;i<NumOfEquations-1;i++)
	{
		if (surfaceSpecies.count(i)==0)//if not a site or a surface species -- surfaceSpecies includes both
		{
			if (Type==0)//if PFR
			{
				ypval[i] = -rvalInitPtr[i];//just rhe negative of the residual
				ypvalSum+=ypval[i];
			}
			else //if not PFR
			{
				ypval[i]= -rvalInitPtr[i]/EffVolume;//just the reverse divided by the volume
				ypvalSum+=ypval[i];

			}
			
		}
		else //if a surface species
		{
			if (Type==0)
				ypval[i]=RCONST(0.0); //arbitrarily setting to zero - does not matter anyways
			else
				ypval[i] = -rvalInitPtr[i]; 
		}

		
		//cout<< ypval[i]<<endl;

	}

	//for a CSTR site balance for sites is obtained by differentiating the original algebraic summation equation
	//Once we cacluate teh yp of all other species involved in teh site balance, we can run the residual function again to get the ypval of the sites
	if (Type==1)
	{
		set<int>::iterator setIt;

		for (setIt=AllSites.begin();setIt!=AllSites.end();setIt++)
		{
			ypval[*setIt]=-getSiteBalanceDynamicResidual(yp,*setIt,&data);
		}
	}

	//ypval[NumOfEquations-2]=-rvalInitPtr[NumOfEquations-2]*Rg*Temp/Volume;
	
	
	//if (Type==0)ypval[FirstFoundProductFlowSpecies]=-ypvalSum;
	
	if (!isSolvingIC) 
	{
		for (int i=0;i<NumOfEquations-1;i++)
			ratefile<<*SpeciesIndex[i]<<"  "<<ypval[i]<<endl;
	}

	//cout<<"rvalInitPtr[NumOfEquations-1] is "<<rvalInitPtr[NumOfEquations-1]<<endl;
	//ypval[NumOfEquations-1]=RCONST(0.0);


	//cout<<"before destroying ypvalInit"<<endl;
	
	
	N_VDestroy_Serial(ypvalInit);

	
	N_VDestroy_Serial(rvalInit);

	
	for (int i=0;i<NumOfEquations;i++)
	{
		if (idval[i]!=0)
			atval[i] = RCONST(AbsTolerance);
		else
			atval[i] = RCONST(AbsTolerance2);
		//sensAtval[i] = RCONST(sensAbsTol);
	}	
	rtol = RCONST(1.0e-4);
	sensrTol = RCONST(1.0e-5);

	t0 = RCONST(0.0);

	//cout<<"after destroying ypvalInit"<<endl;
	
	/*What follows basically is set up of the IDA related inputs */
	/* Call IDACreate and IDAMalloc to initialize IDA memory */
	
	mem = IDACreate();
	if(check_flag((void *)mem, "IDACreate", 0)) return(1);

	retval = IDASetUserData(mem, &data);
	if(check_flag(&retval, "IDASetUserData", 1)) return(1);

	retval = IDASetId(mem, id);
	if(check_flag(&retval, "IDASetId", 1)) return(1);

	if (Type==1)
	{
		retval = IDASetSuppressAlg(mem,true);
		if (check_flag(&retval, "IDASetSuppressAlg", 1)) return (1);
	}

	if (Type==0)
	{
		retval = IDAInit(mem, &Residual, t0, yy, yp);
		if(check_flag(&retval, "IDAInit", 1)) return(1);
	}
	else
	{
		retval = IDAInit(mem, &DynamicCSTRresidual, t0, yy, yp);
		if(check_flag(&retval, "IDAInit", 1)) return(1);
	}

	retval = IDASVtolerances(mem, rtol, avtol);
	if(check_flag(&retval, "IDASVtolerances", 1)) return(1);

	/* Free avtol */
	N_VDestroy_Serial(avtol);

	//retval = IDASetConstraints(mem, constraints);
	//if(check_flag(&retval, "IDASetConstraints", 1)) return(1);
	N_VDestroy_Serial(constraints);/* free constraints */

	/* Call IDADense and set up the linear solver. */
	retval = IDADense(mem, NumOfEquations);
	if(check_flag(&retval, "IDADense", 1)) return(1);

	/*Call IDASpgmr */

	//retval = IDASpgmr(mem, 0);
	//if(check_flag(&retval, "IDASpgmr", 1)) return(1);

	//retval = IDASptfqmr(mem, 0);
    //if(check_flag(&retval, "IDASptfqmr", 1)) return(1);

	/* Specify preconditioner */
   // retval = IDASpilsSetPreconditioner(mem, &PreCondSetup, &PreCondSolve);
   // if(check_flag(&retval, "IDASpilsSetPreconditioner", 1)) return(1);
  
	//Can set jacobian if required
	retval = IDADlsSetDenseJacFn(mem, &Jacobian);
	if(check_flag(&retval, "IDADlsSetDenseJacFn", 1)) return(1);

	retval = IDASetMaxNumSteps (mem, 10000);
	if(check_flag(&retval, "IDASetMaxNumSteps",1)) return (1);

	
	if (Type==0)tfinal = Volume;
	else tfinal = EffVolume/yval[NumOfEquations-1]*5;//5 times the residence time -- more than sufficient to reach steady state

	tout1=RCONST(tfinal)/RCONST(10.0); 
	tout = tout1;



	//cout<<"here 1"<<endl;

	vector<double> pbar;
	for (int i =0;i<NumOfParameters;i++)
	{
		if (useAlternativeKinetics)
			pbar.push_back(abs(AlternativeKinParams[i])+1.0);//abs ensures they are always positive, adding 1.0 ensures they are at least 1.0
		else
			pbar.push_back(abs(KineticsValueInfo[i])+1.0);// KineticsValueInfo will always be nonnegative; adding a 1.0 because I don't want this to ever be zero. 

	}



	
	
	//sensitivity analysis
	if (doSensitivity)
	{


		yS = N_VCloneVectorArray_Serial(NumOfParameters, yy);
		if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
		for (int is=0;is<NumOfParameters;is++) N_VConst(RCONST(0.0), yS[is]);

		//Ith(yS[0],3)=RCONST(1.0);

		ypS = N_VCloneVectorArray_Serial(NumOfParameters, yy);
		if (check_flag((void *)ypS, "N_VCloneVectorArray_Serial", 0)) return(1);
		for (int is=0;is<NumOfParameters;is++) N_VConst(RCONST(0.0), ypS[is]);

		/*sensavtol = N_VCloneVectorArray_Serial(NumOfParameters, yy);
		if (check_flag((void *)ypS, "N_VCloneVectorArray_Serial", 0)) return(1);
		for (int is=0;is<NumOfParameters;is++) N_VConst(RCONST(sensAbsTol), sensavtol[is]);*/

		flag = IDASensInit(mem, NumOfParameters, IDA_STAGGERED, NULL, yS, ypS);
		if(check_flag(&flag, "IDASensInit", 1)) return(1);

		//flag = IDASensSVtolerances(mem, sensrTol, sensavtol);
	   // if(check_flag(&flag, "IDASensSVtolerances", 1)) return(1);

		flag = IDASensEEtolerances(mem);
		if(check_flag(&flag, "IDASensEEtolerances", 1)) return(1);

		flag = IDASetSensErrCon(mem, false);
		if (check_flag(&flag, "IDASetSensErrCon", 1)) return(1);

	
		
		if (!useAlternativeKinetics)
		{	
			if (KineticsValueInfo.size()==0)cout<<"Oh no! the parameter list is empty! This is going to make IDAS's sensitivity analysis go crazy!!"<<endl;
			flag = IDASetSensParams(mem, &KineticsValueInfo[0], &pbar[0], NULL);//passing the first element of KineticsValueInfo is safe in this context - KineticsValueInfo's size is greater than zero and its size is going to fixed! 
			if (check_flag(&flag, "IDASetSensParams", 1)) return(1);
		}
		else
		{
			flag = IDASetSensParams(mem, &AlternativeKinParams[0], &pbar[0], NULL);//passing the first element of KineticsValueInfo is safe in this context - KineticsValueInfo's size is greater than zero and its size is going to fixed! 
			if (check_flag(&flag, "IDASetSensParams", 1)) return(1);

		}

	}
	

	//cout<<"here 2"<<endl;

	//if initial conditions need to be calculated, teh following if block is executed. Initial conditions are required of PFR models requiring site balance
	if (NeedToCheckIC)
	{
		/* TODO: if initial conditions are to be calculated by solving an initial dynamic CSTR to get some initial guess. Check TODO further down
		
		if (Type==0) Solve(1,true);

		for (int i =0;i<NumOfEquations;i++)
		{
			yval[i]=InitialConditionsFromDyCSTR[i];
			ypval[i]=InitialConditionsDerivativeFromDyCSTR[i];
		}*/



		retval = IDASetMaxNumStepsIC(mem, 20);
		if (check_flag(&retval, "IDASetMaxNumStepsIC", 1)) return(1);

		retval = IDASetMaxNumJacsIC(mem,20);
		if (check_flag(&retval, "IDASetMaxNumStepsIC", 1)) return(1);

		retval = IDASetMaxNumItersIC(mem,20);
		if (check_flag(&retval, "IDASetMaxNumItersIC", 1)) return(1);

		tout2 = tout1/RCONST(1000.0);
		retval = IDACalcIC(mem, IDA_YA_YDP_INIT,tout2);
		if(check_flag(&retval, "IDACalcIC", 1))
		{
			long int lstflag,jevals,resForJevals;
			int flag = IDADlsGetLastFlag(mem, &lstflag); 
			cout<<"last flag "<<lstflag<<endl;

			flag = IDADlsGetNumJacEvals(mem, &jevals);
			cout<<"jevals is "<<jevals;
			flag = IDADlsGetNumResEvals(mem, &resForJevals);
			cout<<"res evals for Jac is "<<resForJevals<<endl;

			return(1);
		}
		

	//	cout<<"here 3"<<endl;

		flag = IDAGetConsistentIC(mem, yy, yp);
		if (check_flag(&flag, "IDAGetConsistentIC", 1)) return(1);

	//	cout<<"here 4"<<endl;

		if (doSensitivity)
		{

			flag = IDAGetSensConsistentIC(mem, yS, ypS);
			if (check_flag(&flag, "IDAGetSensConsistentIC", 1)) return(1);
		}
	}


	//cout<<"here 5"<<endl;


    


	ofstream KineticOutputFile;
	ofstream KineticSensOutputFile;
	ofstream KineticsDORCFile;
	
	string kineticFileName(filename);

	if (!isSolvingIC)
	{
	
		KineticOutputFile.open(filename);
		KineticOutputFile<<"Output of the kinetic model at different times"<<endl;
		KineticOutputFile.close();
		printOutputs(RCONST(0.0),yy,yp);
		if (shouldCalculateRates)
			Rates.push_back(getrateVector(yy,&data));

			//Sensitivity output
	
	
		if (doSensitivity)
		{
			KineticSensOutputFile.open("sensitivities.txt");
			KineticSensOutputFile<<"Output of sensitivities at different times"<<endl;
			KineticSensOutputFile.close();
			PrintSensOutput(RCONST(0.0),yS);
			
			if (!useAlternativeKinetics)
			{
				KineticsDORCFile.open("DORC.txt");
				KineticsDORCFile<<"Output of DORC at different times"<<endl;
				KineticsDORCFile.close();			
				PrintDORC(RCONST(0.0),ypS,yp);
			}
			else
			{
				KineticsDORCFile.open("DORC.txt");
				KineticsDORCFile<<"DORC is not evaluated for this option of kinetics specification"<<endl;
			}
			
		}
	}
	
	
	


	//I'm solving it in steps. 

	

	while (tout<=1.01*tfinal)
	{
		//cout<<"solving for the next 10%"<<endl;
		retval = IDASolve(mem, tout, &tret, yy,yp, IDA_NORMAL);
		if(check_flag(&retval, "IDASolve", 1)) return(1);
		
		if (!isSolvingIC)
		{
			printOutputs(tret,yy,yp);
			if (shouldCalculateRates)
				Rates.push_back(getrateVector(yy,&data));
		}

		//Sensitivity part
		if (doSensitivity)
		{
			retval = IDAGetSens(mem, &tret, yS);
			if (check_flag(&retval, "IDAGetSens", 1)) break;
			retval = IDAGetSensDky(mem,tret,1,ypS);
			if (check_flag(&retval, "IDAGetSensDky", 1)) break;
			PrintSensOutput(tret, yS);	
			if (!useAlternativeKinetics)PrintDORC(tret,ypS,yp);
		}

		if (retval == IDA_SUCCESS)
		{
			//if (!isSolvingIC)cout<<"kinetic model solved up to "<<tout<<" ("<<tout/RCONST(tfinal)*100<<" %)"<<endl;
			//else cout<<"Initial conditions calculation solved up to "<<tout<<" ("<<tout/RCONST(tfinal)*100<<" %)"<<endl;
			//retval = IDAReInit(mem, tout, yy, yp);
			//if(check_flag(&retval, "IDAReInit", 1)) return(1);
			tout += tout1;	
		}
	}

	/*TODO this is to save the initial conditions guesses which are basically the end result of a dynamic CSTR run
	if (isSolvingIC)
	{
		InitialConditionsFromDyCSTR.resize(NumOfEquations);
		InitialConditionsDerivativeFromDyCSTR.resize(NumOfEquations);
		for (int i =0;i<NumOfEquations;i++)
		{

			//InitialConditionsFromDyCSTR[i]=yval[i];
			//InitialConditionsDerivativeFromDyCSTR[i] = ypval[i];

			
			if (surfaceSpecies.count(i)>0 || NumOfEquations-1 ==i || NumOfEquations-2 ==i)
			{
				InitialConditionsFromDyCSTR[i]=yval[i];
				InitialConditionsDerivativeFromDyCSTR[i] = ypval[i];
			}
			else
			{
				InitialConditionsFromDyCSTR[i]=yval[i]*yval[NumOfEquations-1];
				InitialConditionsDerivativeFromDyCSTR[i] = ypval[i]*EffVolume;
			}

		}
		

	}*/
	
	if (!isSolvingIC) cout<<"kinetic model outputs written in to "<<kineticFileName<<endl;

	//if (doSensitivity)
	storeOutputs(yy,yp,yS);

	
	//freeing allocated memory
	IDAFree(&mem);
	N_VDestroy_Serial(yy);
	N_VDestroy_Serial(yp);
	N_VDestroy_Serial(id);
//	N_VDestroy_Serial(sensavtol);

	if (doSensitivity)
	{
		N_VDestroyVectorArray_Serial(yS, NumOfParameters);
		N_VDestroyVectorArray_Serial(ypS, NumOfParameters);
		N_VDestroyVectorArray_Serial(sensavtol, NumOfParameters);
	}


    return(0);
}

void KineticModel::storeOutputs(N_Vector y, N_Vector yp, N_Vector* yS)
{
	double* yval;
	double* ypval;
	double* ysval;

	yval = NV_DATA_S(y);
	ypval = NV_DATA_S(yp);//does NOT get used. 

	for (int i=0;i<IndicesOfOutputs.size();i++)
	{
		OutputValues[*SpeciesIndex[IndicesOfOutputs[i]]]=yval[IndicesOfOutputs[i]];
		if (useAlternativeKinetics && doSensitivity)
		{
			
			vector<double> sensitivityOfSpeciesI;
			for (int j=0;j<NumOfParameters;j++)
			{
				ysval = NV_DATA_S(yS[j]);
				//sensitivityOfSpeciesI.push_back(Ith(yS[j],IndicesOfOutputs[i]));
				sensitivityOfSpeciesI.push_back(ysval[IndicesOfOutputs[i]]);
			}

			SensitivityOutputs[*SpeciesIndex[IndicesOfOutputs[i]]] = sensitivityOfSpeciesI;
		}
	}
			

}


int KineticModel::getStoichCoeff(int mol, int rxn)
{
	if (mol>=NumOfEquations-1 || rxn>=Network->AllReactions.size())
		return 0;
	else
		return Network->AllReactions[rxn].occurence(SpeciesIndex[mol]);
}



void KineticModel::GenerateStoichAndRateVectorInfo( map<string,pair<double,int> >& InitialFlowRates, multimap<string,string>& SiteBalance)
{
    map<string*,int,classcomp>::iterator it;
	for (it =MoleculeSetPtr->begin();it!=MoleculeSetPtr->end();it++)
	{
		
		
		SpeciesIndex.insert(pair<int,string*>(SpeciesIndex.size(),it->first));
		
		/* if the string satisfies constraints specified for required outputs, the indices are to be automatically updated */

		if (isRequiredOutput(*(it->first)))
			IndicesOfOutputs.push_back(SpeciesIndex.size()-1);
		

		//Calculate molecular weight of each molecule! 
		Molecule mol(*(it->first),moleculesize(*(it->first)));
		MolWtOfSpecies[SpeciesIndex.size()-1]=mol.MolecularWeight();
		

		
		//pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator> pairIter1;
		
		SetStoichInfoFromMap(MolReactantMapPtr, SpeciesIndex.size()-1);
		SetStoichInfoFromMap(MolProductMapPtr, SpeciesIndex.size()-1);

		
		
		SetRateInfo(SpeciesIndex.size()-1);


		bool isInitialMol = false;
	
		if (InitialFlowRates.count(*(it->first))>0)//since initial reactants are not lumped, I do not have to worry about converting them to lumps! 
		{
			isInitialMol = true;
			pair<double,int> P = InitialFlowRates[*(it->first)];
			InitialMolecFlows.insert(pair<int,pair<double,int> >(SpeciesIndex.size()-1,P));
			if (P.second==1)
			{
				AllSites.insert(SpeciesIndex.size()-1);//add into AllSites if the string is one of the heterogeneous composite sites	
				multimap<string,string>::iterator mm_it;
				pair<multimap<string,string>::iterator, multimap<string,string>::iterator> mmpair_it;

				mmpair_it = SiteBalance.equal_range(*(it->first));
				
				for (mm_it=mmpair_it.first;mm_it!=mmpair_it.second;mm_it++)
				{
					SiteBalanceInfo.insert(pair<int,string>(SpeciesIndex.size()-1,mm_it->second));
				}

			
			}
			
			if (P.second==0)
			{
				if (FirstDcldFlowSpecies==-1)
					FirstDcldFlowSpecies = SpeciesIndex.size()-1;
				//incrementing inlet mass flow rate
				InletMassFlowRate+=P.first*MolWtOfSpecies[SpeciesIndex.size()-1];
			}
			

		}

		bool isSurface = false;
				
		for (int i =0;i<Network->CompositeSites.size();i++)
		{
			if (Network->CompositeSites[i].second == Heterogeneous)
			{
				if (("[" + Network->CompositeSites[i].first + "]" == *(it->first))  || PopulateSiteOccupancyMap(*(it->first), SpeciesIndex.size()-1, Network->CompositeSites[i].first))
				{
					surfaceSpecies.insert(SpeciesIndex.size()-1);//if there exists site atoms of any kind, insert species into surfaceSpecies
					isSurface = true;
				}
			}
		}

		if (!isSurface && !isInitialMol)
			FirstFoundProductFlowSpecies = SpeciesIndex.size()-1;

		

		

		//cout<<"Inlet mass flow rate is "<<InletMassFlowRate<<endl;


		
	}
	


	
}

void KineticModel::SetRateInfo(int mol)
{
	multimap<string*,int>::iterator it;
	pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator> pairIter;
	pairIter = MolReactantMapPtr->equal_range(SpeciesIndex[mol]);
	
	for (it = pairIter.first;it!=pairIter.second;it++)
	{
		RateInfo[it->second-1].insert(mol);
	}
}

void KineticModel::SetStoichInfoFromMap(std::multimap<string*,int> * Map, int mol)
{

	pair<multimap<string*,int>::iterator, multimap<string*,int>::iterator> pairIter;
	pairIter = Map->equal_range(SpeciesIndex[mol]);
	multimap<string*,int>::iterator it;

	for (it = pairIter.first;it!=pairIter.second;it++)
	{
		int rxn = it->second-1;//rxns in mol(Reactant/Product)map 
		int stoich = ReactionsPtr->at(rxn).NetOccurence(SpeciesIndex[mol]);//this way even if net occurence is zero, stoich info will be updated! 
		StoichInfo.insert(pair<int, pair<int,int> >(mol, pair<int,int>(rxn,stoich)));
	}
	
	
}


int KineticModel::check_flag(void *flagvalue, char *funcname, int opt)
{
  
	int *errflag;
	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, 
			    "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
				funcname);
		return(1);
	} else if (opt == 1) {
	/* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
		  fprintf(stderr, 
			      "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
				  funcname, *errflag);
		 return(1); 
    }
	} else if (opt == 2 && flagvalue == NULL) {
	/* Check if function returned NULL pointer - no memory allocated */
		fprintf(stderr, 
			        "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
				  funcname);
		return(1);
	}

  return(0);
}





int KineticModel::getNumEquations()
{
	return NumOfEquations;
}


bool KineticModel::isRequiredOutput(string S)
{
	Molecule mol(S, moleculesize(S));
	mol.unique_smiles();
	if (OutputConstraintsSpecified)
		return (RequiredOutputs)(mol);
	else
		return (specifiedOutputs.count(mol.moleculestring())>0);
}

void KineticModel::printOutputs(double tret, N_Vector y, N_Vector yp)
{
	ofstream kineticFile;
	double* yval;
	double* ypval;

	yval = NV_DATA_S(y);
	ypval = NV_DATA_S(yp);

	kineticFile.open(filename, ios::app);//open in append mode! 
	
	kineticFile<<endl;

	kineticFile<<"At "<<(tret*100/tfinal)<<" % of the total volume/time"<<endl;

	for (int i= 0;i<IndicesOfOutputs.size();i++)
	{
		double massFlow = 0.0;

		if (surfaceSpecies.count(IndicesOfOutputs[i])==0)
			massFlow = yval[IndicesOfOutputs[i]]*MolWtOfSpecies[IndicesOfOutputs[i]];
		kineticFile<<(*SpeciesIndex[IndicesOfOutputs[i]])<<"  "<<yval[IndicesOfOutputs[i]]<<" ("<<ypval[IndicesOfOutputs[i]]<<")"<<"  "<<massFlow<<endl;
	}
	
	kineticFile<<"Volumetric flow rate is "<<yval[NumOfEquations-1]<<endl;
}

void KineticModel::PrintSensOutput(double tret, N_Vector* y)
{
	ofstream sensitivityFile;
		
	sensitivityFile.open("sensitivities.txt", ios::app);//open in append mode! 
	
	sensitivityFile<<endl;

	sensitivityFile<<"At "<<(tret*100/Volume)<<" % of the total volume"<<endl;

	for (int i= 0;i<IndicesOfOutputs.size();i++)
	{
		sensitivityFile<<(*SpeciesIndex[IndicesOfOutputs[i]]);
		for (int j = 0;j<NumOfParameters;j++)
		{
			sensitivityFile<<"  "<<Ith(y[j],IndicesOfOutputs[i]);
		}
		sensitivityFile<<endl;

	}
	sensitivityFile<<"Nu ";
	for (int j = 0;j<NumOfParameters;j++)
	{
		sensitivityFile<<"  "<<Ith(y[j],NumOfEquations-1);
	}
	sensitivityFile<<endl;
	
}

void KineticModel::PrintDORC(double tret, N_Vector* ypS, N_Vector yp)
{
	ofstream DORCFile;
	double* ypval;

	ypval = NV_DATA_S(yp);
		
	DORCFile.open("DORC.txt", ios::app);//open in append mode! 
	
	DORCFile<<endl;

	DORCFile<<"At "<<(tret*100/Volume)<<" % of the total volume"<<endl;

	DORCFile<<"Rxn";
	for (int i= 0;i<IndicesOfOutputs.size();i++)
		DORCFile<<"  "<<(*SpeciesIndex[IndicesOfOutputs[i]]);

	DORCFile<<endl;

	
	for (int j = 0;j<NumOfParameters;j++)
	{
		DORCFile<<j<<"  "<<endl;
		for (int i= 0;i<IndicesOfOutputs.size();i++)
		{
			double DORC = (Ith(ypS[j],IndicesOfOutputs[i]))*KineticsValueInfo[j]/(ypval[IndicesOfOutputs[i]]+RCONST(1.0e-20));
			
			//We also have to add the terms corresponding to the reverse steps so that thermodynamic consistency is maintained

			pair<multimap<int,int>::iterator, multimap<int,int>::iterator > pairIter;
			pairIter = ReverseRxns.equal_range(j);

			for (multimap<int,int>::iterator it = pairIter.first;it!=pairIter.second;it++)
			{
				//cout<<"Adding reverse reactions"<<endl;
				DORC+=(Ith(ypS[it->second],IndicesOfOutputs[i]))*KineticsValueInfo[it->second]/(ypval[IndicesOfOutputs[i]]+RCONST(1.0e-20));
			}

			//DORCFile<<" "<<DORC<<" (Rxn: "<<j<<", Sensi: "<<Ith(ypS[j],IndicesOfOutputs[i])<<")";
			DORCFile<<"  "<<DORC;
		}
		DORCFile<<endl;
	
		
		//cout<<Ith(yS[j],IndicesOfOutputs[i])<<" "<<KineticsInfo[j]<<" "<<yval[IndicesOfOutputs[i]]<<endl;
	}
	

	
	/*DORCFile<<"Nu ";
	for (int j = 0;j<NumOfParameters;j++)
	{
		DORCFile<<"  "<<Ith(yS[j],NumOfEquations-1)*KineticsInfo[j]/(yval[NumOfEquations-1]+RCONST(1.0e-20));
	}
	DORCFile<<endl;*/
}





bool KineticModel::SetKinetics()
{
	ofstream kineticsFile;
	kineticsFile.open("kineticsvalue.txt");
	if (useAlternativeKinetics)
		kineticsFile<<"This file consists of, for each reaction, the kinetic value, rule name, and (if calculated as a reverse of a reaction) the enthalpy, entropy, and delN."<<endl;
	else
		kineticsFile<<"This file consists of, for each reaction, the kinetic value, rule name, pre exponential factor, activation barrier, and (if calculated as a reverse of a reaction) the enthalpy, entropy, and delN."<<endl;

	

	for (int i=0;i<ReactionsPtr->size();i++)
	{
		
		double PreExp, ActE, TempIndex, kinValue;
		int rule = (*ReactionsPtr)[i].get_rule();
		
		
		PreExp = 0.0; ActE=0.0; TempIndex = 0.0;kinValue = 0.0;
		bool calcK = false; bool usesBEP = false; bool usesLFER = false;
		double alpha = 0.0; double beta = 0.0; bool stick = false;

		double dH = Network->calculateThermoOfRxn(EnthalpyType,(*ReactionsPtr)[i],Temp)*1000.0;
		double dS = Network->calculateThermoOfRxn(EntropyType,(*ReactionsPtr)[i],Temp);
		double delN = Network->getDeltaNGasPhase((*ReactionsPtr)[i]);

		int kIndex =-1; int EaIndex = -1; double kRelative = 0.0; double EaRelative = 0.0; double RefTemp = 0.0;
		if (useAlternativeKinetics)
		{
			KineticsInfo kinInfo(*AlternativeKinetics, &AlternativeKinParams ,(*ReactionsPtr)[i],Rg,Temp, dH, dS, delN);
			
			if (!kinInfo.getKineticParameters(kinValue, kIndex, kRelative,EaIndex, EaRelative,RefTemp, TempIndex, calcK))
				return false;

			KineticsEstimates kEstimates;
			kEstimates.kParam = kIndex;
			kEstimates.kRelative = kRelative;
			kEstimates.EaParam = EaIndex;
			kEstimates.EaRelative = EaRelative;
			kEstimates.kRefTemp = RefTemp;
			kEstimates.isReverse = calcK;
			kEstimates.TempExp = TempIndex;

			KineticsEvaluationFromParam kin(kEstimates, &AlternativeKinParams);
			kinValue = kin.getKineticsValue(Temp);


			ActualKineticsUsed.push_back(kin);		

		}
		else
		{
			KineticsInfo kinInfo(*kinetics,(*ReactionsPtr)[i],Rg,Temp, dH, dS, delN, 0.01, Network->firstGasSpecies(i,0), Network->firstGasSpecies(i,1)); 


			if (!kinInfo.getKineticParameters(PreExp,ActE,TempIndex,kinValue,calcK, usesBEP, usesLFER, alpha, beta, stick))
				return false;

			//Here, I am using the KineticsValueInfo[i] that has the appropriate kinetics. 
			//So, I am just creating pretty much a dummy KineticsEvaluationFromParam object
			// by storing kParam as i, I access the correct kinetics, which is, KineticsValueInfo[i]! 
			KineticsEstimates kEstimates;
			kEstimates.kParam = i;
			kEstimates.kRelative = 1.0;
			kEstimates.EaParam = -1;
			kEstimates.EaRelative = 0.0;
			kEstimates.kRefTemp = Temp; //setting the ref temp arbitrarily here.
			kEstimates.isReverse = calcK;
			kEstimates.TempExp = 0.0;
			ActualKineticsUsed.push_back(KineticsEvaluationFromParam(kEstimates,&KineticsValueInfo));
		}
	
			
		if (calcK)
		{
			//taking these directly from kinInfo's implementation - somewhat repetitive - needs clean up!
			double deltaG = dH - dS*Temp;
			double K = exp(-deltaG/(RCONST(Temp*8.314)));
			K=K*pow(Rg*Temp,-delN);
			RxnsWithRevKin[i]=K;//this is the correction that is multiplied to the kinetics in the reverse direction to get our rate for this reaction
			//NOTE: right now, we don't use this K - I am just lugging it around MAYBE remove this  and just keep a set for RxnsWithRevKin
		}

		KineticsValueInfo[i]= log(kinValue);
		

		
		if (!useAlternativeKinetics)
		{
			kineticsFile<<kinValue<<"   "<<Network->Rtlist[Network->AllReactions[i].get_rule()].getRuleName()<<"  "<<PreExp<<"  "<<ActE<<"  ";
			if (calcK)
				kineticsFile<<"calculated as reverse with dH: "<<dH<<"  and dS: "<<dS<<"  and delN: "<<delN;
			kineticsFile<<endl;
		}
		else
		{
			kineticsFile<<kinValue<<"   "<<Network->Rtlist[Network->AllReactions[i].get_rule()].getRuleName();
			if (calcK)
				kineticsFile<<"  kRelative: "<<kRelative<<"; calculated as reverse with dH: "<<dH<<"  and dS: "<<dS<<"  and delN: "<<delN;
			kineticsFile<<endl;
		}
		
		
	}

	



	

	kineticsFile.close();
	

	for (map<int,double>::iterator rev_it = RxnsWithRevKin.begin();rev_it!=RxnsWithRevKin.end();rev_it++)
	{
		findReverseReactions(rev_it->first);
	}

	


	return true;
}

bool KineticModel::PopulateSiteOccupancyMap(string mol, int Index, string site)
{
	if (!Network->CompositeSites.empty())
	{
		Molecule m(mol, moleculesize(mol));
		int NumMatches = Patternmatch(m,Substructure(site,patternsize(site)),0).number_of_matches();
		if (NumMatches>0)
		{
			SiteOccupants.insert(pair<string, pair<int,int> >(site, pair<int,int>(Index, NumMatches)));
			
			return true;
		}		
	}
	return false;

}

void KineticModel::findReverseReactions(int rxn)
{
	for (int i=0;i<Network->AllReactions.size();i++)
	{
		if (Network->AllReactions[rxn].isReverseReaction(Network->AllReactions[i]) && abs(KineticsValueInfo[i]-KineticsValueInfo[rxn]/RxnsWithRevKin[rxn])<0.0001)
		{
			ReverseRxns.insert(pair<int,int>(rxn,i));
			ReverseRxns.insert(pair<int,int>(i,rxn)); //flipping the order and adding too 
		}
	}
}


void KineticModel::calculateRates()
{
	shouldCalculateRates = true;
}

	
/*

int KineticModel::CreateCHEMKINfileAndWriteSpecies()
{
	ofstream file("CHEMKINsurf.txt");
	ofstream file2("CHEMKINspecies.txt");
	ofstream file3("CHEMKINgas.txt");

	FindAllElements();//find all the elements in the system
	file3<<"ELEMENTS"<<endl;

	//cout<<AllElements.size()<<endl;
	for (set<string>::iterator it = AllElements.begin();it!=AllElements.end();it++)
		file3<<"  "<<*it;
	file3<<endl;
	file3<<"END"<<endl;
	file3<<"SPECIES"<<endl;


	

	set<int>::iterator it;
	
	map<string*,int,classcomp>::iterator molIter;

	map<string, int> IdenticalCHEMKINstrings;

	for (molIter = Network->AllMolecules.begin();molIter != Network->AllMolecules.end();molIter++)
	{
		
		Molecule mol(*(molIter->first), moleculesize(*(molIter->first)));
		mol.unique_smiles();
		string CHEMKINstring = mol.GetMF();
		
		if (IdenticalCHEMKINstrings.count(CHEMKINstring)>0)
		{
			IdenticalCHEMKINstrings[CHEMKINstring]+=1;
			CHEMKINstring+="-"+IntToStr(IdenticalCHEMKINstrings[CHEMKINstring]);
		}			
		else IdenticalCHEMKINstrings[CHEMKINstring]=1;
		CHEMKINSpeciesName[molIter->first] = CHEMKINstring;
		file2<<CHEMKINstring<<" --> "<<*molIter->first;

		//write thermo of each species into file 2 (enthalpy, entropy, and Cp at specified temps
		//TempsForCpCalc should be user specified. If not, calculate at 298, 400, 500, 600, 800, 1000K

		if (TempsForCpCalc.size()==0)
		{
			TempsForCpCalc.insert(298);
			TempsForCpCalc.insert(400);
			TempsForCpCalc.insert(500);
			TempsForCpCalc.insert(600);
			TempsForCpCalc.insert(800);
			TempsForCpCalc.insert(1000);
		}
		
		double dH=0.0;
		double dS = 0.0;
		bool isSiteMol = false;

		if (Network->isSiteIntermediate(*molIter->first))isSiteMol = true;

		ThermoGA::calculateDeltaH(mol,isSiteMol,298.0,dH);
		ThermoGA::calculateDeltaS(mol,isSiteMol,298.0,dS);
		file2<<"  "<<dH/4.18;
		file2<<"  "<<dS/4.18;
		
		for (set<int>::iterator SetIt =TempsForCpCalc.begin();SetIt!=TempsForCpCalc.end();SetIt++)
		{
			double Cp = 0.0;

			ThermoGA::calculateCp(mol,isSiteMol,*SetIt,Cp);
			file2<<"  "<<Cp/4.18;
		}
		file2<<endl;		
	}

	for (map<int,string*>::iterator map_it=SpeciesIndex.begin();map_it!=SpeciesIndex.end();map_it++)
	{
		if (surfaceSpecies.count(map_it->first)==0 && AllSites.count(map_it->first)==0)
			file3<<CHEMKINSpeciesName[map_it->second]<<endl;
	}
	file3<<"END"<<endl;


	for (it =AllSites.begin();it!=AllSites.end();it++)
	{
	
		string siteNameWithoutSqBrackets = (*SpeciesIndex[*it]).substr(1,SpeciesIndex[*it]->length()-2);
		
		if (SiteOccupants.count(siteNameWithoutSqBrackets)>0)
		{
			file<<"SITE/SURFACE/      SDEN/      /"<<endl;
			multimap<string, pair<int,int> >::iterator map_it;
			pair<multimap<string, pair<int,int> >::iterator, multimap<string, pair<int,int> >::iterator > pairIter;

			pairIter = SiteOccupants.equal_range(siteNameWithoutSqBrackets);

			for (map_it=pairIter.first;map_it!=pairIter.second;map_it++)
			{

				CHEMKINSpeciesName[SpeciesIndex[map_it->second.first]]+="(S)";
				file<<CHEMKINSpeciesName[SpeciesIndex[map_it->second.first]]<<"/"<<map_it->second.second<<"/"<<endl;
			}
		}
			
		CHEMKINBulkName[SpeciesIndex[*it]] = CHEMKINSpeciesName[SpeciesIndex[*it]]+"(B)";
		CHEMKINSpeciesName[SpeciesIndex[*it]]+="(S)";
		file<<CHEMKINSpeciesName[SpeciesIndex[*it]]<<endl;
		file<<"BULK "<<CHEMKINBulkName[SpeciesIndex[*it]]<<"/     /"<<endl;
		file<<"END"<<endl;
		
	}
	file<<"REACTIONS	MWOFF	KCAL/MOLE    "<<endl;
	file<<endl;
	file<<"!-----------------------------------------------------------------------------	"<<endl;
	file<<"! Reactions begin now and have been grouped according to their families"<<endl;
	file<<"!-----------------------------------------------------------------------------	"<<endl;
	file<<endl;

	file3<<"REACTIONS"<<endl;

	file.close();
	file2.close();
	file3.close();

	return 0;
}

*/

string KineticModel::GetMFofSpecies(int index)
{
	string molstr = *SpeciesIndex[index];
	Molecule mol(molstr, moleculesize(molstr));
	mol.unique_smiles();

	return mol.GetMF();
}

vector<vector<double> > KineticModel::getRates()
{
	return Rates;
}

/*
string KineticModel::generateCHEMKINRxnStrings(vector<Molecule>& reactants, vector<Molecule>& products)
{
	map<string,int> ReactantsStoich;
	map<string,int> ProductsStoich;

	for (int i =0;i<reactants.size();i++)
	{
		string molstring = reactants[i].moleculestring();
		string* molstrptr = StringRegistry::getStringPointer(molstring);
		string CHEMKINname = CHEMKINSpeciesName[molstrptr];
		if (ReactantsStoich.count(CHEMKINname)==0)
			ReactantsStoich[CHEMKINname]=-1;
		else ReactantsStoich[CHEMKINname]-=1;

		if (CHEMKINBulkName.count(molstrptr)==1)
		{
			string bulkName = CHEMKINBulkName[molstrptr];
			if (ProductsStoich.count(bulkName)==0)
				ProductsStoich[bulkName]=1;
			else ProductsStoich[bulkName]+=1;
		}
	}

	for (int i =0;i<products.size();i++)
	{
		string molstring = products[i].moleculestring();
		string* molstrptr = StringRegistry::getStringPointer(molstring);
		string CHEMKINname = CHEMKINSpeciesName[molstrptr];
		if (ProductsStoich.count(CHEMKINname)==0)
			ProductsStoich[CHEMKINname]=1;
		else ProductsStoich[CHEMKINname]+=1;

		if (CHEMKINBulkName.count(molstrptr)==1)
		{
			string bulkName = CHEMKINBulkName[molstrptr];
			if (ReactantsStoich.count(bulkName)==0)
				ReactantsStoich[bulkName]=-1;
			else ReactantsStoich[bulkName]-=1;
		}
	}

	string rxnstr="";

	map<string,int>::iterator it;

	for (it = ReactantsStoich.begin();it!=ReactantsStoich.end();it++)
	{
		if (it!=ReactantsStoich.begin())
			rxnstr+="+";
		if (it->second==-1)
			rxnstr+=it->first;
		else rxnstr=rxnstr+IntToStr(abs(it->second))+it->first;
	}
	rxnstr+="=";
	for (it = ProductsStoich.begin();it!=ProductsStoich.end();it++)
	{
		if (it!=ProductsStoich.begin())
			rxnstr+="+";
		if (it->second==1)
			rxnstr+=it->first;
		else rxnstr=rxnstr+IntToStr(abs(it->second))+it->first;
	}

	return rxnstr;	
}

void KineticModel::generateCHEMKINkineticsStrings(double A, double n, double Ea, bool stick)
{
	if (stick)RxnsWithStickingCoeff.insert(CHEMKINPreExpInfo.size());
	CHEMKINPreExpInfo.push_back(DoubleToString(A));
	CHEMKINIndexInfo.push_back(DoubleToString(n));
	CHEMKINActEInfo.push_back(DoubleToString(Ea));
	
}

int KineticModel::WriteRxnInfoInCHEMKINfile()
{
	ofstream file;
	ofstream file3;
	file.open("CHEMKINsurf.txt", ios::app);
	file3.open("CHEMKINgas.txt", ios::app);


	int maxRxnLength=0;
	for (int i =0;i<CHEMKINRxnStrings.size();i++)
	{
		int length = CHEMKINRxnStrings[i].length();
		if (maxRxnLength<=length)
			maxRxnLength = length;
	}
	int maxPreExpLength =0; 
	int maxIndexLength = 0; 
	int maxEaLength = 0;

	for (int i =0;i<CHEMKINPreExpInfo.size();i++)
	{
		int length = CHEMKINPreExpInfo[i].length();
		if (maxPreExpLength<=length)
			maxPreExpLength = length;
	}

	for (int i =0;i<CHEMKINIndexInfo.size();i++)
	{
		int length = CHEMKINIndexInfo[i].length();
		if (maxIndexLength<=length)
			maxIndexLength = length;
	}

	for (int i =0;i<CHEMKINActEInfo.size();i++)
	{
		int length = CHEMKINActEInfo[i].length();
		if (maxEaLength<=length)
			maxEaLength = length;
	}


	multimap<int, int>::iterator it;
	int prevRule = -1;

	int rxnCount =1;
	for (it=Network->ReactionsMap.begin();it!=Network->ReactionsMap.end();it++)
	{
		
		if (it->first!=prevRule)
		{
			file<<endl;
			file<<"!----------------------------------------------------------"<<endl;
			file<<"            !Reactions of rule ";
			file<<Network->Rtlist[it->first].getRuleName()<<":"<<endl;
			file<<"!----------------------------------------------------------"<<endl;
			file<<endl;

			prevRule = it->first;
		}

		
		
		string str(maxRxnLength+5-CHEMKINRxnStrings[(*it).second-1].length(), ' ');
		string str2(maxPreExpLength+3-CHEMKINPreExpInfo[(*it).second-1].length(), ' ');
		string str3(maxIndexLength+3-CHEMKINIndexInfo[(*it).second-1].length(), ' ');
		string str4(maxEaLength+3-CHEMKINActEInfo[(*it).second-1].length(), ' ');

		if (!isGasPhaseRxn(it->second-1))
		{
			WriteInfoFile(CHEMKINRxnStrings[(*it).second-1],file);
			WriteInfoFile(str,file);
			WriteInfoFile(CHEMKINPreExpInfo[(*it).second-1],file);
			WriteInfoFile(str2,file);
			WriteInfoFile(CHEMKINIndexInfo[(*it).second-1],file);
			WriteInfoFile(str3,file);
			WriteInfoFile(CHEMKINActEInfo[(*it).second-1],file);
			WriteInfoFile(str4,file);
			file<<"  !"<<rxnCount;
			file<<endl;
			if (RxnsWithStickingCoeff.count((*it).second-1)==1)
				file<<"STICK"<<endl;
		}
		else
		{
			WriteInfoFile(CHEMKINRxnStrings[(*it).second-1],file3);
			WriteInfoFile(str,file3);
			WriteInfoFile(CHEMKINPreExpInfo[(*it).second-1],file3);
			WriteInfoFile(str2,file3);
			WriteInfoFile(CHEMKINIndexInfo[(*it).second-1],file3);
			WriteInfoFile(str3,file3);
			WriteInfoFile(CHEMKINActEInfo[(*it).second-1],file3);
			WriteInfoFile(str4,file3);
			file<<"  !"<<rxnCount;
			file3<<endl;
		}
		rxnCount++;

	}
	file<<"END"<<endl;
	file3<<"END"<<endl;

	return 0;


}

void KineticModel::FindAllElements()
{
	for (set<string*,classcomp>::iterator set_it=Network->InitialReactants.begin();set_it!=Network->InitialReactants.end();set_it++)
	{

		//cout<<"initial reactant: "<<*(*set_it)<<endl;
		set<string> elementlist = Molecule(*(*set_it), moleculesize(*(*set_it))).GetElements();
		//cout<<elementlist.size()<<endl;

		for (set<string>::iterator it = elementlist.begin();it!=elementlist.end();it++)
			AllElements.insert(*it);
	}
	//cout<<AllElements.size()<<endl;
}

bool KineticModel::isGasPhaseRxn(int rxn)
{
	//if any of hte reactant is a surface species or a site, then the reaction is a surface reaction. 
	for (set<int>::iterator it = RateInfo[rxn].begin();it!=RateInfo[rxn].end();it++)
	{
		if (surfaceSpecies.count(*it)>0 || AllSites.count(*it)>0)//it is a surface species or a site
			return false;
	}
	return true;
}

*/


UserData::UserData(vector<generated_rxn>* rxns, std::vector<double> * params, map<int,double>* Rev, double * T,
				   double * V, double * P, int * N, std::map<int,string*> * SI, std::multimap<int,pair<int,int> > *Stoich, 
				   std::vector<set<int> > *Rate, double ATol, double R, set<int>* surf, set<int>* sites, 
				   map<int, pair<double,int> >* InitialFlows, multimap<string, pair<int,int> >* Occupants, 
				   multimap<int, string>* siteBalance, multimap<int,int>* RevRxns,int* FstFlow, int* FstFlowPdct, 
				   double* InletMdot,map<int,int>* MolWeights, vector<KineticsEvaluationFromParam>* kineticEvaluators, 
				   vector<double> PrecondInv)
{
	Reactions = rxns;
	Parameters = params;

	
	Temp = T;
	Volume = V;
	Pressure = P;
	NumOfEquations = N;
	SpeciesIndex = SI;
	StoichInfo = Stoich;
	RateInfo = Rate;
	AbsTolerance = ATol;
	Rg = R;
	RevKinSet =Rev; 
	surfaceSpecies = surf;
	AllSites = sites;
	InitialMolecFlows = InitialFlows;
	SiteOccupants = Occupants;
	SiteBalanceInfo = siteBalance;
	ReverseRxns = RevRxns;
	FirstFlowSpecies = FstFlow;
	FirstFoundProductFlowSpecies = FstFlowPdct;
	InletMassFlowRate = InletMdot;
	SpeciesMolWeights = MolWeights;
	kinetics = kineticEvaluators;
	PreconditionerInverse = PrecondInv;
}


KineticsEvaluationFromParam::KineticsEvaluationFromParam(KineticsEstimates& kEst, vector<double>* Parameters)
{
	Param = Parameters;
	estimates = kEst;
}

double KineticsEvaluationFromParam::getKineticsValue(double Temperature)
{
	double kinetics = 1.0;
	kinetics= kinetics*estimates.kRelative;
	if (estimates.kParam>=0)
		kinetics= kinetics*exp(Param->at(estimates.kParam));
	double Ea = estimates.EaRelative;
	//cout<<"kinetics now is "<<kinetics<<"  "<<estimates.kRelative<<"  "<<estimates.kParam<<"  "<<Param->at(estimates.kParam)<<endl;
	
	if (estimates.EaParam>=0)
		Ea+=exp(Param->at(estimates.EaParam));
	
	kinetics=kinetics*exp(-Ea*1000.0/8.314*(1/Temperature - 1/estimates.kRefTemp));

	kinetics = kinetics*pow(Temperature, estimates.TempExp);

	//if (estimates.kParam==1)
	//	cout<<"kinetics "<<estimates.kRelative<<"  "<<exp(Param->at(estimates.kParam))<<"  "<<estimates.EaParam<<"  "<<estimates.EaRelative<<"  "<<Ea<<"  "<<estimates.TempExp<<"  "<<kinetics<<endl;
	//cout<<"kinetics is "<<kinetics<<endl;
	return kinetics;
}

KineticsEstimates::KineticsEstimates()
{
	kParam = -1; 
	kRelative = 0.0;
	EaParam = -1;
	EaRelative = 0.0;
	kRefTemp = 298.0;
	isReverse = false;
	TempExp = 0.0;

}










