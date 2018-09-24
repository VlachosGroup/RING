#include <map>
#include <vector>
#include <utility>
#include <set>
#include <cmath>

#include "additionalfunc.h"
#include "automorphs.h"

using std::vector; using std::pair; using std::map; using std::multimap;
using std::set; 
using std::pow;
/*
1. find the automorphs normally
2. Calculate the internal rotation of leaves (include double bond case as well)
3. Corrections
   1. an atom whose neighbors can be distributed into two sets of symmetry class - reduce by a factor of 2
   2. spiro - reduce by a factor of 2 or retain as appropriate
   3. go to all atoms at distance d=1, and find if any of them has children of the same class as the root class
      figure out if including this root atom, one can still account for symmetry - if not reduce by a factor equal to size of the rootclass
   4. find extended C=C, reduce by an appropriate factor as required
*/
	

Automorphs::Automorphs(const Molecule& molec)
{
	mol = &molec;

	NumAuto = 1;
	Symmetry = 1;
	RootClass = -1;
	RootAtom  = -1;

	//first map the class ranks to their atoms
	MapClassRanksAndAtoms();

	//now gloss over this map to see if there is a class with 2 or more members. Take the first one as root class. 
	SetRootClassAndAtom();
	
	
	//Now begin breadth first search from RootAtom to find the distance of each atom from the root. 
	BFSFromRoot();
	//Set the NeighboringClassSetsForDepthOne with dummy vectors
	vector<set<int> > Dummy;
	for (int i=0; i<mol->getNN(RootAtom);i++)
		NeighboringClassSetsForDepthOne.insert(pair<int, vector<set<int> > >(mol->get_adjacency(RootAtom,i), Dummy));
	
	//now find all the sets of permutation classes
	FindAllPermutationClasses();
	//finally, find the num of automorphs as well as symmetry number
	
	CalculateAutomorphsAndSymmetryNumbers();
	

}

void Automorphs::MapClassRanksAndAtoms()
{
	
	for (int i=0;i<mol->getsize();i++)
	{
		if (ClassRankAtomMap.count(mol->getclassRank(i))!=0)
			ClassRankAtomMap[mol->getclassRank(i)].push_back(i);
		else 
		{
			vector<int> v;
			v.push_back(i);
			ClassRankAtomMap.insert(pair<int, vector<int> >(mol->getclassRank(i),v));

		}
	}
}

void Automorphs::SetRootClassAndAtom()
{
	map<int,vector<int> >::iterator it;

	for (it=ClassRankAtomMap.begin();it!=ClassRankAtomMap.end();it++)
	{
		if (it->second.size()>=2)
		{
			RootClass = it->first;
			break;
		}
	}
	if (RootClass ==-1)RootClass = ClassRankAtomMap.begin()->first;

	RootAtom = ClassRankAtomMap[RootClass].at(0);//root atom is the first atom in the vector of atoms corresponding to the RootClass
	//cout<<"root is "<<RootAtom<<endl;
	
}

void Automorphs::BFSFromRoot()
{
	set<int> s;

	for (int i =0;i<mol->getsize();i++)//initializing
	{
		AtomPredecessorMap.insert(pair<int, set<int> >(i, s));
		DistanceFromRoot.push_back(100000);//a large number
		
	}
	// finding the distance of each atom from the root! 	
	vector<int> NextAtom;
	NextAtom.push_back(RootAtom);
	DistanceFromRoot[RootAtom] = 0;
	while (!NextAtom.empty())
	{
		int CurrentAtom = NextAtom[0];
		NextAtom.erase(NextAtom.begin());
		for(int i=0;i<mol->getNN(CurrentAtom);i++)
		{
			int Neighbor = mol->get_adjacency(CurrentAtom,i);

			if (DistanceFromRoot[Neighbor]>=  DistanceFromRoot[CurrentAtom])//if neighbor is not actually closer to the root than the current one
			{
				AtomPredecessorMap[Neighbor].insert(CurrentAtom);//update the CurrentAtom as predecssor of Neighbor
				NextAtom.push_back(Neighbor);
				DistanceFromRoot[Neighbor] = DistanceFromRoot[CurrentAtom]+1;
			}
			
		}
	}
}


void Automorphs::FindAllPermutationClasses()
{

	map<int,vector<int> >::iterator it;

	
	for (it = ClassRankAtomMap.begin();it!=ClassRankAtomMap.end();it++)
	{
		//for each class, arrange the elements according to their depth.
		multimap<int, int> DepthAndPermutClassAtomsMap;
		for (unsigned int i=0;i<it->second.size();i++)
		{
			int CurrentAtom = it->second.at(i);
			DepthAndPermutClassAtomsMap.insert(pair<int,int> (DistanceFromRoot[CurrentAtom],CurrentAtom));
		}

		//now create the set of permutation classes
		multimap<int, int>::iterator iter;
		int prevDistance = -1;
		set<int> NewClass;
		for (iter = DepthAndPermutClassAtomsMap.begin();iter != DepthAndPermutClassAtomsMap.end();iter++)
		{
			if (iter->first!=prevDistance)
			{
				AllPermutationClasses.push_back(NewClass);
			}
			//NewClass.insert(iter->second);
			AllPermutationClasses.back().insert(iter->second);
			prevDistance = iter->first;
		}

	}

}

void Automorphs::CalculateAutomorphsAndSymmetryNumbers()
{
	map<int,vector<int> >::iterator it;
	for (unsigned int i =0 ;i<AllPermutationClasses.size();i++)
	{
		set<int> Class = AllPermutationClasses[i];
		
		map<int, set<int> > ExclusivePredecessorsMap;
		map<int, set<int> > ExclusivePredessorsBondMap;//keeps track of bonds connecting the atom and their predecessors
		
		set<int>::iterator set_it;
		for (set_it = Class.begin();set_it!=Class.end();set_it++)
		{
			
			//Predecessors.insert(AtomPredecessorMap[*set_it].begin(), AtomPredecessorMap[*set_it].end());
			
			if (AtomPredecessorMap[*set_it].size()==1)
			{
				int Pred=*AtomPredecessorMap[*set_it].begin();
				int BondOrder = mol->find_BO(Pred, *set_it);
				if (ExclusivePredecessorsMap.count(Pred)>=1)//if Pred is already in the map
				{
					ExclusivePredecessorsMap[Pred].insert(*set_it); //just insert the new one
					ExclusivePredessorsBondMap[Pred].insert(BondOrder);
				}
				else 
				{
					set<int> temp;
					temp.insert(*set_it);
					ExclusivePredecessorsMap.insert(pair<int, set<int> >(Pred,temp));
					temp.clear();
					temp.insert(BondOrder);
					ExclusivePredessorsBondMap.insert(pair<int, set<int> >(Pred,temp));						
				}
			}
		}
		map<int, set<int> >::iterator map_it;
		
		for (map_it=ExclusivePredecessorsMap.begin();map_it!=ExclusivePredecessorsMap.end();map_it++)
		{
			//one more check needs to be done to see if all the successors are connected with the same bond order
			if (ExclusivePredessorsBondMap[map_it->first].size()==1 && Class.size()>1)
			{
				NumAuto*=factorial(map_it->second.size());
				Symmetry*= GetSymmetryNumberFromEqSet(map_it->second, map_it->first);
			}
			
			// if the predecessor atom is actually at depth =1 (connected directly to the root atom)
			//then we need to consider it for corrections to calculate symmetry; so keep track of it separately
			if (DistanceFromRoot.at(map_it->first)==1)
				NeighboringClassSetsForDepthOne[map_it->first].push_back(map_it->second);
		}
		
	}
	NumAuto*= ClassRankAtomMap[RootClass].size();
	Symmetry*= ClassRankAtomMap[RootClass].size();

	
	Symmetry = (Symmetry*GetSymmetryNumberForAllLeaves())/(SymmetryCorrections());
	
}

int Automorphs::NumberOfAutomorphs()
{
	return NumAuto;
}

int Automorphs::SymmetryNumber()
{
	return Symmetry;
}

int Automorphs::GetSymmetryNumberForAllLeaves()
{
	int IntSym =1;
	for (int i=0;i<mol->getsize();i++)
	{
		int hyd = mol->getHydrogens(i);
		if (mol->getNN(i)<=1 && hyd>=1)//if the atom is a leaf with nonzero hydrogens, or if all the atoms around it are hydrogens (this is also technically a leaf and is not taken care of by the parent algorithm
			IntSym*=hyd;
		if (mol->getNN(i)==0 && mol->getatom(i)->get_up() + mol->getatom(i)->get_lp()==0 && mol->getatomtype(i)!="H" && hyd>1) // this check is to ensure that the unpaired electron and lone pair is treated as a possible bond and is differentiated from the rest of the hydrogen 
			IntSym*= (hyd-1);
		
	}
	
	return IntSym;
}

int Automorphs::SymmetryCorrections()
{
	/*
    Corrections
    1. an atom whose neighbors can be distributed into two sets of symmetry class - reduce by a factor of 2
    2. spiro - reduce by a factor of 2 or retain as appropriate (if two sets of symmetry classes, then reduce, else leave as it is)
    3. go to all atoms at distance d=1, and find if any of them has children of the same class as the root class
       figure out if including this root atom, one can still account for symmetry - if not reduce by a factor equal to size of the rootclass
    4. find extended C=C, reduce by an appropriate factor as required

	*/

	int CorrectionFactor = 1;

	//first and second correction - both of them basically depend on the same criterion 
	for (int i=0;i<mol->getsize();i++)
	{
		//trying to figure out if there are two or more atoms of each of the classes. 
		map<int,int> UniqueClassesFreqMap = mol->getClassesFreqMapNeighboringAtom(i);
		map<int, int>::iterator it;
		bool NoClassWithOneAtom = true;

		if (UniqueClassesFreqMap.size()>=2)
		{
			for (it=UniqueClassesFreqMap.begin();it!=UniqueClassesFreqMap.end();it++)
			{
				if (it->second<2)
				{
					NoClassWithOneAtom = false;
					break;
				}
			}
			//if indeed no class exists with only one atom and has no hydrogens too, then we need to divide a factor of 2, n times, 
			//where n = one less than the number of classes because the symmetry had been over counted 2^n times! 
			if (NoClassWithOneAtom && mol->getHydrogens(i)==0)
				CorrectionFactor*= (int) pow(double(2),(int) UniqueClassesFreqMap.size()-1);// I am guessing n will be more than 1 only if the atom is S or P or the likes (with > 4 valency)
		}
			
	}
	//third correction

	map<int, vector< set<int> >  >::iterator it;

	//for each atom at depth = 1 from the root
	for (it = NeighboringClassSetsForDepthOne.begin();it!=NeighboringClassSetsForDepthOne.end();it++)
	{
		
		//for each class of atoms neighboring this atom
		for (unsigned int i = 0; i<it->second.size();i++)
		{
			
			//check if the class of this set is the same as root class
			int FirstAtomOfSet = *(it->second.at(i).begin());
			int ClassOfSet = mol->getclassRank(FirstAtomOfSet);
			int BondOrderSet = mol->find_BO(FirstAtomOfSet,it->first);
			if (ClassOfSet == RootClass && mol->find_BO(RootAtom,it->first)==BondOrderSet)//root is of the same class and bonded with the same bond order as the set
			{
				
				set<int> EqClass;
				EqClass.insert(it->second.at(i).begin(),it->second.at(i).end());
				//EqClass.insert(ClassOfSet);
				EqClass.insert(RootAtom);//TODO is this correct or is the prev line?
				if (GetSymmetryNumberFromEqSet(EqClass, it->first)!=EqClass.size())
					CorrectionFactor*=ClassRankAtomMap[RootClass].size();
			}
		}

	}
	//fourth correction
	CorrectionFactor*= GetExtDbBondCorrection();
	
	return CorrectionFactor;

}

int Automorphs::GetSymmetryNumberFromEqSet(set<int> s, int p)
{
	set<int> UniqueOtherAtomsClassSet; // keeps track of the distinct classes of the other neighboring atoms outside of s
	bool isOtherBondHeavy = false;
	for (int i =0;i<mol->getNN(p);i++)
	{
		int neighbor = mol->get_adjacency(p,i);
		
		if (s.count(neighbor)==0) 
		{
			UniqueOtherAtomsClassSet.insert(mol->getclassRank(neighbor));
			//if (!isOtherBondHeavy && BO[p][i]>1 && BO[p][i]<=4)isOtherBondHeavy = true;
		}
		
	}
	if (mol->getatom(p)->get_up() >= 1 ||mol->getatom(p)->get_lp()== 1 )UniqueOtherAtomsClassSet.insert(-1);//inserting a dummy value - to make sure that the unpaired electron or the lone pair is treated as a separate class altogether! 
	if ( (UniqueOtherAtomsClassSet.size()==0) || (UniqueOtherAtomsClassSet.size()==1 && mol->getHydrogens(p)==0))// 0 implies all others are hydrogens or there are none, 1 implies there's just one other type of atoms
	{
		return s.size();
	}
	else return 1;


}

void Automorphs::traverseDoubleBonds(int atom, vector<int>* path)
{
	int totalHeavyBonds = mol->dbcount(atom) + mol->tpcount(atom);
	int parent = -1;
	if (path->size()>0) parent= path->back();
	
	if ((totalHeavyBonds==1 && parent ==-1) ||(totalHeavyBonds>=2 && parent>=0))
	{
		for (int i = 0;i<mol->getNN(atom);i++)//for each neighbor of atom
		{
			if (mol->get_BO(atom,i)==2 || mol->get_BO(atom,i)==3)
			{
				int neighbor = mol->get_adjacency(atom,i);
				if (neighbor!=parent)
				{
					
					path->push_back(atom);
					traverseDoubleBonds(neighbor,path);
					break;
				}
			}
		}
	}
	else
	{
		path->push_back(atom);
	}
}
	


		

int Automorphs::GetExtDbBondCorrection()
{
	int correction = 1;
	set<int> AtomsCovered;

	for (int i=0;i<mol->getsize();i++)
	{
		if (AtomsCovered.count(i)==0)
		{
			
			if (mol->dbcount(i)==1)
			{
				map<int,int> UniqueClassesFreqMap; //the number of times each class occurs in the neighbors
				//UniqueClassesFreqMap = mol->getClassesFreqMapNeighboringAtom(i);

				for (int j=0;j<mol->getNN(i);j++)
				{
					if (mol->get_BO(i,j)!=2)
					{
						int classOfNeighbor = mol->getclassRank(mol->get_adjacency(i,j));
						if (UniqueClassesFreqMap.count(classOfNeighbor)==0)
							UniqueClassesFreqMap.insert(pair<int,int>(classOfNeighbor,1));
						else UniqueClassesFreqMap[classOfNeighbor]++;
					}
				}

				int num_hydrogens = mol->getHydrogens(i);
				
				if ( (UniqueClassesFreqMap.size()==1 && num_hydrogens==0) || (UniqueClassesFreqMap.size()==0 && num_hydrogens>0) )
				{
					
					vector<int> pathTraversed;
					traverseDoubleBonds(i,&pathTraversed);
					
					AtomsCovered.insert(pathTraversed.begin(), pathTraversed.end());
					int endAtom = pathTraversed.back();
					

					if (!(mol->getNN(endAtom)==1 && mol->getHydrogens(endAtom)<=1))//a correction factor is not required if the endAtom is a leaf or has just one hydrogen emanating and one other (parent) atom			
						correction*=2;
				}
			}
		}
	}

	
	return correction;

}



