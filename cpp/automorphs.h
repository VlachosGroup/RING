#ifndef AUTOMORPHS_H
#define AUTOMORPHS_H

#include <vector>
#include <set>
#include <map>

#include "molecule.h"

class Automorphs
{
	protected:
		const Molecule* mol;
		int NumAuto;
		int Symmetry;
		int RootClass;
		int RootAtom;
		
        std::map<int,std::vector<int> > ClassRankAtomMap;
        std::map<int,std::set<int> > AtomPredecessorMap;
		std::vector<int> DistanceFromRoot;
		std::vector< std::set<int> >  AllPermutationClasses;
        std::map<int, std::vector<std::set<int> > > NeighboringClassSetsForDepthOne;// for each atom 'i' at depth 1 from root, this map stores the sets of symmetric classes of atoms that are in the atom i's immediate neighbors
		void MapClassRanksAndAtoms();
		void SetRootClassAndAtom();
		void BFSFromRoot();
		void FindAllPermutationClasses();
		int GetSymmetryNumberForAllLeaves();
		int SymmetryCorrections();
		int GetSymmetryNumberFromEqSet(std::set<int>, int);
		int GetExtDbBondCorrection();
		void traverseDoubleBonds(int, std::vector<int>*);
		void CalculateAutomorphsAndSymmetryNumbers();
		
	public:
		Automorphs(const Molecule&);
		int NumberOfAutomorphs();
		int SymmetryNumber();
};

#endif
