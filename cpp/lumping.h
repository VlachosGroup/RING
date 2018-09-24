#ifndef LUMPING_H
#define LUMPING_H

#include <string>
#include <vector>

#include "molecule.h"

class LumpInfo
{
	protected:
		int size;
		int HydCount;
		int DoubleBondCount;
		int TripleBondCount;
		int AromaticBondCount;
		int MaxRingSize;
		int RingCount;
		bool isCyclic;
		int BranchValue;
		int ringleavesDistance;
		std::string* MoleculeString;
		int rank;//rank of the molecule;
		int leaves;
		int PONAcharacteristic;//Paraffin, Olefin, Naphthenics and aromatics characteristics - To lump all olefins, paraffins and hydrocarbon aromatics of one size together.
			// - 1 don't really care or it is none of P, O, N, or A; 0-it is paraffinic, 1- it is olefinic,  2- cycloalkanes and cycloalkenes; 3-it is hydrocarbon aromatic with alkyl or no substituents; 4- hydrocaron aromatic with alkenyl substitutents 
		int PONAElectronicValue;//stores the Electronic value of POA. 
	public:
		LumpInfo(int,int,bool,int, int, std::string*);
		void setMoleculeString(std::string*);
		std::string getMoleculeString();
		std::string* getMolStringPtr();
		int getSize();
		int getHydrogens();
		int getBranchValue();
		int getringDistanceValue();
		void setBranchValue(int);
		void setringDistanceValue(int);
		int getleaves();
		void setleaves(int);
		void setPONAcharacteristic(int);
		int getPONAcharacteristic();
		int getPONAElectronicValue();
		void setPONAElectronicValue(int);
		void setRank(int);
		int getRank();
		void setBondCounts(int,int,int);//provide double, triple, aromatic
		int getDoubleBondCount();
		int getTripleBondCount();
		int getAromaticBondCount();
		int getMaxRingSize();
		int getRingCount();
		void setMaxRingSize(int);
		void setRingCount(int);
		int getScore();//This returns a score calculated as 3^Doublebonds*5^Triplebonds*7^AromaticBonds
};

class LumpingStrategy
{
	protected:
		bool toLump;
		int chainParameter;//-1 - no lumping of chains, 0-branches closest together, 1-branches farthest apart
		int ringParameter;//-1 - no lumping of rings, 0- substituents closest together, 1-substituents farthest apart
		int paraffinParameter;//-1 no paraffin lumping; 0 - all paraffins lumped to least branched, 1 - all paraffins lumped to most branched
		int olefinParameter;//-1 no olefin lumping; 0 - all olefins lumped to least branched, 1- all olefins lumped to most branched
		int aromaticsParameter;//-1 no aromatics lumping; 0 - all aromatics lumped to least branched, 1- all aromatics lumped to most branched
		int naphthenicsParameter;//-1 no cycloalkanes and cycloalkenes lumping; 1 - all lumped to least branched; 1- all lumped to most branched 
        std::vector<int> MoreLumpingParameter;// -1- all lumped to least branched; 1 - all lumped to most branched
		
		bool HasSetFunctionalLumpingConstraints, HasSetParaffinConstraints, HasSetOlefinConstraints, HasSetNaphthenicsConstraints, HasSetAromaticsConstraints, HasSetMoreLumpingConstraints;
				
		//molecule constraint pointers that describe structural constraints on what molecules can be lumped for initial functional lumping and subsequent lumping of paraffins, olefins, naphthenics, aromatics, and any other additional lumping
		ConstrPtr FunctionalLumpingConstrPtr;
		ConstrPtr ParaffinConstrPtr;
		ConstrPtr OlefinConstrPtr;
		ConstrPtr NaphthenicsConstrPtr;
		ConstrPtr AromaticsConstrPtr;
        std::vector<ConstrPtr> MoreLumpingConstrPtr;

	public:
		LumpingStrategy(bool);
		LumpingStrategy();
		void setParameters(int,int, int, int, int, int);
		bool shoudLump();
		int getchainParameter();
		int getringParameter();
		int getParaffinParameter();
		int getOlefinParameter();
		int getAromaticsParameter();
		int getNaphthenicsParameter();
        std::vector<int> getMoreLumpingParameter();
		void setFunctionalLumpingConstraints(ConstrPtr);
		void setParaffinConstraints(ConstrPtr);
		void setOlefinConstraints(ConstrPtr);
		void setNaphthenicsConstraints(ConstrPtr);
		void setAromaticsConstraints(ConstrPtr);
		void setMoreLumpingConstraints(ConstrPtr,int);//the integer specifies what kind of representative is desired -- see above
		ConstrPtr getParaffinConstraints();
		ConstrPtr getOlefinConstraints();
		ConstrPtr getNaphthenicsConstraints();
		ConstrPtr getAromaticsConstraints(); 
		ConstrPtr getFunctionalLumpingConstraints();
        std::vector<ConstrPtr> getMoreLumpingConstraints();
		bool isThereFunctionalConstraints();
		bool isThereParaffinConstraints();
		bool isThereOlefinConstraints();
		bool isThereNaphthenicsConstraints();
		bool isThereAromaticsConstraints();
		bool isThereMoreLumpingConstraints();

};

class LumpedReaction
{
	protected:
        std::vector<int> reactantsLumpSet;
        std::vector<int> productsLumpSet;
		int ruleIndex;
	public:
		LumpedReaction(std::vector<int>, std::vector<int>, int);
        std::vector<int> getReactantLumps() const; 
        std::vector<int> getProductLumps() const;
		int getRule() const;
		
};

struct LumpedReactionCompare{
  bool operator() (LumpedReaction a, LumpedReaction b) const
  {
	  if (a.getReactantLumps()<b.getReactantLumps()) return true;
	  else if (a.getReactantLumps()>b.getReactantLumps()) return false;
	  else
	  {
		  if (a.getProductLumps()< b.getProductLumps()) return true;
		  else if (a.getProductLumps()>b.getProductLumps()) return false;
		  else return (a.getRule()<b.getRule());
	  }
	  
  }

};

class MoreAdditionalLumpingInfo //Note this class keeps track of all additional lumping stuff the language allows for (that is in addition to PONA)
{
	public:
		int Size;
		int HydCount;
		int DoubleBonds;
		int TripleBonds;
		int AromaticBonds;
		int NumberOfRings;
		int MaxRingSize;
		int WhichMoreLumpingDcl;//keeps the index of which of the possibly many declarations of more lumping
};

struct MoreLumpingCompare{
  bool operator() (MoreAdditionalLumpingInfo a, MoreAdditionalLumpingInfo b) const
  {
	  if (a.WhichMoreLumpingDcl < b.WhichMoreLumpingDcl) return true;
	  else if (a.WhichMoreLumpingDcl > b.WhichMoreLumpingDcl) return false;
	  else
	  {
		  if (a.Size<b.Size)return true;
		  else if (a.Size>b.Size) return false;
		  else
		  {
			 if (a.HydCount<b.HydCount) return true;
			 else if (a.HydCount>b.HydCount) return false;
			 else
			 {
				 if (a.DoubleBonds <b.DoubleBonds) return true;
				 else if (a.DoubleBonds > b.DoubleBonds) return false;
				 else
				 {
					 if (a.TripleBonds < b.TripleBonds) return true;
					 else if (a.TripleBonds > b.TripleBonds) return false;
					 else
					 {
						 if (a.AromaticBonds < b.AromaticBonds) return true;
						 else if (a.AromaticBonds > b.AromaticBonds) return false;
						 else 
						 {
							 if (a.NumberOfRings < b.NumberOfRings) return true;
							 else if (a.NumberOfRings > b.NumberOfRings) return false;
							 else return(a.MaxRingSize < b.MaxRingSize); 
						 }
					 }
				}
			 }
		  }
	  }
  }
};


#endif
