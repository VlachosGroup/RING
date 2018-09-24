#ifndef PATTERNMATCH_H
#define PATTERNMATCH_H

#include <string>
#include <vector>
#include <set>
#include <utility>

#include "molecule.h"
#include "substructure.h"

class Patternmatch //Patternmatch class 
{
	protected:
		std::vector < std::vector<int> > Matches;//stores all the matches
		std::vector < std::vector<int> > rank_matches;//stores the rank of each of the atoms of all the matches
		std::vector < std::vector <int> > M;//the Matrix M
		std::vector <std::vector <std::vector <int> > > Ma;//the 3D matrix that stores all intermediate M Matrices
		void find_matches(const Molecule&,const Substructure&, int);//function that finds the matches
		void refine(const Molecule&,const Substructure&);//function for refinement
		int check_zeros();//function that checks if M has an entire row of zeros
		bool checkatomtypematch(std::string, std::string);
		bool checkconnectivitymatch(int, int);
		bool checkatomenvironmentmatch(int, int);
		bool checkatomflagmatch(int,int);
		bool IsSameCharacteristrics(int, int);
		int Subsize;//size of Substructure 
		int Molsize;//Size of the Molecule
		const Molecule* Mo;
		const Substructure* Su;
		std::vector<int> Hfactor;
		bool atomsetmatch(int,int);
		std::vector< int >unique_matcheslist;//list of indices of unique matches
		std::vector<int> unique_matches_frequency;//frequency of occurences of each
		void InitializeAndSetMa(int);
		void calcHfactor();
		Patternmatch(const Molecule&, const Substructure&,int, int);//constructor - the molecule, the substructure to be found, option - 0 implies find all matches, anything else implies stop at first match. The fourth option is an atomindex of the molecule to which the first atom of the substructure wil be matched.
		//This feature is used specifically for nested SMARTS (atom environment constraints being the presence or absence of entire groups). The index is set to -1 for no specification and a valid index value at other points. 

	public:
		Patternmatch(const Molecule&, const Substructure&,int);//Constructor - the molecule, the substructure to be found, and an option - 0 implies find all matches, anything else implies stop at 1. 
		Patternmatch();
		void list_matches();//lists all the matches in an ordered manner
		void unique_matches();//lists the unique matches - gives matches which are topologically unique - but different unique matches could have the same set of atoms. For example, matching CC in ethane would give two unique matches. To be used for reactions
		int GetDistinctMatches();//lists the distinct matches - here, any two pair of matches of the same pattern should have atleast one distinct atom in each match. To be used for structural constraints
		void print_M();//prints M
		int H_factor(int);//gets the H_factor for a match that has to be multiplied with MatchFrequency to get overall reaction frequency
		int H_factor_unique_match(int);//gets the H_factor for a particular unique match, utility same as above
		void print_Ma(int);//prints out Md which is Ma[d]
		int number_of_matches();//returns the number of matches
		int number_of_unique_matches();//returns the number of unique matches
		std::vector<int>get_unique_matches(int);//gets the particular unique match
		std::vector<int> get_match(int);//get a particular match
		int getMatchFrequency(int);//gets the number of times ith unique match (i is the argument) occurs in the molecule
        std::set<int> AtomsCoveredBySubstr(int);//gets the set of atoms matched by the specified atom of hte substructure from the list of all matches
		bool IsMatchWithinRing(int);//checks if a particular match is completely contained in a ring
        std::pair<int,int> GetDistinctAllAndRingMatches();//gives a pair of all distinct matches and those that are only within a ring
		friend class Molecule;
		//friend class Atomcontainer;
};	

#endif
