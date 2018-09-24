#ifndef GENERATED_RXN_H
#define GENERATED_RXN_H

//#include <iostream>
//#include <fstream>
#include <string>
//#include <cstring>
//#include <sstream>
//#include <std::vector>
//#include <list>
//#include <set>
#include <map>
//#include <queue>
#include <deque>
//#include <algorithm>
//using namespace std;
//
//
class generated_rxn
{
	protected:
		std::deque<std::string*> reactants;
		std::deque<std::string*> products;
		int frequency;
		std::string reactionMF;
		int RuleIndex;
		bool IsRxnIntramolecular;
		std::map<int, int> ParentMolecule;
		std::map<int, std::vector<int> >AtomsFromOtherReactants;//keeps track of which product has come from second and subsequent reactants
		int countSpecies(std::deque<std::string*>&, std::string*);//counts the number of species in the reactant/product std::deque
	public:
		generated_rxn();
		void add_reactants(std::string*);
		void add_products(std::string*);
		std::string* get_reactants(int);
		std::string* get_products(int);
		std::string reactionstring();
		bool IsReactant(std::string*);
		bool IsProduct (std::string*);
		int number_pdcts();
		int number_reactants();
		int get_frequency();
		void set_frequency(int);
		void set_reactionMF(std::string);
		std::string get_reactionMF();
		void set_rule(int);
		int get_rule();
		int occurence(std::string*); //gets the absolute value of the occurrence -> abs(stoich in reactant - stoich in products)
		void setAtomsFromOthers(int, std::vector<int> );
		void setParentMolecule(std::map<int, int>);
		std::vector<int> getAtomsFromOthers(int);
		int getParentMolecule(int);
		int getParentMolecule(std::string*);
		bool isSameNetReaction(generated_rxn&);
		bool isReverseReaction(generated_rxn&);
		bool isReactionIntramolecular();
		void setIntramolecularity(bool);
		std::vector<std::string*> getDaughterMolecules(std::string*);//gets all the molecules in the product having the given molecule as closer (parent).
		int getReactantStoich(std::string*);//gets stoich in reactant 
		int getProductStoich(std::string*);//gets stoich in product
		int NetOccurence(std::string*);//gives the actual value - net negative means more reactants - net positive means more products. 
		int getdeltaN();//gives the change in the number of moles. (reactants - products)
		int NetPatternDiff(std::string); //gives the net difference of the pattern occurrences in the reactants and products - negative means more in reactants! 
		int netMassDiff();//calculates if there is a mass difference between reactants and products

		bool ProductsWithSameRadicalType(); //checks if there are two radicals in the product of the same atomtype
		int IdenticalProductsFactor();//calculates the product of number of different products that are identical to each other. e.g. if products are P, P,and P - the factor is 3; if P P Q Q - factor is 4, but if P P Q, then it's still 1. 
		
};


#endif
