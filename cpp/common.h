#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <string>
#include <utility>

class Path//Path is a collection of connected atoms
{
	protected:
        std::vector < int > atom_array;//array of atoms
	public:
		void add(int);//add an element
		int get(int) const;//get an element
		void reverse();//reverse the contents of the array
		void print();//print the array
		int size() const;//size 
		int getfirst() const;//the first element
		int getlast() const;//the last element
		void remove(int);//remove an element
		void join(Path, int);//join two paths
		int intersection(Path) const;//taking the intersection of two sets of atoms constituting the two paths
		void clear();//remove the contents of the array
		bool contains(int) const;// checks if a particular element is in the Path
		
};

class env_tuple//triplet containing info on atom environment
{
	protected:
		
		int flag;//a flag to indicate presence or absence - 0 is absent and 1 is present
        std::string type;//string that defines the atomtype (in case of specifications involving single atoms) or SMARTS string (in case of environment constraints being a group)
		int freq;//frequency or the number of occurences
		int comp_opr;//comparison operator (0 for equality, 1 for greater than, -1 for less than)
		
	public:
		void set_flag(int);//set the flag value
		int get_flag();//get the flag value
		void set_type(std::string);//set the atomtype
        std::string get_type();//get the atomtype
		int get_freq();
		void set_freq(int);
		void set_comp_opr(int);
		int get_comp_opr();
		
};

class env_set//a set of environment triplets
{
	protected:
        std::vector<env_tuple> environment;//array of triplets
	public:
		void add(env_tuple);//add a triplet
		env_tuple get(int);//get a triplet
		void remove(int);//remove a triplet
		bool set_is_empty();//check if the array is empty
		int size();//number of triplets - the size of environment array
};

class IntPair//IntPair is a class to store pair of elements
{
	public:
		int first;
		int second;
};

class Triplet//Triplet is a class to store element triplets
{
	public:
		int first;
		int second;
		int third;
};

class Ringset//collection of rings
{
	protected:
        std::vector < Path > ring_list;//array of paths
	public:
		void add(Path);//add a path
		Path get(int) const;//get a path
		void print(int);//print particular ring path
		void print_all();//print all rings
		void remove(int);//remove a path
		int size() const;//number of paths
		void sizesort();//sort by the size of the paths from largest to smallest
		bool is_present(int) const;
		int SizeContainingAtom(int) const;
		int CountAtom(int) const;
};

class bond_operations//contains info on the bond changes
{
	public:
		int first_atom_label;//first of the two atoms describing hte bond
		int second_atom_label;//second of the two atoms describing the bond
		int set_bond_order;//flag to set the bond order: -1 implies inactive, 0 is disconnect (break) and 1 is connect (form bond)
		int change_bond_order;//increase or decrease bond order implied by -1 and 1 for reducing the BO by 1/ increasing it by 1. 
							  // A value of 0 would denote no change in bond order. Note that the set_bond_order and change_bond_order cannot hold independent values
							  //Note that setting the bond order to 1 and change_bond_order to 1 is invalid, for example. 
};

struct classcomp {
  bool operator() (std::string* a, std::string* b) const
  {return (*a)<(*b);}
};

enum SiteType {  Homogeneous, Heterogeneous};


typedef std::pair<unsigned int, unsigned int> UnsignedIntPair;
typedef std::pair<std::pair<int,int>,int> PONALumpingInfo; //the inner pair holds size and hydrogen count, while the outer pair's second value is PONAElectronic value



#endif
