#include <iostream>
#include <string>

#include "common.h"

using std::cout; using std::endl; 
using std::string; //using std::endl; 


//Path implementation
void Path::add(int i)
{
	atom_array.push_back(i);
}

int Path::get(int i) const
{
	return atom_array[i];
}

void Path::reverse()
{
	int size=atom_array.size();
	int n;
	if ((size%2==0))
		n=size/2;
	else
		n=(size-1)/2;
	int temp;
	for (int i=0;i<n;i++)
	{
		temp=atom_array[i];
		atom_array[i]=atom_array[size-1-i];
		atom_array[size-1-i]=temp;
	}
}
void Path::print()
{
	for (unsigned int i=0;i<atom_array.size();i++)
	{
		cout<<atom_array[i]<<"  ";
	}
	cout<<endl;
}

int Path::size() const
{
	return atom_array.size();
}

int Path::getfirst() const
{
	return atom_array[0];
}

int Path::getlast() const
{
	int m=atom_array.size();
	return atom_array[m-1];
}

void Path::remove(int i)
{
	atom_array.erase(atom_array.begin()+i);
}

void Path::join(Path P,int i)
{
	if (getfirst()==i)
		reverse();
	if (P.getlast()==i)
		P.reverse();
							
	for (int k=1;k<P.size();k++)
	{
		add(P.get(k));
	}
}

int Path::intersection(Path P) const
{
	int counter;
	counter=0;
	for (int i=0;i<size();i++)
	{
		for (int j=1;j<P.size()-1;j++)
		{
			if (atom_array[i]==P.get(j))
				counter++;
		}
	}
	return counter;
}

void Path::clear()
{
	atom_array.clear();
}

bool Path::contains(int i) const
{
	int flag=0;
	for (unsigned int k=0;k<atom_array.size();k++)
	{
		if (atom_array[k]==i)
		{
			flag=1;
			break;
		}
	}
	if (flag==1) return true;
	else return false;
}

//end of Path implementation

//Ringset implementation
void Ringset::add(Path P)
{
	ring_list.push_back(P);
}

Path Ringset::get(int i) const
{
	return ring_list[i];
}

void Ringset::remove(int i)
{
	ring_list.erase(ring_list.begin()+i);
}

void Ringset::print(int i)
{
	for (int j=0;j<ring_list[i].size();j++)
	{
		cout<<ring_list[i].get(j)<<"  ";
	}
	cout<<endl;
}

void Ringset::print_all()
{
	for ( int i=0;i<ring_list.size();i++)
	{
		print(i);
	}
}
int Ringset::size() const
{
	return ring_list.size();
}

void Ringset::sizesort()
{
	Path temp;
	for (unsigned int i=1;i<ring_list.size();i++)
	{
		for (unsigned int j=0;j<ring_list.size()-i;j++)
		{
			if (ring_list[j].size()<ring_list[j+1].size())
			{
				temp=ring_list[j];
				ring_list[j]=ring_list[j+1];
				ring_list[j+1]=temp;
			}
		}
	}
	
}

bool Ringset::is_present(int a) const
{
	int flag=0;
	for (unsigned int i=0;i<ring_list.size();i++)
	{
		for (unsigned int j=0;j<ring_list[i].size();j++)
		{
			if (ring_list[i].get(j)==a)
			{
				flag=1;
				break;
				
			}
		}
	}
	if (flag==0)
		return false;
	else
		return true;
}

int Ringset::CountAtom(int a) const
{
	int counter = 0;
	for (int i=0;i<ring_list.size();i++)
	{
		for (int j=0;j<ring_list[i].size();j++)
		{
			if (ring_list[i].get(j)==a)
				counter++;
		}
	}

	return counter;
}

int Ringset::SizeContainingAtom(int a) const
{
	int flag =-1;
	for (int i=0;i<ring_list.size();i++)
	{
		for (int j=0;j<ring_list[i].size();j++)
		{
			if (ring_list[i].get(j)==a)
			{
				
				flag=i;
				break;
			}
		}
	}
	if (flag>=0)
		return ring_list[flag].size();
	else return 0;
}
//end of Ringset implementation

//start env_tuple implementation

void env_tuple::set_flag(int i)
{
	flag=i;
}
int env_tuple::get_flag()
{
	return flag;
}
void env_tuple::set_type(string s)
{
	type=s;
}
string env_tuple::get_type()
{
	return type;
}

int env_tuple::get_freq()
{
	return freq;
}

void env_tuple::set_freq(int i)
{
	freq=i;
}

int env_tuple::get_comp_opr()
{
	return comp_opr;
}

void env_tuple::set_comp_opr(int i)
{
	comp_opr=i;
}

//end of env_tuple implementation

//start of env_set implementation

void env_set::add(env_tuple et)
{
	environment.push_back(et);
}
env_tuple env_set::get(int i)
{
	return environment[i];
}

void env_set::remove(int i)
{
	environment.erase(environment.begin()+i);
}

bool env_set::set_is_empty()
{
	if (environment.empty())
		return true;
	else
		return false;
}
int env_set::size()
{
	return environment.size();
}
//end of implementation of env_set

