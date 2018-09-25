# distutils: language = c++

from cpython.ref cimport PyObject
from reactiontype cimport ConstrPtr
from reactiontype cimport CombinedConstrPtr
from reactiontype cimport C_Molecule
from reactiontype cimport C_Patternmatch
from reactiontype cimport C_Substructure
from reactiontype cimport C_ReactionType
from reactiontype cimport c_patternsize
from reactiontype cimport c_SiteType
from lumping cimport C_LumpInfo
from lumping cimport C_LumpingStrategy 


from reactiontype cimport C_RxnNetGen
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool

#cdef vector[C_Molecule] c_GlobalConstraints
GlobalConstraints = []
ReactionReactantConstraints = []
ReactionProductConstraints = []
ReactionCombinedConstraints = []

#cdef public struct PyArrayObject:
#    C_Molecule mol


def patternsize(s1):
    if type(s1) is unicode:
        cpp_str = (<unicode>s1).encode('utf8')
    rval = c_patternsize(cpp_str)
    return rval

cdef class LumpingStrategy:
    cdef C_LumpingStrategy* _ptr
    
    def __cinit__(self, x):
        self._ptr = new C_LumpingStrategy(x)
        
    def should_lump(self):
        return self._ptr[0].shoudLump()


cdef class Molecule:
    cdef C_Molecule* _ptr
    
    def __cinit__(self, x, y):
        if type(x) is unicode:
            x = (<unicode>x).encode('utf8')
        y = <int>y
        self._ptr = new C_Molecule(x,y)
        
    def __init__(self, x, y):
        pass
    
    def printmol(self):
        print("in molecule")
    
cdef class Substructure:
    cdef C_Substructure* _ptr
    
    def __cinit__(self, x, y):
        if type(x) is unicode:
            x = (<unicode>x).encode('utf8')
        self._ptr = new C_Substructure(x,y)
    
    def __init__(self, x, y):
        pass

cdef class Patternmatch:
    cdef C_Patternmatch* _ptr
#    
    def __cinit__(self, molecule, substructure, constraint_index):
        #mol = Molecule("CCC", 3)
        #mol = Molecule()
        c_mol = (<Molecule?>molecule)._ptr
        c_substruct = (<Substructure?>substructure)._ptr
        self._ptr = new C_Patternmatch(c_mol[0], c_substruct[0], 
                                       constraint_index)
    
    def __init__(self, molecule, substructure, constraint_index):
        pass
#    
#    
#
    def __cinit__(self):
        self._ptr = new C_Patternmatch()
#        
    def get_distinct_matches(self):
        return self._ptr[0].GetDistinctMatches()

cdef class ReactionType:
    cdef C_ReactionType* _ptr

    def __cinit__(self):
        self._ptr = new C_ReactionType()
        
#    def __init__(self, constraints = 0):
#        self.RxnConstraints = []
    
    def add_reactant_pattern(self, substructure):
        self._ptr[0].add_reactant_pattern(
                (<Substructure?>substructure)._ptr[0])
        
    def AddReactantCopy(self, x, reactant_duplicate):
        self._ptr[0].AddReactantCopy(x, (<map[int,int]?>reactant_duplicate))
        
    def add_reactant_constraints(self, reactant_const):
        self.RxnConstraints = reactant_const

    def disconnect_bond(self, x, y):
        self._ptr[0].disconnect_bond(x, y)
        
    def connect_bond(self, x, y):
        self._ptr[0].connect_bond(x, y)
        
    def increaseBO(self, x, y, z):
        self._ptr[0].increaseBO(x, y, z)
        
    def decreaseBO(self, x, y, z):
        self._ptr[0].decreaseBO(x, y, z)

    def set_rule_name(self, x):
        if type(x) is unicode:
            x = (<unicode>x).encode('utf8')
        self._ptr[0].setRuleName(x)
         
cdef class RxnNetGen:
    cdef C_RxnNetGen* _ptr
    
    def __cinit__(self):
        self._ptr = new C_RxnNetGen()
    
    def add_initial_reactants(self, reactant_list):
        c_reactant_list = []
        for i in reactant_list:
            if type(i) is unicode:
                i = (<unicode>i).encode('utf8')
                c_reactant_list.append(i)
        self._ptr[0].AddInitialReactants(c_reactant_list)
        
    def add_reaction_rules(self, reaction_types):
        cdef vector[C_ReactionType] c_rxn_types
        for rxn_type in reaction_types:
            c_rxn_type = (<ReactionType?>rxn_type)._ptr[0]
            c_rxn_types.push_back(c_rxn_type)
        self._ptr[0].AddReactionRules(c_rxn_types)
    
    def GenerateNetwork(self):
        self._ptr[0].GenerateNetwork()
        
    def set_all_lumped_reactions(self):
        self._ptr[0].SetAllLumpedRxns()
        
    def print_rxn_list(self):
        self._ptr[0].print_rxnlist()
        
    def set_calc_thermo(self, x):
        self._ptr[0].SetCalcThermo(x)
        
    def add_lumping_strategy(self, x):
        self._ptr[0].AddLumpingStrategy(<C_LumpingStrategy?>x)
        
    def add_composite_sites(self, x):
        self._ptr[0].AddCompositeSites(<vector[pair[string, c_SiteType] ]?>x)
        
    def check_global_constraints(self, constraint):
        global_constraint = constraint
        
cdef public bool check_reactant_constraint0(C_Molecule mol, int x, int y):
#mol is the molecule to check for
#x is the reactiontype in list of reaction rules
#y is going to be 0    
    print("In constraint check!")
    print(x)
    #C_ReactionType = Rtlist[x]
    if x == 6:
        print("found rule 6")
        print(y)
        substruct = Substructure("O1(-C2)", c_patternsize("O1(-C2)"))
        c_substruct = (<Substructure?>substruct)._ptr
        #c_mol = (<Molecule?>mol)._ptr
        pattern = new C_Patternmatch(mol, c_substruct[0], 1)
        value = pattern.GetDistinctMatches()
        #inpu = (<object?>mol)._ptr[0]
        print(value)
#        for f in GlobalConstraints:
#            print("In here")
#            print(type(f))
        if value == 0:
            return True
        else:
            return False
    elif x == 7:
        print("found rule 7")
        print(y)
        substruct = Substructure("O1(-C2)", c_patternsize("O1(-C2)"))
        c_substruct = (<Substructure?>substruct)._ptr
        #c_mol = (<Molecule?>mol)._ptr
        pattern = new C_Patternmatch(mol, c_substruct[0], 1)
        value = pattern.GetDistinctMatches()
        print(value)
        if value == 0:
            return True
        else:
            return False
    else:
        print(" ")
        return True

cdef public bool check_reactant_constraint1(C_Molecule mol, int x, int y):
    return True
    
cdef public bool check_combined_constraint(C_Molecule mol0, C_Molecule mol1, 
                                           int x):
    return True

cdef public bool check_combined_constraint_for_reactiontype(
        C_ReactionType reacttype, C_Molecule mol0, C_Molecule mol1):
    return True

cdef public int check_combined_constraints_for_Rtlist(int x):
    return 1

cdef public bool check_product_constraint(C_Molecule mol, int x):
    return True
        
cdef public struct Bunny: # public type declaration
    int vorpalness

cdef public int spam # public variable declaration

cdef public void grail(Bunny *b): # public function declaration
    print("Hi")
    a=b.vorpalness*2
    print(a)
    
cdef public bool check_global_constraints(C_Molecule mol):
    return True    


