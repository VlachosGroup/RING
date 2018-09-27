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

def patternsize(molstring):
    if type(molstring) is unicode:
        cpp_str = (<unicode>molstring).encode('utf8')
    rval = c_patternsize(cpp_str)
    return rval

cdef class LumpingStrategy:
    cdef C_LumpingStrategy* _ptr
    
    def __cinit__(self, lumpconstraints):
        self._ptr = new C_LumpingStrategy(lumpconstraints)
        
    def should_lump(self):
        return self._ptr[0].shoudLump()


cdef class Molecule:
    cdef C_Molecule* _ptr
    
    def __cinit__(self, molstring, stringsize):
        if type(molstring) is unicode:
            molstring = (<unicode>molstring).encode('utf8')
        stringsize = <int>stringsize
        self._ptr = new C_Molecule(molstring,stringsize)
        
    def __init__(self, molstring,stringsize):
        pass
    
    def printmol(self):
        print("in molecule")
    
cdef class Substructure:
    cdef C_Substructure* _ptr
    
    def __cinit__(self, molstring, stringsize):
        if type(molstring) is unicode:
            molstring = (<unicode>molstring).encode('utf8')
        self._ptr = new C_Substructure(molstring,stringsize)
    
    def __init__(self, molstring, stringsize):
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

    def disconnect_bond(self, indices1, indices2):
        self._ptr[0].disconnect_bond(indices1, indices2)
        
    def connect_bond(self, indices1, indices2):
        self._ptr[0].connect_bond(indices1, indices2)
        
    def increaseBO(self, indices1, indices2, bondorder):
        self._ptr[0].increaseBO(indices1, indices2, bondorder)
        
    def decreaseBO(self, indices1, indices2, bondorder):
        self._ptr[0].decreaseBO(indices1, indices2, bondorder)

    def set_rule_name(self, rulename):
        if type(rulename) is unicode:
            rulename = (<unicode>rulename).encode('utf8')
        self._ptr[0].setRuleName(rulename)
         
cdef class RxnNetGen:
    cdef C_RxnNetGen* _ptr
    
    def __cinit__(self):
        self._ptr = new C_RxnNetGen()
    
    def add_initial_reactants(self, reactant_list):
        c_reactant_list = []
        for reactant in reactant_list:
            if type(reactant) is unicode:
                reactant = (<unicode>reactant).encode('utf8')
                c_reactant_list.append(reactant)
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
        
    def set_calc_thermo(self, boolvalue):
        self._ptr[0].SetCalcThermo(boolvalue)
        
    def add_lumping_strategy(self, lumpconstraint):
        self._ptr[0].AddLumpingStrategy(<C_LumpingStrategy?>lumpconstraint)
        
    def add_composite_sites(self, compositesites):
        self._ptr[0].AddCompositeSites(<vector[pair[string, c_SiteType] ]?>compositesites)
        
    def check_global_constraints(self, constraint):
        global_constraint = constraint
        
cdef public bool check_reactant_constraint0(C_Molecule mol, int reactiontype, int reactantindex):
#mol is the molecule to check for
#x is the reactiontype in list of reaction rules
#y is going to be 0    
    print("In constraint check!")
    print(reactiontype)
    #C_ReactionType = Rtlist[x]
    if reactiontype == 6:
        print("found rule 6")
        print(reactantindex)
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
    elif reactiontype == 7:
        print("found rule 7")
        print(reactantindex)
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

cdef public bool check_reactant_constraint1(C_Molecule mol, int reactiontype, int reactantindex):
    return True
    
cdef public bool check_combined_constraint(C_Molecule mol0, C_Molecule mol1, 
                                           int reactiontype):
    return True

cdef public bool check_combined_constraint_for_reactiontype(
        C_ReactionType reactiontype, C_Molecule mol0, C_Molecule mol1):
    return True

cdef public int check_combined_constraints_for_Rtlist(int reactiontype):
    return 1

cdef public bool check_product_constraint(C_Molecule mol, int reactiontype):
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


