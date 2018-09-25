
# Using a .pxd file gives us a separate namespace for
# the C++ declarations. Using a .pxd file also allows
# us to reuse the declaration in multiple .pyx modules.
from libcpp.list cimport list
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp cimport bool
from lumping cimport C_LumpingStrategy

cdef extern from "additionalfunc.h":
    cdef int c_patternsize "patternsize"(string)

cdef extern from "common.h":
    ctypedef enum c_SiteType "SiteType":
        Homogeneous
        Heterogeneous
    
cdef extern from "molecule.h":
    cdef cppclass C_Molecule "Molecule":
        C_Molecule () except +
        C_Molecule (string, int) except +
    
    ctypedef bool (*ConstrPtr)(const C_Molecule &)
    ctypedef bool (*CombinedConstrPtr)(const C_Molecule &, const C_Molecule &)
    
cdef extern from "substructure.h":
    cdef cppclass C_Substructure "Substructure":
        C_Substructure () except +
        C_Substructure (string, int) except +
        
cdef extern from "patternmatch.h":
    cdef cppclass C_Patternmatch "Patternmatch":
        C_Patternmatch() except +
        C_Patternmatch (const C_Molecule&, const C_Substructure&, int) except +
        int GetDistinctMatches()

cdef extern from "reaction.h":
    cdef cppclass C_ReactionType "Reactiontype":
        C_ReactionType() except +
        vector [C_Substructure] reactant_pattern
        void add_reactant_pattern(C_Substructure)
        void AddReactantCopy(int, map[int,int])
        void add_reactantconstraint(ConstrPtr)
        void add_productconstraint(ConstrPtr)
        void add_combined_constraint(CombinedConstrPtr)
        void disconnect_bond(int,int)
        void connect_bond(int,int)
        void increaseBO(int,int, int)
        void decreaseBO(int,int, int)   
        void setRuleName(string)
        
cdef extern from "rng.h":
    cdef cppclass C_RxnNetGen "rxn_net_gen":
        C_RxnNetGen() except +
        C_RxnNetGen(list[string]&, vector[C_ReactionType]&, ConstrPtr, 
                    C_LumpingStrategy&, vector[string]&, vector[string]&, 
                    bool) except +
        void GenerateNetwork()
        void AddInitialReactants(list[string]&)
        void AddReactionRules(vector[C_ReactionType]&)
        void AddGlobalConstraints(ConstrPtr)
        void AddLumpingStrategy(C_LumpingStrategy&)
        void AddCompositeAtoms(vector[string]&)
        void AddCompositeSites(vector[pair[string, c_SiteType] ]&)
        void SetAllLumpedRxns()
        void print_rxnlist()
        void SetCalcThermo(bool)
 