/* Generated by Cython 0.28.5 */

#ifndef __PYX_HAVE__pyring__react
#define __PYX_HAVE__pyring__react

struct Bunny;

/* "pyring/reactiontype.pyx":233
 *     return True
 * 
 * cdef public struct Bunny: # public type declaration             # <<<<<<<<<<<<<<
 *     int vorpalness
 * 
 */
struct Bunny {
  int vorpalness;
};

#ifndef __PYX_HAVE_API__pyring__react

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

__PYX_EXTERN_C bool check_reactant_constraint0(Molecule, int, int);
__PYX_EXTERN_C bool check_reactant_constraint1(Molecule, int, int);
__PYX_EXTERN_C bool check_combined_constraint(Molecule, Molecule, int);
__PYX_EXTERN_C bool check_combined_constraint_for_reactiontype(Reactiontype, Molecule, Molecule);
__PYX_EXTERN_C int check_combined_constraints_for_Rtlist(int);
__PYX_EXTERN_C bool check_product_constraint(Molecule, int);
__PYX_EXTERN_C void grail(struct Bunny *);
__PYX_EXTERN_C bool check_global_constraints(Molecule);

__PYX_EXTERN_C int spam;

#endif /* !__PYX_HAVE_API__pyring__react */

/* WARNING: the interface of the module init function changed in CPython 3.5. */
/* It now returns a PyModuleDef instance instead of a PyModule instance. */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initreact(void);
#else
PyMODINIT_FUNC PyInit_react(void);
#endif

#endif /* !__PYX_HAVE__pyring__react */