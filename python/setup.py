import os
from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize, build_ext

react_cpp_sources = list(map(lambda x: os.path.join('../cpp', x), [
           'common.cpp',
           'atom.cpp', 
           'atomcontainer.cpp', 
           'stringreg.cpp',
           'molecule.cpp',
           'substructure.cpp',
           'patternmatch.cpp',
           'reaction.cpp',
           'generated_rxn.cpp',
           'lumping.cpp',
           'automorphs.cpp',
           'groupadditivity.cpp',
           'rng.cpp',
           'additionalfunc.cpp',
           'element.cpp'
           ]))
react_srcs = react_cpp_sources + [
        os.path.join('pyring' , 'reactiontype.pyx' )
        ]


react_ext = Extension('pyring.react', 
                      react_srcs,
                      include_dirs = ['../cpp']
                     )

#react_ext = Extension('pyring.reaction', 
#                      sources = [
#                                 'alliedclassimp.cpp',
#                                 'atomicconst.cpp',
#                                 'registryimp.cpp',
#                                 'reactiontypeimp.cpp',
#                                 'pyreactiontype.pyx'
#                      ],
#                      include_dirs = 'src'
#)


setup(
      name = 'pyring',
      version = '0.1',
      author = 'Udit Gupta, Bharat Medasani',
      author_email = 'uditgupta0912@gmail.com',
      ext_modules = cythonize([react_ext], annotate=True)
      #cmd_cls = {'build_ext', built_ext}
)
