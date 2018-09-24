
IDA_LIB=idas

CPPFLAGS = -g -I$(IDA_LIB)/include

RING_SRC_FILES = additionalfunc.cpp atom.cpp atomcontainer.cpp automorphs.cpp common.cpp element.cpp groupadditivity.cpp kineticModel.cpp lumping.cpp Mechanisms.cpp molecule.cpp pathways.cpp patternmatch.cpp qsarimplementation.cpp reaction.cpp generated_rxn.cpp stringreg.cpp rng.cpp substructure.cpp ThirdPartyFiles.cpp KineticsInfo.cpp

IDA_SRC_FILES = idaa.c idaa_io.c idas_band.c idas_bbdpre.c idas.c idas_dense.c idas_direct.c idas_ic.c idas_io.c idas_spbcgs.c idas_spgmr.c idas_spils.c idas_sptfqmr.c 
SUNDIALS_SRC_FILES = sundials_band.c sundials_dense.c sundials_direct.c sundials_iterative.c sundials_math.c sundials_nvector.c sundials_spbcgs.c sundials_spgmr.c sundials_sptfqmr.c 
NVECTOR_SRC_FILES = nvector_serial.c fnvector_serial.c


RING_OBJ_FILES = $(RING_SRC_FILES:.cpp=.o)
IDA_OBJ_FILES = $(IDA_SRC_FILES:%.c=$(IDA_LIB)/src/idas/%.o)
SUNDIALS_OBJ_FILES = $(SUNDIALS_SRC_FILES:%.c=$(IDA_LIB)/src/sundials/%.o)
NVECTOR_OBJ_FILES = $(NVECTOR_SRC_FILES:%.c=$(IDA_LIB)/src/nvec_ser/%.o)

all: $(RING_OBJ_FILES) $(IDA_OBJ_FILES) $(SUNDIALS_OBJ_FILES) $(NVECTOR_OBJ_FILES)
	@cp $(IDA_LIB)/src/idas/*.o .
	@cp $(IDA_LIB)/src/sundials/*.o .
	@cp $(IDA_LIB)/src/nvec_ser/*.o .

clean:
	rm -f *.o $(IDA_LIB)/src/idas/*.o $(IDA_LIB)/src/sundials/*.o $(IDA_LIB)/src/nvec_ser/*.o *~

# I have no idea why it needs the extra \\ when we do \\\\%
WIN = $(RING_SRC_FILES:.cpp=.obj) $(IDA_SRC_FILES:%.c=$(IDA_LIB)\\src\\idas\\\\%.obj) $(SUNDIALS_SRC_FILES:%.c=$(IDA_LIB)\\src\\sundials\\\\%.obj) $(NVECTOR_SRC_FILES:%.c=$(IDA_LIB)\\src\\nvec_ser\\\\%.obj)

echowindows:
	@echo $(WIN)

