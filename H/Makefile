
include make.inc

ALL=	hmat

all: $(ALL)

#================================================================

#OBJhmat = structures.o cs_matrix.o cs_gmat.o cs_matrixmath.o cs_housekeeper.o cs_apy.o cs_interface.o main_hmat.o
OBJhmat = main_hmat.o cs_interface.o cs_apy.o cs_housekeeper.o cs_matrixmath.o cs_gmat.o structures.o

#================================================================

hmat: $(OBJhmat)
	$(CXX) $(CppLINK) -Xlinker --start-group -Xlinker $(addprefix $(MYLIB),$(OBJhmat)) -Xlinker --end-group -o $(BIN)/hmat.exe

-include $(OBJhmat:%.o=%.d)

#================================================================

.SUFFIXES:	.cpp
.cpp.o :
	${CXX} ${CppFLAGS} -c $^ -o $(addprefix $(MYLIB),$@)

#================================================================

CLEANOBJ = src/*.o src/*.lst lib/*.o lib/*.d bin/*.exe
clean:
	rm -f *.o core *.lst *.mod
	rm -f $(CLEANOBJ)
