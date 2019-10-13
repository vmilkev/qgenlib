.RECIPEPREFIX = >
# define the C compiler to use
CXX = icc

MATRLIB = ../../mathlibs/lib/matrix.a

# define any compile-time flags
CXXFLAGS = -g -O0 -qopenmp -traceback -check-pointers=rw -debug all -debug parallel -mcmodel=large -I${MATRLIB} -I$(CATCH_SINGLE_INCLUDE) -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits
# CXXFLAGS = -O3 -qopenmp -mcmodel=large -I${MATRLIB} -I$(CATCH_SINGLE_INCLUDE) -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits

# define any directories containing header files other than /usr/include
#CPPFLAGS = -I./include
CPPFLAGS = -cxxlib -static-libgcc -static-libstdc++ -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl

#define the directive for object files
OBJDIR = ./lib
SRCDIR = ./src
BINDIR = ./bin

# define the C source files
SRCS = $(SRCDIR)/catchmain.cpp $(SRCDIR)/test_matrix_float.cpp $(SRCDIR)/test_matrix_double.cpp

# define the C object files 
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# define the executable file 
MAIN = $(BINDIR)/test.exe

.PHONY: depend 

all: $(MAIN)
>@echo Program compiled

$(MAIN): $(OBJS)
>$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(MAIN) $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
>@echo "Compiling: " $@
>$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
>$(RM) $(OBJDIR)/*.o *~ $(MAIN)

depend: $(SRCS)
>makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it