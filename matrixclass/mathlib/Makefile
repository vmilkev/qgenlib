.RECIPEPREFIX = >
# define the C compiler to use
CXX = icc

# define any compile-time flags
CXXFLAGS = -g -O0 -qopenmp -traceback -check-pointers=rw -debug all -debug parallel -mcmodel=large -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits
# CXXFLAGS = -O3 -qopenmp -mcmodel=large -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits

# define any directories containing header files other than /usr/include
#CPPFLAGS = -I./include

#define the directive for object files
OBJDIR = ./lib
SRCDIR = ./src

# define the C source files
SRCS = $(SRCDIR)/cs_matrix.cpp

# define the C object files 
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# define the executable file 
MAIN = $(OBJDIR)/matrix.a

.PHONY: depend 

all: $(MAIN)
>@echo Program compiled

$(MAIN): $(OBJS)
>ar rcs $(MAIN) $(OBJS)
#    $(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(MAIN) $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
>@echo "Compiling: " $@
>$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
>$(RM) $(OBJDIR)/*.o *~ $(MAIN)

depend: $(SRCS)
>makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it