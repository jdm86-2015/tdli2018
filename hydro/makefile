
# compilation flags and compiler
CXX = gfortran
INC = -I./include
WARNINGS = -Wall
OTHER = -DPPT
TEST = "test.f90"
CXXFLAGS = -O3 $(WARNINGS) $(OTHER) #-lhdf5 -lhdf5_fortran

# specify all files related to the distribution of the code.
AUXFILES := makefile
PROJDIRS := "./include" "./src"

# specify executable
ifdef TEST
EXE := test
EFILE := "main.f90"
SRCFILES := $(shell find $(PROJDIRS) -type f ! -name $(EFILE) -name "*.f90")
else
EXE := hydro
TEST = "test.f90"
SRCFILES := $(shell find $(PROJDIRS) -type f ! -name $(TEST) -name "*.f90")
endif

HDRFILES := $(shell find $(PROJDIRS) -type f -name "*.h")

# specify the objects
OBJFILES := $(patsubst %.f90,%.o,$(SRCFILES))

# specify the dependencies
DEPFILES := $(patsubst %.f90,%.d,$(SRCFILES))

# specify all included files for the distribution
ALLFILES := $(SRCFILES) $(HDRFILES) $(AUXFILES) "./src/main.f90" "./src/test.f90"

.PHONY : all clean

all : $(EXE)

$(EXE) :  $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean :
	$(RM) $(wildcard $(OBJFILES) $(DEPFILES) $(EXE) hydro.tar.gz)

dist :
	tar czvf hydro.tar.gz $(ALLFILES)

-include $(DEPFILES)

%.o: %.f90 makefile
	$(CXX) $(CXXFLAGS) -cpp -MMD -MP -c $< -o $@
