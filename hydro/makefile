
# compilation flags and compiler
CXX = gfortran
INC = -I./include
WARNINGS = -Wall -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -fbacktrace
OTHER = -DPPT
CXXFLAGS = -O3 $(WARNINGS) $(OTHER) #-MMD -MP #-lhdf5 -lhdf5_fortran

SRCFILES := ./src/boundary.f90 \
			./src/conservedVars.f90 \
			./src/eos.f90 \
			./src/primatives.f90 \
			./src/reconstructor.f90 \
			./src/riemann.f90 \
			./src/sweepFunc.f90 \
			./src/timeStepCalc.f90 \
			./src/rkStepper.f90 \
			./src/setupGrid.f90 \
			./src/setSod.f90 \
			./src/output.f90 \
			./src/hydroDriver.f90 \

# specify the objects
OBJFILES := $(patsubst %.f90,%.o,$(SRCFILES))

# specify the dependencies
#DEPFILES := $(patsubst %.f90,%.d,$(SRCFILES))

# specify all included files for the distribution
ALLFILES := $(SRCFILES) "./src/main.f90" "./src/test.f90" makefile
 
EXE := hydro
TEST := testScript

.PHONY : all clean

all : $(EXE)

test : $(TEST)

$(EXE) :  $(OBJFILES) ./src/main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TEST) : $(OBJFILES) ./src/test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean :
	$(RM) $(wildcard ./src/*.o *.mod ./src/*.d hydro.tar.gz)

dist :
	tar czvf hydro.tar.gz $(ALLFILES)

-include $(DEPFILES)

%.o: %.f90 makefile
	$(CXX) $(CXXFLAGS) -cpp -c $< -o $@

./src/main.o: ./src/main.f90 makefile
	$(CXX) $(CXXFLAGS) -cpp -c $< -o $@

./src/test.o: ./src/test.f90 makefile
	$(CXX) $(CXXFLAGS) -cpp -c $< -o $@

