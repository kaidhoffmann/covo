# Compiler and flags
CXX      = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp
LDFLAGS  = -fopenmp

# Object files
OBJ = main.o parameters.o directory.o healpix.o catalogue.o correlation.o

#------------------------------------------------------
# Build targets
#------------------------------------------------------

# Main executable
covo: $(OBJ)
	$(CXX) -o ../covo $(OBJ) $(LDFLAGS)

# Compile .cpp files to .o files
%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)

#------------------------------------------------------
# Utility targets
#------------------------------------------------------

.PHONY: clean test test-baseline

# Remove all build artifacts
clean:
	rm -f $(OBJ) ../covo

# Run regression tests (builds first if needed)
test: covo
	@cd .. && bash tests/run_regression_tests.sh

# Generate baseline for regression tests
test-baseline:
	@cd .. && bash tests/run_regression_tests.sh || true
	@echo "Baseline generation complete. Review tests/baseline/ before committing."
