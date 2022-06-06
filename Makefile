#============================================================
#
# Makefile to compile this code
#
#============================================================

# Use bash since sh on datastar does not recognize ">&" used in dep: target
SHELL       = /bin/bash 

# Now let us add all code files
CPARAMS  = define.h
CCLASSES = Classes/cellclass.h Classes/fluxclass.h Classes/internomials.h Classes/Splineclass.h
CPHYSICS = Physics/gasPhysics.cpp Physics/thermodynamics.cpp Physics/hydrodynamics.cpp
CSCHEMES = Schemes/timeEvolution.cpp Schemes/boundaries.cpp
CAMR	 = Schemes/AMR_aux.cpp Schemes/AMR_main.cpp
CIO      = IO/read_input.cpp IO/set_ambient.cpp IO/writeFile.cpp
COTHER   = IO/ConversorString.cpp Schemes/interpolations.cpp Schemes/refinement_criteria.cpp

CALL = $(CPARAMS) $(CCLASSES) $(CPHYSICS) $(CSCHEMES) $(CAMR) $(CIO) $(COTHER)

#-----------------------------------------------------------------------
# add C++ standard used
#-----------------------------------------------------------------------

LDFLAGS  += -std=c++11
CXXFLAGS += -std=c++11

#=======================================================================
# Now the compilation

main: $(CALL) main.o
	
	@echo "Compiling..."
	-@$(CXX) -std=c++11 -fPIC $(LDFLAGS)  -o main.exe main.o $(LIBS) $(GRACKLE_LIB)
	@echo "Done. Creating Outputs folder..."
	-@mkdir Outputs
	@echo "Finished!"

remove:

	@echo "Removing text files..."
	-@rm *.txt
	@echo "Done!"

clean:

	@echo "Removing all .o and .exe"
	-@rm *.o *.exe *.txt
	@echo "Done!"

# That's all
#=======================================================================
