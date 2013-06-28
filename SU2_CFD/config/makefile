#
# Overall makefile for SU2_CFD
# Author: F. Palacios
# Version: 1.0
# 

MAIN=			SU2_CFD

include $(SU2_HOME)/$(MAIN)/config/makefile.in

INC_COMMON=		$(SU2_HOME)/Common/include
SRC_COMMON=		$(SU2_HOME)/Common/src
INCLUDES_COMMON=	$(INC_COMMON)/config_structure.hpp $(INC_COMMON)/dual_grid_structure.hpp \
			$(INC_COMMON)/geometry_structure.hpp $(INC_COMMON)/primal_grid_structure.hpp \
			$(INC_COMMON)/sparse_structure.hpp $(INC_COMMON)/linear_solvers_structure.hpp \
			$(INC_COMMON)/option_structure.hpp $(INC_COMMON)/grid_movement_structure.hpp
SOURCES_COMMON=		$(SRC_COMMON)/config_structure.cpp $(SRC_COMMON)/dual_grid_structure.cpp \
			$(SRC_COMMON)/geometry_structure.cpp $(SRC_COMMON)/primal_grid_structure.cpp \
			$(SRC_COMMON)/sparse_structure.cpp $(SRC_COMMON)/linear_solvers_structure.cpp \
      $(SRC_COMMON)/grid_movement_structure.cpp

SRC_CFD=		$(SU2_HOME)/$(MAIN)/src
INC_CFD=		$(SU2_HOME)/$(MAIN)/include
INCLUDES_CFD=		$(wildcard $(INC_CFD)/*.hpp)
SOURCES_CFD=		$(wildcard $(SRC_CFD)/*.cpp)

OBJ=			$(SU2_HOME)/$(MAIN)/obj
OBJECTS_COMMON=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_COMMON)))
OBJECTS_CFD=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_CFD)))

BIN=			  $(SU2_HOME)/$(MAIN)/bin
BINARIES=		$(BIN)/$(MAIN)

.SECONDARY:
.SUFFIXES:
.PHONY: default help clean all

default: help

help:
	@echo 'You must specify a target between: make all, make clean '

all: $(BINARIES)

$(OBJ)/%.o: $(SRC_COMMON)/%.cpp $(INCLUDES_COMMON) makefile
	@echo "Compiling \""$@"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_CFD) $(INC_CGNS) -o $@ -c $<

$(OBJ)/%.o: $(SRC_CFD)/%.cpp $(INCLUDES_CFD) makefile
	@echo "Compiling \""$@"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_CFD) $(INC_CGNS) -o $@ -c $<

$(BIN)/%: $(OBJECTS_COMMON) $(OBJECTS_CFD)
	@echo "Creating "$@
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_CFD) $(INC_CGNS) -o $@ $(OBJECTS_COMMON) $(OBJECTS_CFD) $(LIB_CGNS)
	@echo "Copying \""$@"\" to Python directory"
	@-cp $@ $(SU2_HOME)/SU2Py/$(MAIN)

clean:
	@echo "Deleting temporary files"
	@-rm -f core 
	@-rm -f `(find . | grep "~")`
	@echo "Deleting objects"
	@-rm -f *.o $(OBJECTS_CFD)
	@-rm -f *.o $(OBJECTS_COMMON)
	@echo "Deleting binaries"
	@-rm -f $(BINARIES)
