#
# Overall makefile for SU2_GPC
# Author: F. Palacios
# Version: 1.0
# 

MAIN=			SU2_GPC

include $(SU2_HOME)/$(MAIN)/config/makefile.in

INC_COMMON=		$(SU2_HOME)/Common/include
SRC_COMMON=		$(SU2_HOME)/Common/src
INCLUDES_COMMON=	$(INC_COMMON)/config_structure.hpp $(INC_COMMON)/dual_grid_structure.hpp \
			$(INC_COMMON)/geometry_structure.hpp  $(INC_COMMON)/primal_grid_structure.hpp  \
			$(INC_COMMON)/grid_movement_structure.hpp $(INC_COMMON)/sparse_structure.hpp \
			$(INC_COMMON)/linear_solvers_structure.hpp $(INC_COMMON)/option_structure.hpp
SOURCES_COMMON=		$(SRC_COMMON)/config_structure.cpp $(SRC_COMMON)/dual_grid_structure.cpp \
			$(SRC_COMMON)/geometry_structure.cpp $(SRC_COMMON)/primal_grid_structure.cpp \
			$(SRC_COMMON)/grid_movement_structure.cpp $(SRC_COMMON)/sparse_structure.cpp \
			$(SRC_COMMON)/linear_solvers_structure.cpp

SRC_GPC=		$(SU2_HOME)/$(MAIN)/src
INC_GPC=		$(SU2_HOME)/$(MAIN)/include
INCLUDES_GPC=		$(wildcard $(INC_GPC)/*.hpp)
SOURCES_GPC=		$(wildcard $(SRC_GPC)/*.cpp)

OBJ=			$(SU2_HOME)/$(MAIN)/obj
OBJECTS_COMMON=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_COMMON)))
OBJECTS_GPC=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_GPC)))

BIN=			$(SU2_HOME)/$(MAIN)/bin
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
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_GPC) $(INC_CGNS) -o $@ -c $<

$(OBJ)/%.o: $(SRC_GPC)/%.cpp $(INCLUDES_GPC) makefile
	@echo "Compiling \""$@"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_GPC) $(INC_CGNS) -o $@ -c $<

$(BIN)/%: $(OBJECTS_COMMON) $(OBJECTS_GPC)
	@echo "Creating "$@
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_GPC) $(INC_CGNS) $(LIB_CGNS) -o $@ $(OBJECTS_COMMON) $(OBJECTS_GPC) $(LDFLAGS)
	@echo "Copying \""$@"\" to Python directory"
	@-cp $@ $(SU2_HOME)/SU2Py/$(MAIN)

clean:
	@echo "Deleting temporary files"
	@-rm -f core 
	@-rm -f `(find . | grep "~")`
	@echo "Deleting objects"
	@-rm -f *.o $(OBJECTS_GPC)
	@-rm -f *.o $(OBJECTS_COMMON)
	@echo "Deleting binaries"
	@-rm -f $(BINARIES)
