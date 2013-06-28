#
# Overall makefile for SU2_MAC
# Author: Francisco Palacios
# Version: 1.0
#

MAIN=			SU2_MAC

include $(SU2_HOME)/$(MAIN)/config/makefile.in

INC_COMMON=		$(SU2_HOME)/Common/include
SRC_COMMON=		$(SU2_HOME)/Common/src
INCLUDES_COMMON=	$(INC_COMMON)/config_structure.hpp $(INC_COMMON)/dual_grid_structure.hpp \
			$(INC_COMMON)/geometry_structure.hpp $(INC_COMMON)/grid_adaptation_structure.hpp \
			$(INC_COMMON)/primal_grid_structure.hpp $(INC_COMMON)/option_structure.hpp
SOURCES_COMMON=		$(SRC_COMMON)/config_structure.cpp $(SRC_COMMON)/dual_grid_structure.cpp \
			$(SRC_COMMON)/geometry_structure.cpp $(SRC_COMMON)/grid_adaptation_structure.cpp \
			$(SRC_COMMON)/primal_grid_structure.cpp

SRC_MAC=		$(SU2_HOME)/$(MAIN)/src
INC_MAC=		$(SU2_HOME)/$(MAIN)/include
INCLUDES_MAC=		$(wildcard $(INC_MAC)/*.hpp)
SOURCES_MAC=		$(wildcard $(SRC_MAC)/*.cpp)

OBJ=			$(SU2_HOME)/$(MAIN)/obj
OBJECTS_COMMON=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_COMMON)))
OBJECTS_MAC=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_MAC)))

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
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_MAC) -o $@ -c $<

$(OBJ)/%.o: $(SRC_MAC)/%.cpp $(INCLUDES_MAC) makefile
	@echo "Compiling \""$@"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_MAC) -o $@ -c $<

$(BIN)/%: $(OBJECTS_COMMON) $(OBJECTS_MAC)
	@echo "Creating "$@
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_MAC) -o $@ $(OBJECTS_COMMON) $(OBJECTS_MAC) $(LDFLAGS)
	@echo "Copying \""$@"\" to Python directory"
	@-cp $@ $(SU2_HOME)/SU2Py/$(MAIN)

clean:
	@echo "Deleting temporary files"
	@-rm -f core 
	@-rm -f `(find . | grep "~")`
	@echo "Deleting objects"
	@-rm -f *.o $(OBJECTS_MAC)
	@-rm -f *.o $(OBJECTS_COMMON)
	@echo "Deleting binaries"
	@-rm -f $(BINARIES)
