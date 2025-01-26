# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
I := $(patsubst %/,%,$(dir $(_mkfile_path)))

ifneq ($(words $(MAKECMDGOALS)),1)
.DEFAULT_GOAL = all
%:
	@$(MAKE) $@ --no-print-directory -rRf $(firstword $(MAKEFILE_LIST))
else
ifndef ECHO
T := $(shell $(MAKE) $(MAKECMDGOALS) --no-print-directory \
     -nrRf $(firstword $(MAKEFILE_LIST)) \
     ECHO="COUNTTHIS" | grep -c "COUNTTHIS")
N := x
C = $(words $N)$(eval N := x $N)
ECHO = python3 $(I)/makefile_progress.py --stepno=$C --nsteps=$T
endif

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all all_libs exe_targets clean

.SILENT:

# Targets are grouped here

all: all_libs exe_targets
	@echo "All done"

exe_targets: SigmalizedResiduals 
all_libs: 	 CppToolsLib ROOTToolsLib PBarLib CalPhenixLib

# CppTools target groups

CppToolsLib: 	 	 ErrorHandler StrTools IOTools Box FitTools
ErrorHandler: 	 	 $(CPP_TOOLS_LIB_DIR)/libErrorHandler.so
StrTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/libStrTools.so
IOTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/libIOTools.so
Box: 			  	 	 $(CPP_TOOLS_LIB_DIR)/libBox.so

# ROOTTools target groups

ROOTToolsLib: 	 	 TCanvasPrinter
TCanvasPrinter: 	 $(ROOT_TOOLS_LIB_DIR)/libTCanvasPrinter.so
FitTools: 	       $(ROOT_TOOLS_LIB_DIR)/libFitTools.so

# ProgressBar target groups

PBarLib: 		 	 PBar
PBar: 			 	 $(PBAR_LIB_DIR)/libPBar.so

# current repository target groups (Analysis)

CalPhenixLib: 		 InputReader
InputReader:	 	 lib/libInputReader.so

# CppTools sumbodule targets

$(CPP_TOOLS_LIB_DIR): 
	mkdir -p $@

$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o: $(CPP_TOOLS_SRC_DIR)/ErrorHandler.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/StrTools.o: $(CPP_TOOLS_SRC_DIR)/StrTools.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/IOTools.o: $(CPP_TOOLS_SRC_DIR)/IOTools.cpp | ErrorHandler StrTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools

$(CPP_TOOLS_LIB_DIR)/Box.o: $(CPP_TOOLS_SRC_DIR)/Box.cpp | IOTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/Table.o: $(CPP_TOOLS_SRC_DIR)/Table.cpp | IOTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/lib%.so: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# ROOTTools sumbodule targets

$(ROOT_TOOLS_LIB_DIR):
	mkdir -p $@

$(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.o: $(ROOT_TOOLS_SRC_DIR)/TCanvasPrinter.cpp | \
												    $(ROOT_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

$(ROOT_TOOLS_LIB_DIR)/FitTools.o: $(ROOT_TOOLS_SRC_DIR)/FitTools.cpp | \
											 $(ROOT_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

$(ROOT_TOOLS_LIB_DIR)/lib%.so: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# ProgressBar sumbodule targets

$(PBAR_LIB_DIR):
	mkdir -p $@

$(PBAR_LIB_DIR)/PBar.o: $(PBAR_SRC_DIR)/PBar.cpp | $(PBAR_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(PBAR_LIB_DIR)/lib%.so: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# Current repository targets (Analysis)

lib:
	mkdir -p $@

lib/InputReader.o: src/InputReader.cpp | lib ErrorHandler
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(JSON_INCLUDE) $(JSON_LIB) \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB)

lib/lib%.so: lib/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

SigmalizedResiduals: src/SigmalizedResiduals.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ALL_INCLUDE) $(ALL_LIB)

# other

bin:
	mkdir $@

clean: 
	@echo Cleaning
	rm -rf bin ; \
	rm -rf CppTools/lib ; \
	rm -rf ROOTTools/lib ; \
	rm -rf ProgressBar/lib ; \
	rm -rf lib

endif
