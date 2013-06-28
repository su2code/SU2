########################################################################
# Pointwise - Proprietary software product of Pointwise, Inc.
#           Copyright (c) 1995-2009 Pointwise, Inc.
#           All rights reserved.
#
# Plugin SDK Makefile
########################################################################


# default to BUILD=Debug
ifndef BUILD
    BUILD := Debug
endif

# Common defines
ifeq ($(BUILD),Debug)
	DBG_SUFFIX := d
else
	DBG_SUFFIX :=
endif

SHLIB_SUFFIX := so
TOP_DIR := .
DIST_DIR := ./dist
DIST_BIN_DIR := $(DIST_DIR)/$(machine)/bin
DIST_LIB_DIR := $(DIST_DIR)/$(machine)/lib


########################################################################
# startup environment checking
#
KNOWN_BUILD_TYPES := Release Debug
KNOWN_SDK_MACHINES := linux linux_x86_64 macosx

ERROR_SHOWUSAGE = "FALSE"

ifeq ($(findstring $(machine),$(KNOWN_SDK_MACHINES)),)
    $(warning *** INVALID "machine=$(machine)")
    ERROR_SHOWUSAGE = "TRUE"
endif

ifeq ($(findstring $(BUILD),$(KNOWN_BUILD_TYPES)),)
    $(warning *** INVALID "BUILD=$(BUILD)")
    ERROR_SHOWUSAGE = "TRUE"
endif


ifeq ($(ERROR_SHOWUSAGE),"TRUE")
    $(warning )
    $(warning USAGE:)
    $(warning )
    $(warning gmake [machine=MID] [BUILD=BID] [target])
    $(warning )
    $(warning where,)
    $(warning MID is one of: $(KNOWN_SDK_MACHINES))
    $(warning BID is one of: $(KNOWN_BUILD_TYPES))
    $(warning )
    $(warning You can also set the "machine" and "BUILD" values in the shell)
    $(warning environment before running gmake.)
    $(warning )
    $(error *** MAKE TERMINATED ***)
endif



.PRECIOUS: %.cxx %.d %_moc.cxx %.moc

########################################################################
# Tool Macros
#
SHELL := /bin/sh
MKDIR := /bin/mkdir -p
RM := /bin/rm
RMF := $(RM) -f
RMR := $(RM) -rf
CP := /bin/cp -f
CPR := /bin/cp -fR


########################################################################
# Order is important.  Place in order of dependency.
ALL_MODULES := \
    plugins \
    $(NULL)


########################################################################
# Expand the names of the modules to their directories.
#
MODULES := $(patsubst plugins,./src/plugins, $(ALL_MODULES))


########################################################################
# Macros containing a list targets that only perform maintenance. No
# binaries are generated.
########################################################################

# Targets that remove the dependency and object files
CLEAN_TARGETS := $(patsubst %,%_clean,$(ALL_MODULES))

# Targets that remove the dependency files
CLEANDEP_TARGETS := $(patsubst %,%_cleandep,$(ALL_MODULES))

# Targets that remove the dependency, object, and library files
DISTCLEAN_TARGETS := $(patsubst %,%_distclean,$(ALL_MODULES))


.SUFFIXES: .i

########################################################################
# include platform specific build settings
include FLAGS.$(machine)

########################################################################
# The deafult goal for PluginSDK builds
#
.DEFAULT_GOAL := plugins

ifeq ($(MAKECMDGOALS),)
    CMDGOALS := $(.DEFAULT_GOAL)
else
    CMDGOALS := $(MAKECMDGOALS)
endif


########################################################################
# build the list of targets that do NOT match "print.*"
NON_PRINT_TARGETS := $(filter-out print.%,$(CMDGOALS))

# build the list of targets that DO match "print.*"
PRINT_TARGETS := $(filter-out $(NON_PRINT_TARGETS),$(CMDGOALS))


########################################################################
# Check on whether or not to include the dependencies.  If we're just doing
# one of the maintenance operations where we delete files, we don't want to
# include the dependencies which may trigger building them.
MAINT_TARGETS := \
    clean \
    cleandep \
    distclean \
    clobber \
    nuke \
    clobber-dr \
    nuke-dr \
    %-rd \
    %-dr \
    print.% \
    $(CLEAN_TARGETS) \
    $(CLEANDEP_TARGETS) \
    $(DISTCLEAN_TARGETS) \
    $(NULL)

.PHONY: \
    $(MAINT_TARGETS) \
    $(ALL_MODULES) \
    plugins \
    default \
    all \
    dependencies \
    install_validate \
    %_install \
    %_uninstall \
    $(NULL)


########################################################################
# Strip maintenance targets from the command-line goals.

NONCLEANGOALS := $(filter-out $(MAINT_TARGETS) $(PRINT_TARGETS),$(CMDGOALS))


########################################################################
# Include the description for each module

include $(patsubst %,%/module.mk,$(MODULES))


########################################################################
#                        GENERIC BUILD TARGETS
########################################################################

default: plugins

all: plugins

########################################################################
#                        WILDCARD BUILD TARGETS
########################################################################

#----------------------------------------------------------------------------
# The targets below are *not* BUILD dependent! Map them to the real target.
# If you add a new override target here, be sure to also add it to
# MAINT_TARGETS or .PHONY as appropriate. If you don't, gmake will match
# the %-dr target instead of your new override target.

clobber-dr: clobber
nuke-dr: nuke

# End of targets that are *not* BUILD dependent! See comments above...
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# The wildcard "xxxx-dr" target handler. Builds target xxxx with BUILD=Release
# and then Builds target xxxx with BUILD=Debug. For example;
#
#   gmake CaeTest_clean-dr CaeTest-dr
#
# is equivalent to:
#
#   gmake BUILD=Release CaeTest_clean
#   gmake BUILD=Debug CaeTest_clean
#   gmake BUILD=Release CaeTest
#   gmake BUILD=Debug CaeTest
# 
%-dr:
	@echo ""
	@echo "----- Release $* -----"
	-@$(MAKE) --no-print-directory BUILD=Release $*
	@echo ""
	@echo "----- Debug $* -----"
	-@$(MAKE) --no-print-directory BUILD=Debug $*

# typo handler redirects to %-dr target
%-rd:
	-@$(MAKE) --no-print-directory $*-dr


########################################################################
# Load (OPTIONAL) locally defined make targets. If a plugin developer wants
# to extend the global make scheme, they should create a local.mk file
# in the plugin SDK base folder. To provide for future SDK upgrades, the
# standard Makefile file should NOT be directly edited.
#
ifneq ($(wildcard local.mk),)
include local.mk
endif


########################################################################
#                           MAINT TARGETS
########################################################################

clean: $(CLEAN_TARGETS)

cleandep:  $(CLEANDEP_TARGETS)

distclean:  $(DISTCLEAN_TARGETS)
	-$(RM) -rf dist/$(machine)/bin

clobber:
	-$(RM) -rf {dist,src}/$(machine)

nuke: clobber
	-$(RM) -rf dist


########################################################################
#                           HELPER TARGETS
########################################################################

# This may be set in FLAGS.$(machine)
PLUGIN_INSTALL_RELPATH ?= ./plugins
PLUGIN_INSTALL_FULLPATH := $(PW_RUNTIME_ROOT)/$(PLUGIN_INSTALL_RELPATH)

install_validate:
ifeq ($(PW_RUNTIME_ROOT),)
	@echo "Cannot install: PW_RUNTIME_ROOT is not defined!"
else
 ifeq ($(wildcard $(PW_RUNTIME_ROOT)),)
	@echo "Cannot install: PW_RUNTIME_ROOT is invalid!"
 else
  ifeq ($(wildcard $(PLUGIN_INSTALL_FULLPATH)),)
	@echo "Cannot install: PLUGIN_INSTALL_RELPATH is invalid!"
  endif
 endif
endif

PLUGIN_INSTALL_TARGETS := $(patsubst %,%_install,$(ALL_PLUGIN_MODULES))
PLUGIN_UNINSTALL_TARGETS := \
    $(patsubst %,%_uninstall,$(ALL_PLUGIN_MODULES))

install: $(PLUGIN_INSTALL_TARGETS)

%_install:
	@echo "Ignoring invalid install target '$@'"

uninstall: $(PLUGIN_UNINSTALL_TARGETS)

%_uninstall:
	@echo "Ignoring invalid uninstall target '$@'"

$(DIST_BIN_DIR):
	$(MKDIR) $(DIST_BIN_DIR)

$(DIST_LIB_DIR):
	$(MKDIR) $(DIST_LIB_DIR)

dependencies:
	@echo "Dependencies up to date."

# Used to echo a macro value. usage: gmake print.DIST_BIN_DIR
print.%:
	@echo "$*=[$($*)]"
