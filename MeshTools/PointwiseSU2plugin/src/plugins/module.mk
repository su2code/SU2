########################################################################
# Pointwise - Proprietary software product of Pointwise, Inc.
#             Copyright (c) 1995-2012 Pointwise, Inc.
#             All rights reserved.
#
# the master module.mk for all src\plugins\ targets
########################################################################

########################################################################
#
# THIS FILE IS USED FOR *NIX BUILDS ONLY!
#
# You may need equivalent changes in .../Externals/*.mk for winOS builds!
#
########################################################################

########################################################################
# Module defines
#
PLUGINS_LOC := src/plugins
PLUGINS_OBJ_LOC := src/$(machine)/$(BUILD)/plugins

# KNOWN_SDK_MACHINES is only defined when building from SDK distro.
# This will cause the ifeq to fail. When building from SDK, we do NOT
# want to copy the plugins to the macos bundle. As an enhancement, we
# can support a way to specify where the Pointwise bundle is located so
# the plugin can be automatically installed after building.
#
ifeq ("$(machine)$(KNOWN_SDK_MACHINES)","macosx")
    PLUGINS_DIST_DIR := $(DIST_DIR)/$(machine)/$(APP_NAME)/Contents/Plugins
else
    ifeq ($(BUILD),Debug)
	PLUGINS_DIST_DIR := $(DIST_DIR)/$(machine)/plugins/debug
    else
	PLUGINS_DIST_DIR := $(DIST_DIR)/$(machine)/plugins
    endif
endif


########################################################################
# The list of plugin targets. Add targets to the modulelocal.mk file.
#
#          See PluginSDK/src/plugins/README.txt
#
ALL_PLUGIN_MODULES := \
	$(NULL)

# Get (required) locally defined make targets.
ifneq ($(wildcard $(PLUGINS_LOC)/modulelocal.mk),)
include $(PLUGINS_LOC)/modulelocal.mk
endif


########################################################################
# Runtime helper module defines
#
# list of all shared runtime API helper modules
PLUGINS_RT_MODULES := \
	CAEP \
	PWGM \
	PWP \
	$(NULL)

# root-level shared template folder
PLUGINS_RT_LOC := $(PLUGINS_LOC)/shared

PLUGINS_RT_INCL := \
	-I$(PLUGINS_LOC) \
	$(patsubst %,-I$(PLUGINS_RT_LOC)/%,$(PLUGINS_RT_MODULES)) \
	$(NULL)

PLUGINS_STDDEFS = \
	-DBUILD_PWPLUGIN_DYNLIB

########################################################################
# Each API's shared src files. Do NOT put the subfolder prefix on
# the base file names. The subfolder should be added in the corresponding
# XX_RT_XXXSRC macro expansion. These file lists are used by the each
# plugins module.mk file to properly build their required .o files.
#
PLUGINS_RT_PWPLOC := $(PLUGINS_RT_LOC)/PWP
PLUGINS_RT_PWPFILES := \
	apiPWP.c \
	apiPWPUtils.c \
	pwpPlatform.c \
	$(NULL)

PLUGINS_RT_PWPSRC := \
	$(patsubst %,$(PLUGINS_RT_PWPLOC)/%,$(PLUGINS_RT_PWPFILES))

PLUGINS_RT_PWGMLOC := $(PLUGINS_RT_LOC)/PWGM
PLUGINS_RT_PWGMFILES := \
	apiGridModel.c \
	$(NULL)

PLUGINS_RT_PWGMSRC := \
	$(patsubst %,$(PLUGINS_RT_PWGMLOC)/%,$(PLUGINS_RT_PWGMFILES))

PLUGINS_RT_CAEPLOC := $(PLUGINS_RT_LOC)/CAEP
PLUGINS_RT_CAEPFILES := \
	apiCAEP.c \
	apiCAEPUtils.c \
	$(NULL)

PLUGINS_RT_CAEPSRC := \
	$(patsubst %,$(PLUGINS_RT_CAEPLOC)/%,$(PLUGINS_RT_CAEPFILES))


########################################################################
# Expand the names of the modules to their directories.
#
PLUGIN_MODULE_PATHS := $(patsubst %,$(PLUGINS_LOC)/%,$(ALL_PLUGIN_MODULES))


########################################################################
# The master target
#
plugins: $(ALL_PLUGIN_MODULES)

# Remove the dependency and object files
PLUGIN_CLEAN_TARGETS := $(patsubst %,%_clean,$(ALL_PLUGIN_MODULES))

# Remove the dependency files
PLUGIN_CLEANDEP_TARGETS := $(patsubst %,%_cleandep,$(ALL_PLUGIN_MODULES))

# Remove the dependency, object, and library files
PLUGIN_DISTCLEAN_TARGETS := $(patsubst %,%_distclean,$(ALL_PLUGIN_MODULES))

plugins_cleandep: $(PLUGIN_CLEANDEP_TARGETS)

plugins_clean: $(PLUGIN_CLEAN_TARGETS)

plugins_distclean: $(PLUGIN_DISTCLEAN_TARGETS)

########################################################################
# The SDK targets
#
SDK_SRCDIR := $(PLUGINS_LOC)/PluginSDK
SDK_DISTDIR := $(DISTSDK_DIR)/PluginSDK
SDK_REL_DISTDIR := $(SDK_DISTDIR)/src/plugins

SDK_REL_SRCDIRS := \
	$(PLUGINS_LOC)/shared/CAEP \
	$(PLUGINS_LOC)/shared/PWGM \
	$(PLUGINS_LOC)/shared/PWP \
	$(PLUGINS_LOC)/templates/CAEP \
	$(PLUGINS_LOC)/templates/PWP \
	$(PLUGINS_LOC)/CaeUnsXML \
	$(PLUGINS_LOC)/CaeStrXML
SDK_REL_DISTDIRS := $(SDK_REL_SRCDIRS:$(PLUGINS_LOC)/%=$(SDK_REL_DISTDIR)/%)

SDK_VS_SRCDIR := $(PLUGINS_LOC)/../Pointwise
SDK_VS_DISTDIR := $(SDK_DISTDIR)/src/Pointwise

SDK_DISTDIRS := $(SDK_VS_DISTDIR) $(SDK_REL_DISTDIRS)

$(DISTSDK_DIR):
	@echo "***"
	@echo "*** Create Plugin SDK Distribution"
	@echo "***"
	-@$(MKDIR) $@

$(SDK_DISTDIRS):
	-@$(MKDIR) $@

# copy PluginSDK files
SDK_SRC_FILES := \
	$(filter-out $(SDK_SRCDIR)/src,$(wildcard $(SDK_SRCDIR)/*)) \
	$(wildcard $(SDK_SRCDIR)/src/plugins/*)
SDK_DIST_FILES := $(SDK_SRC_FILES:$(SDK_SRCDIR)/%=$(SDK_DISTDIR)/%)

$(SDK_DIST_FILES): $(SDK_DISTDIR)/%: $(SDK_SRCDIR)/%
	$(CP) $< $@
	@chmod ug+w $@

# copy vsprops files
SDK_VS_SRC_FILES := $(wildcard $(SDK_VS_SRCDIR)/Win??Target.vsprops)
SDK_VS_DIST_FILES := $(SDK_VS_SRC_FILES:$(SDK_VS_SRCDIR)/%=$(SDK_VS_DISTDIR)/%)

$(SDK_VS_DIST_FILES): $(SDK_VS_DISTDIR)/%: $(SDK_VS_SRCDIR)/%
	$(CP) $< $@
	@chmod ug+w $@

# copy SDK source files
SDK_REL_SRC_FILES := $(PLUGINS_LOC)/module.mk \
	$(PLUGINS_LOC)/README.txt \
	$(PLUGINS_LOC)/structured.vsprops \
	$(PLUGINS_LOC)/unstructured.vsprops \
	$(foreach dir,$(SDK_REL_SRCDIRS),$(wildcard $(dir)/*))
SDK_REL_DIST_FILES := $(SDK_REL_SRC_FILES:$(PLUGINS_LOC)/%=$(SDK_REL_DISTDIR)/%)

$(SDK_REL_DIST_FILES): $(SDK_REL_DISTDIR)/%: $(PLUGINS_LOC)/%
	$(CP) $< $@
	@chmod ug+w $@

# SDK target
PluginSDK: \
	$(DISTSDK_DIR) \
	$(SDK_DISTDIRS) \
	$(SDK_DIST_FILES) \
	$(SDK_REL_DIST_FILES) \
	$(SDK_VS_DIST_FILES)

# SDK forced copy
PluginSDK_force:
	@echo ""
	@echo "***"
	@echo "*** Constructing $@ Distribution in $(DISTSDK_DIR)/"
	@echo "***"
	$(RMR) $(DISTSDK_DIR)/PluginSDK
	$(MKDIR) $(DISTSDK_DIR)
	$(CPR) $(PLUGINS_LOC)/PluginSDK $(DISTSDK_DIR)/PluginSDK
	$(MKDIR) $(DISTSDK_DIR)/PluginSDK/src/Pointwise
	$(CP) $(PLUGINS_LOC)/../Pointwise/Win??Target.vsprops $(DISTSDK_DIR)/PluginSDK/src/Pointwise/.
	$(CP) $(PLUGINS_LOC)/module.mk $(DISTSDK_DIR)/PluginSDK/src/plugins/.
	$(CP) $(PLUGINS_LOC)/README.txt $(DISTSDK_DIR)/PluginSDK/src/plugins/.
	$(CP) $(PLUGINS_LOC)/structured.vsprops $(DISTSDK_DIR)/PluginSDK/src/plugins/.
	$(CP) $(PLUGINS_LOC)/unstructured.vsprops $(DISTSDK_DIR)/PluginSDK/src/plugins/.
	$(CPR) $(PLUGINS_LOC)/shared $(DISTSDK_DIR)/PluginSDK/src/plugins/shared
	$(CPR) $(PLUGINS_LOC)/templates $(DISTSDK_DIR)/PluginSDK/src/plugins/templates
	$(CPR) $(PLUGINS_LOC)/CaeUnsXML $(DISTSDK_DIR)/PluginSDK/src/plugins/CaeUnsXML
	$(CPR) $(PLUGINS_LOC)/CaeStrXML $(DISTSDK_DIR)/PluginSDK/src/plugins/CaeStrXML
	find $(DISTSDK_DIR)/. -type f -print | xargs chmod 0664
	find $(DISTSDK_DIR)/. -type d -print | xargs chmod 0774
	chmod ug+x $(DISTSDK_DIR)/PluginSDK/mkplugin
	chmod ug+x $(DISTSDK_DIR)/PluginSDK/depend.sh
	@echo "***"
	@echo "*** END $@ Distribution"
	@echo "***"

PluginSDK_distclean:
	$(RMR) $(DISTSDK_DIR)/PluginSDK

########################################################################
# Required standard targets implemented in a given plugin's module.mk file
#
CLEANDEP_TARGETS += $(PLUGIN_CLEANDEP_TARGETS)
CLEAN_TARGETS += $(PLUGIN_CLEAN_TARGETS)
DISTCLEAN_TARGETS += $(PLUGIN_DISTCLEAN_TARGETS)

########################################################################
# Full list of standard plugin maintenance targets filtered out of
# PLUGIN_NONCLEANGOALS below. PLUGIN_MAINT_TARGETS can/may be
# modified in a given plugin's module.mk file included below.
#
PLUGIN_MAINT_TARGETS = \
	PluginSDK \
	PluginSDK_distclean \
	$(PLUGIN_CLEANDEP_TARGETS) \
	$(PLUGIN_CLEAN_TARGETS) \
	$(PLUGIN_DISTCLEAN_TARGETS) \
	$(NULL)

.PHONY: \
    default \
	$(PLUGIN_MAINT_TARGETS) \
	$(NULL)


########################################################################
# *all* targets that will trigger a plugin build target. These targets require
# the loading of dependencies.
#
ALL_PLUGIN_TARGETS := plugins $(ALL_PLUGIN_MODULES) $(PLUGIN_MAINT_TARGETS)


# Strip plugin goals from CMDGOALS leaving only non-plugin related goals.
# We consider 'default all' to be a plugin goals here so that the dependencies
# will be loaded.
#
NON_PLUGIN_CMDGOALS := $(filter-out default all $(ALL_PLUGIN_TARGETS), $(CMDGOALS))


# strip NON_PLUGIN_CMDGOALS from CMDGOALS leaving only *plugin* related goals
#
PLUGIN_CMDGOALS := $(filter-out $(NON_PLUGIN_CMDGOALS), $(CMDGOALS))


########################################################################
# Each plugin module adds to this for calculating the dependencies. Each object
# file will result in a corresponding dependency file.
#
PLUGIN_OBJ :=


########################################################################
# Include the module.mk for each plugin target
include $(patsubst %,%/module.mk,$(PLUGIN_MODULE_PATHS))


# Strip out plugin targets from top-level NONCLEANGOALS. They do NOT ever need
# Pointwise external dependencies. plugins/module.mk handles the loading of
# plugin dependencies.
#
NONCLEANGOALS := $(filter-out $(ALL_PLUGIN_TARGETS),$(NONCLEANGOALS))

# these targets require the plugin dependencies
PLUGIN_NONCLEANGOALS := $(filter-out $(PLUGIN_MAINT_TARGETS), $(PLUGIN_CMDGOALS))

# Include the plugin dependencies if there are any targets specified that are
# not related to cleaning up files
ifeq ($(PLUGIN_NONCLEANGOALS),)
    # only maint targets, do nothing!
else
#   $(warning *** PLUGIN_NONCLEANGOALS="$(PLUGIN_NONCLEANGOALS)")
#   $(warning *** PWD="$(PWD)")
#   $(warning *** D_TARGETS="$(PLUGIN_OBJ:.o=.d)")
    -include $(PLUGIN_OBJ:.o=.d)
endif
