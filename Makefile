###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES         := no

CHEMICAL_REACTIONS  := no
CONTACT             := no
FLUID_PROPERTIES    := no
HEAT_CONDUCTION     := no
MISC                := no
NAVIER_STOKES       := no
PHASE_FIELD         := no
RDG                 := no
RICHARDS            := no
SOLID_MECHANICS     := no
STOCHASTIC_TOOLS    := no
TENSOR_MECHANICS    := yes
WATER_STEAM_EOS     := no
XFEM                := no
POROUS_FLOW         := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
SUPPORT_DIR        := /projects/nathanm/micromorphic/micromorphic_library/micromorphic_element/src/cpp
APPLICATION_NAME   := tardigrade
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk

ex_srcfiles := $(shell find $(APPLICATION_DIR) -name "*.C")
ex_deps     := $(patsubst %.C, %.$(obj-suffix).d, $(ex_srcfiles))

###############################################################################
# Additional special case targets should be added here

ADDITIONAL_CPPFLAGS := -lmicromat -lmicrobalance -fmax-errors=5
ADDITIONAL_CPPFLAGS += -I$(SUPPORT_DIR)
EXTERNAL_FLAGS      := -L$(SUPPORT_DIR) -lmicromat -lmicrobalance
