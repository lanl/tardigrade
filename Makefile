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
MISC                := yes
NAVIER_STOKES       := no
PHASE_FIELD         := no
RDG                 := no
RICHARDS            := no
SOLID_MECHANICS     := no
STOCHASTIC_TOOLS    := no
TENSOR_MECHANICS    := yes
XFEM                := no
POROUS_FLOW         := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
MICROMORPHIC_DIR   := /home/nathan/research
MICROMORPHIC_COMPILER_PATH := /usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu/
ANACONDA_INCLUDE   := /home/nathan/anaconda3/include
ANACONDA_LIB       := /home/nathan/anaconda3/lib
APPLICATION_NAME   := tardigrade
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here

#ADDITIONAL_INCLUDES += -I$(MICROMORPHIC_DIR)/voro++ -I$(MICROMORPHIC_DIR)/micromorphic_element/src/cpp
#ADDITIONAL_INCLUDES += -I$(MICROMORPHIC_DIR)/vector_tools/src/cpp
#ADDITIONAL_INCLUDES += -I$(MICROMORPHIC_DIR)/overlap -I$(MICROMORPHIC_DIR)/quickhull

ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/voro++/voro++ -I$(MICROMORPHIC_DIR)/micromorphic_element/src/cpp
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/vector_tools/src/cpp -I$(MICROMORPHIC_DIR)/error_tools/src/cpp
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/overlap_coupling/src/cpp -I$(MICROMORPHIC_DIR)
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/solver_tools/src/cpp
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/xdmf -I$(MICROMORPHIC_DIR)/xdmf/build
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/xdmf/core -I$(MICROMORPHIC_DIR)/xdmf/build/core
ADDITIONAL_CPPFLAGS += -I$(ANACONDA_INCLUDE)/libxml2
ADDITIONAL_CPPFLAGS += -I$(BOOST_ROOT)
ADDITIONAL_CPPFLAGS += -I$(MICROMORPHIC_DIR)/yaml-cpp/include

ADDITIONAL_LIBS     += -L$(MICROMORPHIC_DIR)/voro++/voro++ -L$(MICROMORPHIC_DIR)/micromorphic_element/src/cpp
ADDITIONAL_LIBS     += -L$(MICROMORPHIC_DIR)/overlap_coupling/src/cpp
ADDITIONAL_LIBS     += -L$(MICROMORPHIC_DIR)/solver_tools/src/cpp
ADDITIONAL_LIBS     += -lmicromat -lmicrobalance -lvoro++# -loverlap
ADDITIONAL_LIBS     += -L$(MICROMORPHIC_COMPILER_PATH) -lresolv -DDEBUG_MODE
ADDITIONAL_LIBS     += -L$(MICROMORPHIC_DIR)/xdmf/build/lib -lXdmf -lXdmfCore -lXdmfUtils
ADDITIONAL_LIBS     += -L$(ANACONDA_LIB)/libxml2 -lxml2
ADDITIONAL_LIBS     += -L$(ANACONDA_LIB) -lhdf5 -ltiff
ADDITIONAL_LIBS     += -L$(MICROMORPHIC_DIR)/yaml=cpp/build

ifeq ($(METHOD),dbg)
    ADDITIONAL_CPPFLAGS += -DDEBUG_MODE
    ADDITIONAL_LIBS += -DDEBUG_MODE
endif
