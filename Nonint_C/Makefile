# ####################################################
# 
# The make file for Non-integer sparse FFT code
# 
# Author: 	lixin
# Date: 	11/13/2012
#
# ###################################################

# compiler spec
# Note that for condor, it has to be g++-4.4
CXX = g++

# options
# I don't want to include pedantic, because this guy gives me warnings for complex imaginary
CXXFLAGS = -Wall -Wextra -Wredundant-decls -Wdisabled-optimization -Wnon-virtual-dtor -Woverloaded-virtual -Wsign-promo -Wno-long-long -Wno-unused-local-typedefs #-pedantic 

# I would like to add warning for failure of inline
# because the code forces a lot of inline functions
# remove it if this annoys you
# CXXFLAGS += -Winline

# gcc optimization options
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -frounding-math -O3 -funroll-loops

UPDATE_PINV_FLAG = -D UPDATING_PSEUDOINVERSE
NOTUPDATE_PINV_FLAG = -U UPDATING_PSEUDOINVERSE

# for gbd debugger
CXXDEBUGFLAGS = -g

# turn on/off debug
#MACRO_DEBUG = -D DEBUG_ON
MACRO_RELEASE = -U DEBUG_ON

CXXDEBUGFLAGS += $(MACRO_DEBUG)

CXXRELASEFLAGS += $(MACRO_RELEASE)
# I am really permissive for inline function
# CXXRELASEFLAGS += -finline-limit=100000 --param inline-unit-growth=5000 --param large-function-growth=1000000

# include directory
INCDIR = inc
# build directory
BUILDDIR = build
# source directory
SRCDIR = src
# document directory
DOCDIR = doc
# library directory
#LIBDIR = lib
# the object name
object = nonint_sfft
# the doxygen file
doxy = Doxyfile
# the eigen library directory
EIGENDIR = eigen
# dependency path
VPATH = $(INCDIR):$(SRCDIR)

# include directory setting
CXXFLAGS += -I$(INCDIR) -I$(EIGENDIR)

# linking flags
LDFLAGS = -L/usr/lib/

# fftw
LDFLAGS += -lfftw3 -lm -lfftw3f -lconfig++
#LDFLAGS += -lfftw3 -lfftw3f -lconfig++

createbuild = $(shell test -d $(BUILDDIR) || mkdir -p $(BUILDDIR))
createdoc = $(shell test -d $(DOCDIR) || mkdir -p $(DOCDIR))

#SRC_MAIN = main.cc

#INC_COMMON = common.h
#INC_PARAMS = params.h
#INC_PEAKS = peaks.h
#INC_PROJECTION = projection.h
#INC_UTILS = utils.h
#INC_INITIALIZATION = initialization.h
#INC_GRADIENT = gradient.h
#INC_SFFT = sfft.h
#INC_DATA = data_io.h
#INC_NONINT_RECOVER_GENERAL = nonint_recover_general.h
#INC_NONINT_RECOVER = nonint_recover.h
#SRC_COMMON = common.cc
#SRC_PEAKS = peaks.cc
#SRC_PARAMS = params.cc
#SRC_GRADIENT = gradient.cc
#SRC_UTILS = utils.cc

#INC_ALL = $(INC_COMMON) $(INC_PARAMS) $(INC_PEAKS) $(INC_PROJECTION) $(INC_UTILS) $(INC_INITIALIZATION) $(INC_GRADIENT) $(INC_NONINT_RECOVER_GENERAL) $(INC_NONINT_RECOVER) $(INC_SFFT) $(INC_DATA)
INC_ALL = $(wildcard $(INCDIR)/*.h)
#SRC_ALL = $(SRC_MAIN) $(SRC_COMMON) $(SRC_PARAMS) $(SRC_GRADIENT)
SRC_ALL = $(wildcard $(SRCDIR)/*.cc)

.PHONY: clean all doc cleandoc release

$(BUILDDIR)/$(object) : $(SRC_ALL) $(INC_ALL)
	$(createbuild)
	$(CXX) $(CXXFLAGS) $(CXXRELASEFLAGS) -o $@ $^ $(LDFLAGS) $(UPDATE_PINV_FLAG)

$(BUILDDIR)/$(object)_noupdate : $(SRC_ALL) $(INC_ALL)
	$(createbuild)
	$(CXX) $(CXXFLAGS) $(CXXRELASEFLAGS) -o $@ $^ $(LDFLAGS) $(NOTUPDATE_PINV_FLAG)

release : $(BUILDDIR)/$(object)

all : release

clean:
	rm -rf $(BUILDDIR)

doc: $(INC_ALL) $(SRC_ALL)
	$(createdoc)
	doxygen $(doxy)

cleandoc: 
	rm -rf $(DOCDIR)
