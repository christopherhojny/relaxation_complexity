#@file    Makefile
#@brief   Makefile for computing the relaxation complexity
#@author  Christopher Hojny

#-----------------------------------------------------------------------------
# path
#-----------------------------------------------------------------------------

SCIPDIR         =       $(SCIP_PATH)

#----------------------------------------------------------------------------
# include default project Makefile from SCIP (need to do this twice, once to
# find the correct binary, then, after getting the correct flags from the
# binary (which is necessary since the ZIMPL flags differ from the default
# if compiled with the SCIP Optsuite instead of SCIP), we need to set the
# compile flags, e.g., for the ZIMPL library, which is again done in make.project
#----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)
SCIPVERSION			:=$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override ARCH		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ARCH=\([^@]*\).*/\1/')
override EXPRINT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* EXPRINT=\([^@]*\).*/\1/')
override GAMS		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GAMS=\([^@]*\).*/\1/')
override GMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPSCHECK	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT 	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override PARASCIP	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* PARASCIP=\([^@]*\).*/\1/')
override READLINE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override ZIMPL		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
include $(SCIPDIR)/make/make.project
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

LIBPATH		=	lib
USRCFLAGS	+=      -isystem$(LIBPATH)
USRCXXFLAGS	+=      -isystem$(LIBPATH)

OPT		=	opt
LPS		=	spx2
VERSION		=	1.0
CONTINUE        =       false
LOCK            =       false
TEST            =       simple
SETTINGS        =       default
TIME            =       3600
NODES           =       2100000000
MEM             =       20000
DISPFREQ        =       10000
PERMUTE         =       0
ONLYPRE         =       false
SETCUTOFF       =       false

#-----------------------------------------------------------------------------
# Programs
#-----------------------------------------------------------------------------

COMPUTERC	=	computeRC

COMPUTERCOBJ	=	auxiliary_cdd.o \
			branch_ryanfoster.o \
			compute_symmetries.o \
			cons_samediff.o \
			cons_conflict.o \
			convex_hull.o \
			hiding_sets.o \
			maximal_separation.o \
			pricer_pattern.o \
			probdata_rc_cg.o \
			probdata_rc_compact.o \
			probdata_rc_conflict.o \
			problem_rc.o \
			prop_convexity.o \
			prop_intersection.o \
			rcParams.o \
			rcPlugins.o \
			vardata_binpacking.o \
			vardata_compact.o \
			main.o \
			readInstance.o

COMPUTERCSRC	=	auxiliary_cdd.c \
			branch_ryanfoster.c \
			compute_symmetries.cpp \
			cons_samediff.c \
			convex_hull.c \
			hiding_sets.c \
			maximal_separation.c \
			pricer_pattern.c \
			probdata_rc_cg.c \
			probdata_rc_compact.c \
			problem_rc.c \
			prop_convexity.c \
			rcParams.c \
			rcPlugins.c \
			vardata_binpacking.c \
			vardata_compact.c \
			main.cpp \
			readInstance.cpp

COMPUTERCOBJFILES =	$(addprefix $(OBJDIR)/,$(COMPUTERCOBJ))
COMPUTERCSRCFILES =	$(addprefix $(SRCDIR)/,$(COMPUTERCSRC))
COMPUTERCDEP	  =	$(SRCDIR)/depend.$(COMPUTERC)
COMPUTERCFILE	  =	$(BINDIR)/$(COMPUTERC).$(BASE).$(LPS)$(EXEEXTENSION)

LIBCDD		=	-L$(CDD_LIB_PATH) -lcdd
LIBBLISS	=	-L$(BLISS_LIB_PATH) -lbliss


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(COMPUTERCFILE) $(COMPUTERCOBJFILES)
endif

.PHONY: all
all:            $(LIBDIR) $(SCIPDIR) $(COMPUTERCFILE)



# ----------- others -------------
.PHONY: doc
doc:
		cd doc; $(DOXY) computeRC.dxy

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(LIBDIR):
		@-mkdir -p lib;
		@-ln -s $(CDD_PATH) lib/cddlib;
		@-ln -s $(BLISS_PATH) lib/bliss;

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(COMPUTERCSRCFILES) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(COMPUTERCDEP)'

# include dependencies
-include	$(COMPUTERC)



# link targets
$(COMPUTERCFILE): $(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(COMPUTERCOBJFILES)
		@echo "-> linking $@"
ifdef LINKCCSCIPALL
		-$(CXX) $(COMPUTERCOBJFILES) $(LINKCCSCIPALL) $(LIBCDD) $(LIBBLISS) -o $@
else
		-$(CXX) $(COMPUTERCOBJFILES) -L$(SCIPDIR)/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) -l$(NLPILIB) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS)  $(LIBCDD) $(LIBBLISS) -o $@
endif


# ----------------------------------------------------------------------------
$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(CC_c) $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@

# --- EOF ---------------------------------------------------------------------
