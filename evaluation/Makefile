#@file    Makefile
#@brief   Makefile for C-project "Evaluation for delay-robust VRPTW"
#@author  Andreas M. Tillmann, TU Braunschweig, Lukas Schürmann University Bonn

# uses SCIP not because it is truly needed but to facilitate incorporating the
# functions into the MIP solver later (and to exploit SCIP's memory management). 


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       $(SCIP_PATH)

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)


#-----------------------------------------------------------------------------
# default settings for test runs
#-----------------------------------------------------------------------------
#VERSION		=	0.1
#LOCK		=	false
#VALGRIND	=	false
#TEST    	= 	testset
#SETTINGS    	=   	default
#USEBESTSOL	=	true
#TIME		=	3600
#NODE		=	2100000000
#MEM 		=	12000
#DISPFREQ	=	10000
#THREADS		=	1
#PERMUTE         =       0
#SEEDS           =       0
#GLBSEEDSHIFT    =       0
#FEASTOL         =       default
#CONTINUE        =       false
#LOCK            =       false
#DEBUGTOOL       =       none
#CLIENTTMPDIR    =       /tmp
#REOPT           =       false
#OPTCOMMAND      =       optimize
#SETCUTOFF       =       0
#VISUALIZE       =       false
#MAXJOBS         =       1


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------
MAINNAME	=	evaluator
MAINOBJ		=	cmain.o tools_evaluating.o tools_data.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------
ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -i $(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -I/usr/include -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^


.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h;

.PHONY: test
test:
	cd check; \
	$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(DISPFREQ) \
	$(LOCK) $(VERSION) $(LPS) $(USEBESTSOL) $(VALGRIND);

.PHONY: testcluster
testcluster:
	cd check; \
	$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) \
	$(DISPFREQ) $(CONTINUE) $(VERSION) $(LPS) $(USEBESTSOL) $(QUEUE) $(PPN) $(CLIENTTMPDIR) $(NOWAITCLUSTER);


#-----------------------------------------------------------------------------
# Linking and Compiling
#-----------------------------------------------------------------------------
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
ifeq ($(VERBOSE),true)
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
else
		@$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
