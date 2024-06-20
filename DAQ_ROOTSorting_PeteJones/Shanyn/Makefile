# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

#------------------------------------------------------------------------------

CLOVERTREE2O       = labrsort4.$(ObjSuf)
CLOVERTREE2S       = labrsort4.$(SrcSuf)
CLOVERTREE2        = labrsort4$(ExeSuf)

#EVENTTREEO       = eventtree.$(ObjSuf)
#EVENTTREES       = eventtree.$(SrcSuf)
#EVENTTREE        = eventtree$(ExeSuf)

#INDEXO       = index.$(ObjSuf)
#INDEXS       = index.$(SrcSuf)
#INDEX        = index$(ExeSuf)

#ANAL3O       = anal3.$(ObjSuf)
#ANAL3S       = anal3.$(SrcSuf)
#ANAL3        = anal3$(ExeSuf)

#KODEO       = kode.$(ObjSuf)
#KODES       = kode.$(SrcSuf)
#KODE        = kode$(ExeSuf)

TREESTUFFO       = treestuff.$(ObjSuf)
TREESTUFFS       = treestuff.$(SrcSuf)
TREESTUFF        = treestuff$(ExeSuf)

OBJS          = $(EVENTO) $(MAINEVENTO) $(CLOVERTREE2O) $(EVENTTREEO) $(INDEXO) \
                $(TSTRINGO) $(TCOLLEXO) $(VVECTORO) $(VMATRIXO) $(VLAZYO) \
                $(HELLOO) $(ACLOCKO) $(STRESSO) $(TBENCHO) $(BENCHO) \
                $(TCOLLBMO) $(STRESSGEOMO) $(STRESSLO) $(TESTBITSO) \
                $(CTORTUREO) $(QPRANDOMO) $(THREADO) $(TSjumpO) $(treestuffO) $(KODEO) $(ANAL3O)

PROGRAMS      = $(EVENT) $(CLOVERTREE2) $(EVENTTREE) $(INDEX) $(TSTRING) \
                $(TCOLLEX) $(TCOLLBM) $(VVECTOR) $(VMATRIX) $(VLAZY) \
                $(HELLOSO) $(ACLOCKSO) $(STRESS) $(TBENCHSO) $(BENCH) \
                $(STRESSGEOM) $(STRESSL) $(TESTBITS) $(CTORTURE) $(QPRANDOM) \
                $(THREADS) $(TSjump) $(treestuff) $(KODE) $(ANAL3)

OBJS         += $(GUITESTO) $(GUIVIEWERO) $(TETRISO)
PROGRAMS     += $(GUITEST) $(GUIVIEWER) $(TETRISSO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

$(EVENTSO):     $(EVENTO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $(EVENTO) $(OutPutOpt) $(EVENTSO)
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"

$(EVENT):       $(EVENTSO) $(MAINEVENTO)
		$(LD) $(LDFLAGS) $(MAINEVENTO) $(EVENTLIB) $(LIBS) \
		   $(OutPutOpt)$(EVENT)
		@echo "$@ done"

$(CLOVERTREE2):      $(CLOVERTREE2O)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EVENTTREE):      $(EVENTTREEO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(INDEX):    $(INDEXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(KODE):    $(KODEO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(ANAL3):    $(ANAL3O)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(TREESTUFF):    $(TREESTUFFO)
		$(LD) $(LDFLAGS) -lTreePlayer $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(HSIMPLE):     $(HSIMPLEO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(MINEXAM):     $(MINEXAMO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(TSTRING):     $(TSTRINGO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(TCOLLEX):     $(TCOLLEXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(TCOLLBM):     $(TCOLLBMO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(VVECTOR):     $(VVECTORO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(VMATRIX):     $(VMATRIXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(VLAZY):       $(VLAZYO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(STRESSL):     $(STRESSLO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(TESTBITS):    $(TESTBITSO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(THREADS):     $(THREADSO)
ifeq ($(HASTHREAD),yes)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) ../lib/libThread.lib $(OutPutOpt)$@
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lThread $(OutPutOpt)$@
endif
		@echo "$@ done"
else
		@echo "This version of ROOT has no thread support, $@ not built"
endif

$(QPRANDOM):    $(QPRANDOMO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) ../lib/libQuadp.lib $(OutPutOpt)$@
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lQuadp $(OutPutOpt)$@
endif
		@echo "$@ done"

$(GUITEST):     $(GUITESTO)
		$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(GUIVIEWER):   $(GUIVIEWERO)
		$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(STRESS):      $(STRESSO) $(EVENT)
		$(LD) $(LDFLAGS) $(STRESSO) $(EVENTLIB) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(STRESSGEOM):  $(STRESSGEOMO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) ../lib/libGeom.lib $(OutPutOpt)$@
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lGeom $(OutPutOpt)$@
endif
		@echo "$@ done"

$(BENCH):       $(BENCHO) $(TBENCHSO)
		$(LD) $(LDFLAGS) $(BENCHO) $(TBENCHO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

Hello:          $(HELLOSO)
$(HELLOSO):     $(HELLOO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
endif
endif
endif

Aclock:         $(ACLOCKSO)
$(ACLOCKSO):    $(ACLOCKO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
endif
endif
endif

Tetris:         $(TETRISSO)
$(TETRISSO):    $(TETRISO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(GLIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(GLIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 	   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
endif
endif
endif

$(TBENCHSO):    $(TBENCHO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ -lstdc++ $(OutPutOpt) $@
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
endif
endif
endif

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###

Event.$(ObjSuf): Event.h
MainEvent.$(ObjSuf): Event.h

EventDict.$(SrcSuf): Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

Hello.$(ObjSuf): Hello.h
HelloDict.$(SrcSuf): Hello.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

Aclock.$(ObjSuf): Aclock.h
AclockDict.$(SrcSuf): Aclock.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

Tetris.$(ObjSuf): Tetris.h
TetrisDict.$(SrcSuf): Tetris.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

TBench.$(ObjSuf): TBench.h
TBenchDict.$(SrcSuf): TBench.h benchLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

guiviewer.$(ObjSuf): guiviewer.h
guiviewerDict.$(SrcSuf): guiviewer.h guiviewerLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
