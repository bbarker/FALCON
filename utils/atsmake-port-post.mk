#
# For building Makefiles for portable projects in ATS 
#

######

ifdef ATSHOME
ifndef MYTARGET 
else
SATS_C := $(patsubst %.sats, %_sats.c, $(SOURCES_SATS))
DATS_C := $(patsubst %.dats, %_dats.c, $(SOURCES_DATS))

%_sats.c: %.sats 
	$(ATSCC) -cc $(INCLUDE_ATS) $<
#	$(ATSOPT) $(INCLUDE_ATS) -o $@ -s $<
%_dats.c: %.dats 
	$(ATSCC) -cc $(MALLOCFLAG) $(INCLUDE_ATS) $<
#	$(ATSOPT) $(INCLUDE_ATS) $(MALLOCFLAG) -o $@ -d $<
endif
endif

ifndef MYTARGET
else
SATS_O := $(patsubst %.sats, %_sats.o, $(SOURCES_SATS))
DATS_O := $(patsubst %.dats, %_dats.o, $(SOURCES_DATS))

%_sats.o: %_sats.c
	$(CC) $(INCLUDE_ATS_PC) $(CFLAGS) $< -o $@
%_dats.o: %_dats.c
	$(warning "TESTING .o") 
	$(CC) $(INCLUDE_ATS_PC) $(MALLOCFLAG) $(CFLAGS) $< -o $@
endif

######
#
-include .depend
ifdef ATSHOME
#
depend:: ; $(RMF) -f .depend
#
ifndef SOURCES_SATS
else
depend:: ; $(ATSOPT) --output-a .depend --depgen -s $(SOURCES_SATS)
endif
#
ifndef SOURCES_DATS
else
depend:: ; $(ATSOPT) --output-a .depend --depgen -d $(SOURCES_DATS)
endif
endif
#
######

ifdef ATSHOME
portdepINIT: $(SATS_C) $(DATS_C)
#	$(eval INCLUDE_ATS_C := -I$(ATSHOME) -I$(ATSHOME)/ccomp/runtime) 
#	$(eval INCLUDE_ATS_PC := -I$(ATSDEPDIR) -I$(ATSDEPDIR)/ccomp/runtime)
	$(RMF) .atscdep .atscdep_tmp
portdepMKDEP: portdepINIT
#	$(CC) $(INCLUDE_ATS_C) $(MALLOCFLAG) $(CFLAGS) *_dats.c \ 
#          $(CCDEPFLAG) .atscdep_tmp; \
#	$(ATSCC) $(INCLUDE_ATS) $(MALLOCFLAG) $(CFLAGS) -c $< \
#          $(CCDEPFLAG) .atscdep_tmp
	$(ATSCC) $(INCLUDE_ATS) $(MALLOCFLAG) $(CFLAGS) -c *_dats.c \
          $(CCDEPFLAG) .atscdep_tmp; \
	tail -n +2 .atscdep_tmp >> .atscdep
portdepDEPSLIST: portdepMKDEP
	$(eval CATSDEPS := $(shell cat .atscdep | tr -d '\\\n'))
portdepDEPSSORT: portdepDEPSLIST
	$(eval CATSDEPS := $(sort $(CATSDEPS)))
portdepFROMDIRS: portdepDEPSSORT
	$(eval FROMDIRS := $(foreach cdep, $(CATSDEPS), $(shell dirname $(cdep))))

# Create necessary directory structure in local directory.
# Add created directories to INCLUDE_ATS_PC.
portdepTO: portdepFROMDIRS
	$(eval CATSDEPSTODIRS := $(subst $(ATSHOME)/, , $(FROMDIRS)))
	$(eval CATSDEPSTOFILES := $(subst $(ATSHOME)/, , $(CATSDEPS)))
portdepMKDIR: portdepTO
	$(foreach todir, $(CATSDEPSTODIRS), \
          $(shell mkdir -p $(ATSDEPDIR)/$(todir)))
	$(eval INCLUDE_ATS_PC := $(INCLUDE_ATS_PC) \
          $(foreach todir, $(CATSDEPSTODIRS), \
          -I./$(ATSDEPDIR)/$(todir)))
portdepCP: portdepMKDIR
	$(foreach tofil, $(CATSDEPSTOFILES), \
	$(shell cp $(ATSHOME)/$(tofil) $(ATSDEPDIR)/$(tofil)))
else
.PHONY: portdepCP
endif

######

#
#object intermediates: not working.
#
# all:: portdepCP $(SATS_O) $(DATS_O)
# 	$(CC) $(LDFLAGS) $(MALLOCFLAG) $(CFLAGS) $(INCLUDE_ATS_PC) *ats.o \
#           -o $(MYTARGET)

#
#skiping object intermediates
#

all:: portdepCP $(SATS_C) $(DATS_C)
	$(CC) $(MALLOCFLAG) $(CFLAGS) $(INCLUDE_ATS_PC) *ats.c \
          -o $(MYTARGET)

######

RMF=rm -f
RMFR=rm -fr

######

cleanport:: ; $(RMF) *~
cleanport:: ; $(RMF) *_?ats.o
cleanport:: ; $(RMF) .atscdep_tmp
cleanport:: ; $(RMF) *.exe $(MYTARGET)

######

cleanats:: ; $(RMF) .portdepTO .portdepMKDIR
cleanats:: ; $(RMF) *~
cleanats:: ; $(RMF) *_?ats.c
cleanats:: ; $(RMF) .atscdep .atscdep_tmp
cleanats:: ; $(RMF) atsmake-port-pre.mk atsmake-port-post.mk
cleanats:: ; $(RMFR) $(ATSDEPDIR) 

######

clean: cleanport

######

cleanall:: cleanats cleanport
cleanall:: ; $(RMF) .depend

###### end of [atsmake-post.mk] ######
