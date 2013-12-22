#
#
# Makefile for K&R examples in Chapter 1
#
#

######

ATSUSRQ="$(ATSHOME)"
ifeq ($(ATSUSRQ),"")
ATSUSRQ="/usr"
endif

ATSPREL=$(ATSHOME)/ccomp/runtime/ats_prelude.c
######

CC=gcc
ATSCC=$(ATSUSRQ)/bin/atscc
ATSOPT=$(ATSUSRQ)/bin/atsopt

######

CFLAGS_ATS=-I$(ATSHOME) -I$(ATSHOME)/ccomp/runtime/
LDFLAGS_ATS=-L$(ATSHOME)/ccomp/lib/
LDFLAGS = $(LDFLAGS_ATS)
CFLAGS = $(CFLAGS_ATS)
ATSGC = $(ATSHOME)/ccomp/runtime/GCATS/gc.o
######

.PHONY: all
all:: checkall
# all:: cleanall

######

checkall::
cleanall:: clean
cleanall:: ; $(RMF) *_?ats.html 

######


checkall:: minDisj
cleanall:: ; $(RMF) minDisj
minDisj: minDisjNoCov_dats.c sstream_dats.c sstream_sats.c
	$(CC) -O2 -D_ATS_GCATS $(CFLAGS) $(LDFLAGS) -o minDisj $(ATSGC) \
        minDisjNoCov_dats.c sstream_dats.c sstream_sats.c $(ATSPREL) -lm -lats

minDisjNoCov_dats.c: minDisjNoCov.dats
	$(ATSOPT) --output minDisjNoCov_dats.c --dynamic minDisjNoCov.dats
sstream_dats.c: sstream.dats
	$(ATSOPT) --output sstream_dats.c --dynamic sstream.dats
sstream_sats.c: sstream.sats
	$(ATSOPT) --output sstream_sats.c --static sstream.sats

######

html:: ; $(ATSOPT) --posmark_html -d minDisjNoCov.dats > minDisjNoCov_dats.html

######

RMF = rm -f

######

clean:
	$(RMF) *~
	$(RMF) *_?ats.c *_?ats.o

###### end of [Makefile] ######
