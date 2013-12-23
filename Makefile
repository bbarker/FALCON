
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
# Need to remove this and generate all source:
ATSGC = $(ATSHOME)/ccomp/runtime/GCATS/gc.o
######

.PHONY: all
all:: checkall
# all:: cleanall

######

checkall::
cleanall:: clean
cleanall:: ; $(RMF) *_?ats.html 
cleanall:: ; $(RMF) *.c
######


checkall:: minDisj
cleanall:: ; $(RMF) minDisj
minDisj: minDisjNoCov_dats.c sstream_dats.c sstream_sats.c
	$(CC) -O2 -D_ATS_GCATS $(CFLAGS) $(LDFLAGS) -o minDisj $(ATSGC) \
        minDisjNoCov_dats.c sstream_dats.c sstream_sats.c $(ATSPREL) -lm -lats

%_sats.c: %.sats; $(ATSOPT) --output $@ --static  $<
%_dats.c: %.dats; $(ATSOPT) --output $@ --dynamic $<

######

html:: ; $(ATSOPT) --posmark_html -d minDisjNoCov.dats > minDisjNoCov_dats.html

######

RMF = rm -f

######

clean:
	$(RMF) *~
	$(RMF) *_?ats.o

###### end of [Makefile] ######
