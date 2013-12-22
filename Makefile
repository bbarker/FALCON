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

######

ATSCC=$(ATSUSRQ)/bin/atscc
ATSOPT=$(ATSUSRQ)/bin/atsopt

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
minDisj: minDisjNoCov.dats sstream.dats sstream.sats
	$(ATSCC) -O2 -D_ATS_GCATS -o minDisj minDisjNoCov.dats sstream.dats sstream.sats

######

html:: ; $(ATSOPT) --posmark_html -d minDisjNoCov.dats > minDisjNoCov_dats.html

######

RMF = rm -f

######

clean:
	$(RMF) *~
	$(RMF) *_?ats.c *_?ats.o

###### end of [Makefile] ######
