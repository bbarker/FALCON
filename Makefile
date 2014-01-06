######

#
# A minimal makefile for generating portable C code
#

######       Environment Configuration      ######

CC = gcc
CCDEPFLAG = -c -MM -MF
ATSDEPDIR = ATSDEPCOPIES

CFLAGS += -O2
#LDFLAGS += -lm 

######   End of Environment Configuration   ######

######        Project Configuration         ######

SOURCES_DATS += minDisjNoCov.dats sstream.dats

SOURCES_SATS += sstream.sats

######

MYTARGET=mindisj

######    End of Project Configuration      ######

.PRECIOUS: *_?ats.c $(ATSDEPDIR)

all:: $(ATSDEPDIR)

$(ATSDEPDIR):
	mkdir $(ATSDEPDIR)	

######

include utils/atsmake-port-pre.mk   # Mostly environment setup
include utils/atsmake-port-post.mk  # Mostly build rules

###### end of [Makefile] ######
