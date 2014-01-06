#
# For building Makefiles for portable projects in ATS 
#

######

ifndef ATSHOME
  ATSHOMEQ="$(ATSHOME)"
else
  ATSHOMEQ="$(ATSHOME)"
endif

#ATSHOMERELOC - should be set to ATSHOME for git repo?
ifndef ATSHOMERELOC
  ATSHOMERELOCQ="$(ATSHOME)"
else
  ATSHOMERELOCQ="$(ATSHOMERELOC)"
endif

######

ATSCC=$(ATSHOMEQ)/bin/atscc
ATSOPT=$(ATSHOMEQ)/bin/atsopt
ATSLIB=$(ATSHOMEQ)/ccomp/atslib/lib
ATSLIB64=$(ATSHOMEQ)/ccomp/atslib/lib64

######

CFLAGS += -D_GNU_SOURCE -std=c99 -D_XOPEN_SOURCE

######

LDFLAGS += -Xlinker --allow-multiple-definition 
LDFLAGS += -L$(ATSLIB) -L$(ATSLIB64) -latslib

######

MALLOCFLAG := -D_ATS_GCATS

######

ifndef ATSHOME
MYPORTDIR=$(ATSDEPDIR)
else
MYPORTDIR=$(ATSHOME)
endif

INCLUDE_ATS += -IATS $(MYPORTDIR)/contrib
INCLUDE_ATS_C := -I$(MYPORTDIR) -I$(MYPORTDIR)/ccomp/runtime 
INCLUDE_ATS_PC := -I$(MYPORTDIR) -I$(MYPORTDIR)/ccomp/runtime

######

all::
cleanats::
cleanall::

###### end of [atsmake-pre.mk] ######
