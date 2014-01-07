#
# For building Makefiles for projects in ATS 
# in a portable fashion.
#
# By portable, we mean non-dependent on an ATS installation;
# certain settings in the makefiles may need to be adjusted,
# but these should primarily be limited to this file.


######

ifndef ATSHOME
  ATSHOMEQ="$(ATSHOME)"
else
  ATSHOMEQ="$(ATSHOME)"
endif
ifndef ATSHOMERELOC
  ATSHOMERELOCQ="$(ATSHOMERELOC)"
else
  ATSHOMERELOCQ="$(ATSHOMERELOC)"
endif

######

ATSCC=$(ATSHOMEQ)/bin/atscc
ATSOPT=$(ATSHOMEQ)/bin/atsopt

#The linker should decide between these
ATSLIB=$(ATSHOMEQ)/ccomp/atslib/lib
ATSLIB64=$(ATSHOMEQ)/ccomp/atslib/lib64

######

CFLAGS += -D_GNU_SOURCE

######

ifndef ATSHOME
else
LDFLAGS += -L$(ATSLIB)
LDFLAGS += -L$(ATSLIB64)
LDFLAGS += -latslib
endif

######

MALLOCFLAG := -DATS_MEMALLOC_LIBC

######

ifndef ATSHOME
else
MYPORTDIR=$(ATSHOME)
endif

INCLUDE += -I$(strip $(MYPORTDIR))
INCLUDE_ATS += -IATS $(strip $(MYPORTDIR))
INCLUDE += -I$(strip $(MYPORTDIR))/contrib
INCLUDE_ATS += -IATS $(strip $(MYPORTDIR))/contrib
INCLUDE += -I$(strip $(MYPORTDIR))/ccomp/runtime
INCLUDE_ATS += -IATS $(strip $(MYPORTDIR))/ccomp/runtime

######

all::
cleanats::
cleanall::

###### end of [atsmake-pre.mk] ######

#Notes: 
# Not sure why strip is necessary sometimes for MYPORTDIR. e.g.
# INCLUDE += -I$(strip $(MYPORTDIR))/contrib

# Some rules may not work with new versions of make (search for 
# "3.82" in the makefiles.
