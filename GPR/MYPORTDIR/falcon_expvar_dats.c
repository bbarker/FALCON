/*
**
** The C code is generated by ATS/Postiats
** The compilation time is: 2014-4-17: 16h:22m
**
*/

/*
** include runtime header files
*/
#ifndef _ATS_CCOMP_HEADER_NONE
#include "pats_ccomp_config.h"
#include "pats_ccomp_basics.h"
#include "pats_ccomp_typedefs.h"
#include "pats_ccomp_instrset.h"
#include "pats_ccomp_memalloc.h"
#ifndef _ATS_CCOMP_EXCEPTION_NONE
#include "pats_ccomp_memalloca.h"
#include "pats_ccomp_exception.h"
#endif // end of [_ATS_CCOMP_EXCEPTION_NONE]
#endif /* _ATS_CCOMP_HEADER_NONE */


/*
** include prelude cats files
*/
#ifndef _ATS_CCOMP_PRELUDE_NONE
//
#include "prelude/CATS/basics.cats"
#include "prelude/CATS/integer.cats"
#include "prelude/CATS/pointer.cats"
#include "prelude/CATS/bool.cats"
#include "prelude/CATS/char.cats"
#include "prelude/CATS/integer_ptr.cats"
#include "prelude/CATS/integer_fixed.cats"
#include "prelude/CATS/float.cats"
#include "prelude/CATS/memory.cats"
#include "prelude/CATS/string.cats"
#include "prelude/CATS/strptr.cats"
//
#include "prelude/CATS/filebas.cats"
//
#include "prelude/CATS/list.cats"
#include "prelude/CATS/option.cats"
#include "prelude/CATS/array.cats"
#include "prelude/CATS/arrayptr.cats"
#include "prelude/CATS/arrayref.cats"
#include "prelude/CATS/matrix.cats"
#include "prelude/CATS/matrixptr.cats"
//
#endif /* _ATS_CCOMP_PRELUDE_NONE */
/*
** for user-supplied prelude
*/
#ifdef _ATS_CCOMP_PRELUDE_USER
//
#include _ATS_CCOMP_PRELUDE_USER
//
#endif /* _ATS_CCOMP_PRELUDE_USER */

/*
staload-prologues(beg)
*/
/*
staload-prologues(end)
*/
/*
typedefs-for-tyrecs-and-tysums(beg)
*/
typedef
struct {
atstkind_t0ype(atstype_double) atslab__0; 
atstkind_t0ype(atstype_double) atslab__1; 
} postiats_tyrec_0 ;
/*
typedefs-for-tyrecs-and-tysums(end)
*/
/*
dynconlst-declaration(beg)
*/
/*
dynconlst-declaration(end)
*/
/*
dyncstlst-declaration(beg)
*/
ATSdyncst_extfun(_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__fprint_expvar, (atstkind_type(atstype_ptrk), postiats_tyrec_0), atsvoid_t0ype) ;
ATSdyncst_mac(atspre_FILE_stdout)
ATSdyncst_mac(atspre_fprint_string)
ATSdyncst_mac(atspre_fprint_double)
/*
dyncstlst-declaration(end)
*/
/*
dynvalist-implementation(beg)
*/
/*
dynvalist-implementation(end)
*/
/*
exnconlst-declaration(beg)
*/
#ifndef _ATS_CCOMP_EXCEPTION_NONE
extern void the_atsexncon_initize (atstype_exncon *d2c, char *exnmsg) ;
#endif // end of [_ATS_CCOMP_EXCEPTION_NONE]
/*
exnconlst-declaration(end)
*/
/*
assumelst-declaration(beg)
*/
/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 93(line=11, offs=1) -- 131(line=12, offs=32)
*/
ATSassume(_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_type) ;
/*
assumelst-declaration(end)
*/
/*
extypelst-declaration(beg)
*/
/*
extypelst-declaration(end)
*/
#if(0)
ATSglobaldec()
postiats_tyrec_0
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_make(atstkind_t0ype(atstype_double), atstkind_t0ype(atstype_double)) ;
#endif // end of [QUALIFIED]

#if(0)
ATSglobaldec()
atstkind_t0ype(atstype_double)
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_exp(postiats_tyrec_0) ;
#endif // end of [QUALIFIED]

#if(0)
ATSglobaldec()
atstkind_t0ype(atstype_double)
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_var(postiats_tyrec_0) ;
#endif // end of [QUALIFIED]

#if(0)
ATSglobaldec()
atsvoid_t0ype
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__print_expvar(postiats_tyrec_0) ;
#endif // end of [QUALIFIED]

#if(0)
ATSglobaldec()
atsvoid_t0ype
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__fprint_expvar(atstkind_type(atstype_ptrk), postiats_tyrec_0) ;
#endif // end of [QUALIFIED]

/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 180(line=18, offs=3) -- 210(line=18, offs=33)
*/
/*
local: 
global: expvar_make$0$0(level=0)
local: 
global: 
*/
ATSglobaldec()
postiats_tyrec_0
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_make(atstkind_t0ype(atstype_double) arg0, atstkind_t0ype(atstype_double) arg1)
{
/* tmpvardeclst(beg) */
ATStmpdec(tmpret0, postiats_tyrec_0) ;
/* tmpvardeclst(end) */
/* funbodyinstrlst(beg) */
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 166(line=17, offs=1) -- 210(line=18, offs=33)
*/
__patsflab_expvar_make:
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 196(line=18, offs=19) -- 210(line=18, offs=33)
*/
ATSINSstore_fltrec_ofs (tmpret0, postiats_tyrec_0, atslab__0, arg0) ;
ATSINSstore_fltrec_ofs (tmpret0, postiats_tyrec_0, atslab__1, arg1) ;
/* funbodyinstrlst(end) */
ATSreturn(tmpret0) ;
} /* end of [_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_make] */

/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 260(line=23, offs=16) -- 269(line=23, offs=25)
*/
/*
local: 
global: expvar_get_exp$1$0(level=0)
local: 
global: 
*/
ATSglobaldec()
atstkind_t0ype(atstype_double)
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_exp(postiats_tyrec_0 arg0)
{
/* tmpvardeclst(beg) */
ATStmpdec(tmpret1, atstkind_t0ype(atstype_double)) ;
/* tmpvardeclst(end) */
/* funbodyinstrlst(beg) */
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 245(line=23, offs=1) -- 269(line=23, offs=25)
*/
__patsflab_expvar_get_exp:
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 266(line=23, offs=22) -- 269(line=23, offs=25)
*/
ATSINSmove(tmpret1, ATSSELfltrec(arg0, postiats_tyrec_0, atslab__0)) ;
/* funbodyinstrlst(end) */
ATSreturn(tmpret1) ;
} /* end of [_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_exp] */

/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 295(line=25, offs=16) -- 304(line=25, offs=25)
*/
/*
local: 
global: expvar_get_var$2$0(level=0)
local: 
global: 
*/
ATSglobaldec()
atstkind_t0ype(atstype_double)
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_var(postiats_tyrec_0 arg0)
{
/* tmpvardeclst(beg) */
ATStmpdec(tmpret2, atstkind_t0ype(atstype_double)) ;
/* tmpvardeclst(end) */
/* funbodyinstrlst(beg) */
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 280(line=25, offs=1) -- 304(line=25, offs=25)
*/
__patsflab_expvar_get_var:
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 301(line=25, offs=22) -- 304(line=25, offs=25)
*/
ATSINSmove(tmpret2, ATSSELfltrec(arg0, postiats_tyrec_0, atslab__1)) ;
/* funbodyinstrlst(end) */
ATSreturn(tmpret2) ;
} /* end of [_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__expvar_get_var] */

/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 352(line=30, offs=14) -- 387(line=30, offs=49)
*/
/*
local: 
global: print_expvar$3$0(level=0)
local: 
global: 
*/
ATSglobaldec()
atsvoid_t0ype
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__print_expvar(postiats_tyrec_0 arg0)
{
/* tmpvardeclst(beg) */
ATStmpdec_void(tmpret3, atsvoid_t0ype) ;
/* tmpvardeclst(end) */
/* funbodyinstrlst(beg) */
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 339(line=30, offs=1) -- 387(line=30, offs=49)
*/
__patsflab_print_expvar:
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 358(line=30, offs=20) -- 387(line=30, offs=49)
*/
ATSINSmove_void(tmpret3, _057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__fprint_expvar(atspre_FILE_stdout, arg0)) ;

/* funbodyinstrlst(end) */
ATSreturn_void(tmpret3) ;
} /* end of [_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__print_expvar] */

/*
/cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 415(line=33, offs=15) -- 475(line=34, offs=50)
*/
/*
local: 
global: fprint_expvar$4$0(level=0)
local: 
global: 
*/
ATSglobaldec()
atsvoid_t0ype
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__fprint_expvar(atstkind_type(atstype_ptrk) arg0, postiats_tyrec_0 arg1)
{
/* tmpvardeclst(beg) */
ATStmpdec_void(tmpret4, atsvoid_t0ype) ;
ATStmpdec_void(tmp5, atsvoid_t0ype) ;
ATStmpdec_void(tmp6, atsvoid_t0ype) ;
ATStmpdec_void(tmp7, atsvoid_t0ype) ;
ATStmpdec_void(tmp8, atsvoid_t0ype) ;
/* tmpvardeclst(end) */
/* funbodyinstrlst(beg) */
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 401(line=33, offs=1) -- 475(line=34, offs=50)
*/
__patsflab_fprint_expvar:
/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 428(line=34, offs=3) -- 475(line=34, offs=50)
*/
ATSINSmove_void(tmp5, atspre_fprint_string(arg0, ATSPMVstring("(exp="))) ;

/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 428(line=34, offs=3) -- 475(line=34, offs=50)
*/
ATSINSmove_void(tmp6, atspre_fprint_double(arg0, ATSSELfltrec(arg1, postiats_tyrec_0, atslab__0))) ;

/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 428(line=34, offs=3) -- 475(line=34, offs=50)
*/
ATSINSmove_void(tmp7, atspre_fprint_string(arg0, ATSPMVstring(", var="))) ;

/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 428(line=34, offs=3) -- 475(line=34, offs=50)
*/
ATSINSmove_void(tmp8, atspre_fprint_double(arg0, ATSSELfltrec(arg1, postiats_tyrec_0, atslab__1))) ;

/*
emit_instr: loc0 = /cygdrive/c/cygwin/home/brand_000/FBA/FALCON/GPR/falcon_expvar.dats: 428(line=34, offs=3) -- 475(line=34, offs=50)
*/
ATSINSmove_void(tmpret4, atspre_fprint_string(arg0, ATSPMVstring(")"))) ;

/* funbodyinstrlst(end) */
ATSreturn_void(tmpret4) ;
} /* end of [_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_056_sats__fprint_expvar] */

/*
** for initialization(dynloading)
*/
atsvoid_t0ype
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_expvar_056_dats__dynload()
{
ATSdynload1(
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_expvar_056_dats__dynloadflag
) ;
ATSif(
ATSCKiseqz(
_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_expvar_056_dats__dynloadflag
)
) ATSthen() {
ATSdynloadset(_057_cygdrive_057_c_057_cygwin_057_home_057_brand_000_057_FBA_057_FALCON_057_GPR_057_falcon_expvar_056_dats__dynloadflag) ;
/*
dynexnlst-initize(beg)
*/
/*
dynexnlst-initize(end)
*/
} /* ATSendif */
ATSreturn_void() ;
} /* end of [*_dynload] */