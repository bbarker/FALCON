(*                                                          *               
 *  Author: Brandon Barker                                  *
 *  February, 2013                                          * 
 *  Partly adapted from MATLAB code by                      *
 *  Kieran Smallbone and C# code by http://csliu.com        *
 *  Also, the parser was based on ML code from              *
 *  Scott Smith: Scott Smith                                *
 *  http://www.cs.jhu.edu/~scott/pl/lectures/parsing.html   * 
 *                                                          *
 *                                                          *
 *  Thanks to Hongwei Xi and Chris Double for their         *
 *  help in learning and using ATS and debugging            *
 *  this code, and to Narayanan Sadagopan for helping to    *
 *  test and debug this code.                               *)

// Compile:
// atscc -O2 -D_ATS_GCATS -o minDisj minDisj.dats sstream.dats sstream.sats -lm
 
(* TODO:

  In ATS2 port, consider using Hongwei's cstream.

  Add fn* to mutually recursive functions in parser and link by "and"s to
  make sure they are recognized as being mutually tail recursive.

  Although it is less likely that genes within a disjunction
  are correlation as much as genes within a conjunction,
  assuming this may be a problem; we should check actual 
  covariance values and take this into account; in this case
  we'd need to map from a set of two genes to a double value:
  one possibility is to sort and append gene names so we can
  use the same, faster data structure gDMap.

  Consider adding a rule for pleiotropy to prevent
  divvying up of expression among complexes that would
  have zero expression anyway. This would seeminly 
  require two passes through the algorithm, because
  we must first find the minimum disjunction before
  we can ascertain whether or not it is already 
  unexpressed altogether.

*)

(*
** please do not change unless you know what you do
*)

(*
//
staload _(*anon*) = "libc/SATS/stdio.sats"
//
staload _(*anon*) = "prelude/DATS/array.dats"
staload _(*anon*) = "prelude/DATS/array0.dats"
//
staload _(*anon*) = "prelude/DATS/list.dats"
staload _(*anon*) = "prelude/DATS/list0.dats"
staload _(*anon*) = "prelude/DATS/list_vt.dats"
//
staload _(*anon*) = "prelude/DATS/matrix.dats"
staload _(*anon*) = "prelude/DATS/matrix0.dats"
//
staload _(*anon*) = "prelude/DATS/option.dats"
staload _(*anon*) = "prelude/DATS/option0.dats"
//
staload _(*anon*) = "prelude/DATS/pointer.dats"
//
staload _(*anon*) = "prelude/DATS/reference.dats"
//

*)

(* ****** ****** *)
//
#include
"share/atspre_staload.hats"
//
(* ****** ****** *)

(* ********************************* Begin CODE ******************************** *)



staload
UN = "prelude/SATS/unsafe.sats"

staload "prelude/SATS/string.sats"
staload "prelude/SATS/filebas.sats"
// Not in ATS2:
//staload "prelude/SATS/printf.sats"
staload "libc/SATS/stdio.sats"
staload "prelude/SATS/integer.sats"
staload "prelude/SATS/list.sats"

staload  "libc/SATS/math.sats"
//dynload  "libc/SATS/math.sats"
staload _ = "libc/SATS/math.sats"

//staload "libats/SATS/linset_avltree.sats"   
//staload _ = "libats/DATS/linset_avltree.dats"   

//dynload "libats/DATS/linset_listord.dats"   

staload "sstream.sats"
dynload "sstream.dats"

exception UnforseenLexeme of ()
exception ParserParensUnbalanced of ()
exception ParserIllegalToken of ()
exception ParserPrematureEND of ()
exception InvalidData of ()
exception InvalidCase of ()
exception EmptyList of ()
exception MapKeyNotFound of ()
exception TestCase of ()

#define NAN 0.0/0.0

absviewtype genes

extern
fun genes_choose(xs: !genes, gene: &string? >> opt (string, b)
):<!wrt> #[b:bool] bool (b)
 
extern
fun genes_copy (xs: !genes): genes

extern
fun genes_free (xs: genes): void

extern 
fun genes_listize(xs: !genes): List(string)
#define NUL '\000'

extern 
fun genes_make_sing(x: string): genes

extern
fun genes_size(xs: !genes):<> size_t

extern 
fun genes_union (xs: genes, ys: genes): genes
overload + with genes_union


local 

staload LS = "libats/ATS1/SATS/linset_listord.sats"   
staload _ = "libats/ATS1/DATS/linset_listord.dats"   
assume genes = $LS.set(string)

in  // in of [local]

implement
genes_union(xs,ys) = $LS.linset_union(xs,ys, lam(x,y) 
  => compare_string_string(x,y))

implement
genes_choose(xs, gene) = $LS.linset_choose(xs, gene)

implement
genes_copy(xs) = $LS.linset_copy(xs)

implement
genes_free(xs) = $LS.linset_free(xs)

implement
genes_listize(xs) = list_of_list_vt($LS.linset_listize1 (xs))

implement
genes_make_sing(x) = $LS.linset_make_sing(x)

implement
genes_size(xs) = $LS.linset_size(xs)

end // end of local

absviewtype gDMap

extern
fun gDMap_find(mp: &gDMap, k: string): double

extern
fun gDMap_free(mp: gDMap):<!wrt> void

extern
fun gDMap_insert(mp: &gDMap, gene: string, dval: double
):<!wrt> bool

extern
fun gDMap_make_nil(): gDMap

local

staload LM = "libats/SATS/linmap_avltree.sats"   
staload _ = "libats/DATS/linmap_avltree.dats"   
assume gDMap = $LM.map(string,double)

in // in of [local]

implement 
gDMap_make_nil() = $LM.linmap_make_nil {string, double} () 

implement
gDMap_find (mp, k): double = let
  var res: double?
  val b = $LM.linmap_search (mp, k, res)
  in 
    if b then let
      prval () = opt_unsome {double} (res)
      in res end 
    else let 
      prval () = opt_unnone {double} (res) 
      in ( (* print(k + "not found\n"); $raise MapKeyNotFound;*) NAN) end
  end // end of [gDMap_find]

implement
gDMap_free(mp) = $LM.linmap_free(mp)

implement
gDMap_insert(mp,gene,dval) = let
  var rdval: double? 
  val b = $LM.linmap_insert<string,double>
    (mp, gene, dval, rdval)
  //How does this work?:
  prval () = opt_clear (rdval)
  in b end

end // end of [local]



(* These abstract view types are used to keep track
   of pleitropic reactions; we keep track of each
   reaction (line) where a given gene appears by adding
   it to that gene's line# set. Expression is later divided up
   equally among all reactions in which the gene appears.        *)

absviewtype gIntSetMap
absviewtype IntSet

extern
fun gIntSetMap_find(mp: &gIntSetMap, k: string): IntSet

extern
fun gIntSetMap_free(mp: gIntSetMap):<> void

extern
fun gIntSetMap_insert(mp: &gIntSetMap, gene: string, dval: double
):<> bool

extern
fun gIntSetMapp_make_nil(): gIntSetMap

local

staload LS = "libats/SATS/linset_listord.sats"   
staload _ = "libats/DATS/linset_listord.dats"   
staload LM = "libats/SATS/linmap_avltree.sats"   
staload _ = "libats/DATS/linmap_avltree.dats"   
assume IntSet = $LS.set(int) //constrain to be >=1 ?
assume gIntSetMap = $LM.map(string, IntSet)

in

(*
implement
gIntSetMap_find (mp, k): IntSet = let
  var res: IntSet?
  val b = $LM.linmap_search (mp, k, lam(x,y) 
    => compare_string_string(x,y), res)
  in 
    if b then let
      prval () = opt_unsome {IntSet} (res)
      in res end 
    else let 
      prval () = opt_unnone {IntSet} (res) 
      in (print(k + "not found\n"); $raise MapKeyNotFound; ~1.0) end
  end // end of [gIntSetMap_find]
*)

end // end of [local]


fn stropt_is_GE1 {i:int} (stropt: stropt i):<> bool (i >= 1) 
  = if stropt_is_some(stropt) then let
        val slen = string_length(stropt_unsome(stropt))
      in
        slen >= 1
      end
    else false

datatype GRTOK =
  | TKgene of string
  | TKand of ()
  | TKor of ()
  | TKlpar of ()
  | TKrpar of ()
  | TKEND of ()

datavtype GREXP = 
  | GRgenes of genes
  | GRconj of genes
  | GRconj of (GREXP,GREXP)
  | GRdisj of genes
  | GRdisj of (GREXP,GREXP)

typedef GRtokenizer = '{
  peek = () -<cloref1> GRTOK,
  next = () -<cloref1> void
}

fun GREXP_copy(gr: !GREXP): GREXP = case+ gr of
  | GRgenes (g) => GRgenes (genes_copy(g))
  | GRconj (g) => GRconj (genes_copy(g))
  | GRdisj (g) => GRdisj (genes_copy(g))
  | GRconj (lx, rx) => GRconj(GREXP_copy(lx), GREXP_copy(rx))
  | GRdisj (lx, rx) => GRdisj(GREXP_copy(lx), GREXP_copy(rx))

fun GREXP_free(gr: GREXP): void = case+ gr of 
  | ~GRgenes(g) => genes_free(g)
  | ~GRconj(g) => genes_free(g)
  | ~GRdisj(g) => genes_free(g)
  | ~GRconj(lx,rx) => let
    val () = GREXP_free(lx)
    val () = GREXP_free(rx)
    in () end
  | ~GRdisj(lx,rx) => let 
    val () = GREXP_free(lx)
    val () = GREXP_free(rx)
    in () end
                                   
extern
fun toCaps(c:char): char

extern
fun minConj(gr: GREXP, emap: &gDMap): GREXP

extern
fun dCd(ex1: !GREXP, ex2: !GREXP): string

extern
fun eCd(ex1: !GREXP, ex2: !GREXP): string

extern
fun dCe(ex1: !GREXP, ex2: !GREXP): string

extern
fun grexp_to_string(e0: !GREXP): string 
//extern
//fun stringInOpIze (gset: !genes, inop:string): string

extern
fun print_tks(ss: sstream): void 

extern
fun print_grexp (grexp: GREXP): void 

extern
fun print_pretty_grexp (grexp: !GREXP): void 

extern
fun subsstr_caps(ss: sstream, p1: size_t, p2: size_t): string

extern
fun whileCharTst
  (ss: sstream, chtst: char -<fun1> bool): void

//fn string_of_string1 {n:nat} (str : string n): string = str

(* Recursive Descent parser functions *)
extern
fun parseD(Tzr: GRtokenizer): GREXP 

extern //called first
fun parseE(Tzr: GRtokenizer): GREXP  

extern
fun parseF(Tzr: GRtokenizer): GREXP 

extern
fun parseGR (Tzr: GRtokenizer):  GREXP 
//End of parser function declarations

extern
fun toCNF (bexp: GREXP, emap: &gDMap): GREXP


extern
fun getToken (ss: sstream): GRTOK

extern
fun isAlpha(c:char): bool
extern
fun isWhiteSpace(c:char): bool

implement
whileCharTst
  (ss, chtst) = let
  val c = sstream_get (ss)
in
  if chtst (c) then (sstream_inc (ss); whileCharTst (ss, chtst))
end


implement
isAlpha(c:char): bool = (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z')
  || (c >= 'a' && c <= 'z') || c = '_' || c = '.' || c = '-'|| c = '&' || c = '|' 

implement
isWhiteSpace(c:char): bool = not (c=NUL orelse isAlpha(c) orelse (c = '\(') orelse (c = ')') )

extern
fun GRtokenizer_make(ss: sstream): GRtokenizer


//implement 
fun stringInOpIze (gset: !genes, inop: string): string = let
    val glist:List(string) = genes_listize (gset)
    fun loop(astr:string, alist:List(string)):<cloref1> string = case+ alist of
      | list_nil () => astr
      | list_cons(x,xs as list_nil()) => loop(astr + x, xs)
      | list_cons(x,xs) => loop(x + inop + astr, xs)
  in       
    loop("",glist)    
  end


//call 2x for standard conjuntion X set, make a diff function calling it
//for set X set
fn conjunctivize (cj: GREXP, cjgenes: !genes, emap: &gDMap): GREXP = let
  val glist:List(string) = genes_listize (cjgenes)
  fun loop(cj: GREXP, alist:List(string), emap: &gDMap): GREXP = case+ alist of
    | list_cons(x,xs) => (case+ xs of 
        | list_cons(xx,xxs) => let
          val cjj = GREXP_copy(cj)
          val setx:genes = genes_make_sing(x)
          val newdisj = GRdisj(cjj,GRgenes (setx))
          in 
            GRconj(toCNF(newdisj,emap), loop(cj,xs,emap))
          end
        | list_nil () => toCNF (GRdisj(cj, GRgenes(genes_make_sing(x))),emap)
      ):GREXP 
    | list_nil () => ($raise EmptyList; cj)
    val retGR = loop(cj,glist,emap) 
    //val _ = linset_free(cjgenes)       
  in retGR end

fn conj1(ex1:GREXP, ex2:GREXP, gs:genes, emap: &gDMap): GREXP = let
    val lx1g = conjunctivize (ex1, gs, emap)
    val lx2g = conjunctivize (ex2, gs, emap)
    val _ = genes_free(gs)      
  in GRconj(lx1g,lx2g) end 

//Hope this isn't called too much...:  
fn conj2(gs1: genes, gs2: genes, emap: &gDMap): GREXP = let
  val glist: List(string) = genes_listize (gs1)
  fun loop(alist: !List(string), gs2: !genes, emap: &gDMap)
  : GREXP = case+ alist of
    | list_cons(x,xs) => (case+ xs of
      | list_cons(xx,xxs) =>
        GRconj(conjunctivize(GRgenes (genes_make_sing(x)), gs2, emap), loop(xs,gs2,emap))
      | list_nil () => conjunctivize(GRgenes (genes_make_sing(x)), gs2, emap)
    ):GREXP
    | list_nil () => ($raise EmptyList; GRgenes (genes_make_sing("ERROR")))      
  val retGR = loop(glist,gs2,emap) 
  val _ = genes_free(gs1)  
  val _ = genes_free(gs2)  
  in retGR end


implement
dCd(ex1, ex2): string = 
  "(" + grexp_to_string(ex1) + ") and (" + grexp_to_string(ex2) + ")"

implement
eCd(ex1, ex2): string = 
  grexp_to_string(ex1) + " and (" + grexp_to_string(ex2) + ")"

implement
dCe(ex1, ex2): string = 
  "(" + grexp_to_string(ex1) + ") and " + grexp_to_string(ex2)

macdef p_disj(e0,e1,e2) = x where { 
      prval () = fold@ ,(e1) and () = fold@ ,(e2)
      val x =  dCd(,(e1), ,(e2))
      prval () = fold@ ,(e0) }      

implement
grexp_to_string(e0): string = case+ e0 of
  | GRconj (e1, e2) => (case+ (e1, e2) of
    (GRdisj (x), GRdisj (z)) => p_disj(e0, e1, e2)
    | (GRdisj (x, y), GRdisj (z)) => p_disj(e0, e1, e2)   
    | ( GRdisj (z), GRdisj (x, y)) => p_disj(e0, e1, e2)    
    | ( GRdisj (z, w), GRdisj (x, y)) => p_disj(e0, e1, e2)    
    | (_, GRdisj (x)) => eCd(e1, e2)
    | (_, GRdisj (x, y)) => eCd(e1, e2)
    | (GRdisj (x), _) => dCe(e1, e2)
    | (GRdisj (x, y), _) => dCe(e1, e2)  
    | (_, _) =>> grexp_to_string(e1) + " and " + grexp_to_string(e2)
  )
  | GRconj (e1) =>> stringInOpIze(e1, " and ")
  | GRdisj (e1, e2) => (case+ (e1, e2) of 
    | (GRgenes (s1), GRgenes (s2)) => 
      stringInOpIze(s1, "") + " or " + stringInOpIze(s2, "")
    | (GRgenes (s1), _) => stringInOpIze(s1, "") + " or " + grexp_to_string(e2)
    | (_, GRgenes (s2)) => grexp_to_string(e1) + " or " + stringInOpIze(s2, "")
    | (_, _) => grexp_to_string(e1) + " or " + grexp_to_string(e2) 
  )
  | GRdisj (e1) =>> stringInOpIze(e1, " or ")
  // Could be a single gene required for this reaction
  | GRgenes(s) => let 
    val ssize = genes_size(s)
    val () = if ssize != 1 then $raise InvalidCase;
    val g1: string = ""
    var gene: string 
    val () = assertloc (genes_choose(s, gene))      
    prval any = opt_unsome {string} gene
    in gene end

implement
print_pretty_grexp (grexp): void = () where {
(* Idea is to put parentheses around a sequence
   of disjunctions, to improve the ease of identifying 
   conjunctions and disjunctions, and improve readability in general. *)
   
  val grstr = grexp_to_string(grexp)
  val grstr = grstr + "\n"
  val _ = print(grstr)

}


implement
print_tks(ss: sstream): void = 
  let
    fun loop():<cloref1> void = 
      let
        val tk = getToken(ss) 
      in case+ tk of 
        | TKEND() => () (* print("END\n") *)
        | TKand() => (print ("AND "); loop())
        | TKor() => (print("OR "); loop())
        | TKlpar() => (print("( "); loop())
        | TKrpar() => (print(") "); loop())
        | TKgene(g) => (print(g+" "); loop())
      end
  in
    loop()      
  end

implement
toCaps(c): char = if (c >= 'a' && c <= 'z') then int2char0(char2int0(c)-32) else c
  
implement 
subsstr_caps(ss: sstream, p1: size_t, p2): string = let
    fun loop(i: size_t):<cloref1> string = if i < p2 then let       
        val x = sstream_substr(ss,i,i+1)
        val c = toCaps(string_test_at(g1ofg0(x), 0))
       in
        tostring(c) + loop(i+1)
      end
      else ""
  in
    loop(p1)
  end

//Scanner & Lexer
implement
getToken (ss) = let
//
val () = whileCharTst (ss, isWhiteSpace)
//
val p0 = sstream_pos (ss)
val c0 = sstream_getinc (ss)

//
in
//
case+ 0 of
| _ when c0 = NUL => TKEND ()
| _ when c0 = '\(' => TKlpar ()
| _ when c0 = ')' => TKrpar ()
| _ => let
    val () = whileCharTst (ss, isAlpha)
    val p1 = sstream_pos (ss)
    val word = sstream_substr (ss, p0, p1)
    val wcaps = subsstr_caps(ss, p0, p1)
  in
    case+ 0 of
    | _ when wcaps = "AND" || wcaps = "&" || wcaps = "&&" => TKand ()
    | _ when wcaps = "OR" || wcaps = "|" || wcaps = "||" => TKor ()
    | _ => TKgene (word)
  end // end of [_]
//
end // end of [getToken]



implement 
GRtokenizer_make(ss: sstream): GRtokenizer = let
  val ssr = ref<sstream> ss  
  val curTK = ref<GRTOK> (getToken(!ssr))
  val fpeek = lam() : GRTOK =<cloref1> !curTK
  val fnext = lam() : void =<cloref1> 
    !curTK := getToken(!ssr)
in '{peek = fpeek, next = fnext}
end // end of [GRtokenizer_make]


(* An unambiguous grammar ***
//E: Expression; partial sum (disjunction)
//F: Factor (Conjunction)
//D: Disjunctive subterm (which could be a terminal, i.e. just a gene)
//gid: Gene ID
//As it happens, this grammar appears isomorphic to the one in the 
//SML tutorial.
E -> E or F | F :: E -> F {or F}* //Non-recursive version
F -> F and D | D
D -> gid | (E)    *)


(* *********************** Parser function implementations *************************)


//Does not appear to be tail recursive due to TKrpar check
implement 
parseD(Tzr: GRtokenizer): GREXP = case+ Tzr.peek () of
  | TKgene (g) => (Tzr.next (); GRgenes(genes_make_sing(g)))  (* skip by the token and build datatype result *)
  | TKlpar () => let
      val Dtok = (Tzr.next (); parseE(Tzr))
      val _ = case+ Tzr.peek () of 
        | TKrpar() => Tzr.next ()
        | _ => ($raise ParserParensUnbalanced; ())
    in Dtok end
  | _ => $raise ParserIllegalToken


implement
parseF (Tzr: GRtokenizer): GREXP = let
  fun loop(term: GREXP):<cloref1> GREXP = case+ Tzr.peek () of
    | TKand () => (Tzr.next(); loop( GRconj(term, parseD(Tzr)) )) (* another Disjunction *) 
    | _ => term 
  in loop (parseD (Tzr)) end (* first, parse the first D in the list *)
                                     
implement                                       
parseE (Tzr: GRtokenizer):  GREXP = let
  fun loop(term: GREXP):<cloref1> GREXP = case+ Tzr.peek () of
      | TKor () => (Tzr.next(); loop( GRdisj(term, parseF(Tzr)) ))
      | _ => term
  in loop( parseF(Tzr) ) end

implement
parseGR (Tzr: GRtokenizer): GREXP = let 
  val term = parseE(Tzr)
  val _ = case+ Tzr.peek () of 
    | TKEND () => ()
    | _ =>  ($raise ParserPrematureEND; ())
  in
    term
  end


fun list_min (inlist: List string, emap: &gDMap): string = let 
  fun loop(cmin: string, rlist: List string, emap: &gDMap): string = case+ rlist of
    | list_cons (x, xs) => let
        val xval = gDMap_find(emap, x)
        val cminval = gDMap_find(emap, cmin)
        //Negative values mean the gene wasn't in the dataset:
        val xval = if isnan(xval) != 0 then cminval else xval
      in
        if cminval < xval then loop(cmin,xs,emap) else loop(x,xs,emap) 
      end
    | list_nil () => cmin    
in case+ inlist of
  | list_cons (x, xs) => loop(x,xs,emap) 
  | list_nil () => ($raise EmptyList; " ")  
end // end of [list_min]

(* Computes the (imputed) sum and variance of a gene list (in a disjunctive set).
   For missing values, we assume the mean value.           
   This function and its cousins could benefit from proofs. *)
     
fun dlist_sum_var (inset: !genes, emap: &gDMap, smap: &gDMap): (double, double) = let 
  val num_genes = g0uint2int_size_int(genes_size(inset))
  val inlist = genes_listize(inset)
  fun loop(rlist: List(string), emap: &gDMap, smap: &gDMap, miss: int, csum: double,
  cvar: double): (int, double, double) = case+ rlist of
    | list_cons (x, xs) => let
        val xval = gDMap_find(emap, x)
        val sval = pow(gDMap_find(smap, x), 2.0)
        //Negative values mean the gene wasn't in the dataset:
        val miss = if isnan(xval) != 0 then miss+1 else miss 
        val csum = if isnan(xval) != 0 then csum else csum + xval
        //Assume independent variables for now
        val cvar = if isnan(xval) != 0 then cvar else (cvar + sval)                 
      in
        loop(xs,emap,smap,miss,csum,cvar)
      end
    | list_nil () => (miss, csum, cvar)    
in case+ inlist of
  | list_cons (x, xs) => let 
      val (miss, csum, cvar) = loop(inlist, emap, smap, 0, 0.0, 0.0)
    in  ( g0int2float_int_double(num_genes)*csum/(g0int2float_int_double(num_genes - miss)), cvar) end  
  | list_nil () => ( $raise EmptyList; (0.0, 0.0))  
end // end of [dlist_sum_var]


// Need to split this in to two functions: recursive function to merge 
// disjunctions, then a function to compute values. 
// This function is designed to work on expressions in CNF.


fun disj_vals(ex: !GREXP, emap: &gDMap, smap: &gDMap): (double, double) = let
  fun cg2d(gr:GREXP): GREXP = case+ gr of
    | GRdisj(_) => (fold@ gr; gr)
    | ~GRconj(gs) => let 
        val sz = genes_size(gs)
        in if size_of_int(1) = sz then GRdisj(gs)
        else ($raise InvalidCase; GRconj(gs)) end
    | ~GRgenes(gs) => let 
        val sz = genes_size(gs)
        in if size_of_int(1) = sz then GRdisj(gs)
        else ($raise InvalidCase; GRgenes(gs)) end
    | X => X
  fun loop(ex: GREXP, emap: &gDMap, smap: &gDMap): GREXP = let
    val glist = (case+ ex of
      | ~GRdisj(ex1,ex2) => (case+ (ex1,ex2) of
        | (~GRgenes(s1), ~GRgenes(s2)) => loop(GRdisj (s1+s2),emap,smap)
        | (~GRdisj(s1), ~GRgenes(s2)) => loop(GRdisj (s1+s2),emap,smap)
        | (~GRgenes(s1), ~GRdisj(s2)) => loop(GRdisj (s1+s2),emap,smap)
        | (~GRdisj(s1), ~GRdisj(s2)) => loop(GRdisj (s1+s2),emap,smap)
        | (~GRdisj(ex11, ex12), ex2) => 
          loop(GRdisj(loop(GRdisj(ex11, ex12),emap,smap),ex2),emap,smap)
        | (ex1, ~GRdisj(ex21, ex22)) => 
           loop(GRdisj(ex1, loop(GRdisj(ex21, ex22),emap,smap)),emap,smap) 
        | (LL,RR) => GRdisj(LL,RR)
        )                 
      | ~GRconj(ex1,ex2) => let
        val ex1dv = disj_vals(ex1, emap, smap)
        val ex2dv = disj_vals(ex2, emap, smap) 
        val minex = (if ex1dv.0 < ex2dv.0 then GREXP_copy(ex1)
        else GREXP_copy(ex2)): GREXP
        val _ = GREXP_free(ex1)
        val _ = GREXP_free(ex2)
        in minex end
      | X => X 
      ):GREXP
  in glist end  // end of [loop]
  val exc = cg2d(GREXP_copy(ex))
  val exl = cg2d(loop(exc,emap,smap))
  val dvals = case+ exl of
    | ~GRdisj(gs) => let 
      val dv = dlist_sum_var (gs, emap, smap)
      val _ = genes_free(gs)
    in dv end
    | X => (print("TEST0\t"); print_pretty_grexp(X); 
           GREXP_free(X); $raise TestCase; (0.0,0.0)) 
in 
  dvals
end // end of [disj_vals]

(*        | (~GRconj(s1), ~GRgenes(s2)) => minConj(GRconj (s1+s2),emap)
        | (~GRgenes(s1), ~GRconj(s2)) => minConj(GRconj (s1+s2),emap)
        | (~GRconj(s1), ~GRconj(s2)) => minConj(GRconj (s1+s2),emap)     *)


implement
minConj(gr, emap): GREXP = let
  fun g2c(gr:GREXP): GREXP = case+ gr of
    | GRconj(_) => (fold@ gr; gr)
    | ~GRgenes(gs) => GRconj(gs)
    | X => X
  in (case+ gr of 
      | ~GRconj(ex1,ex2) => let
      val (ex1,ex2) = (g2c(ex1),g2c(ex2))
      in (case+ (ex1,ex2) of
        | (~GRconj(s1), ~GRconj(s2)) => minConj(GRconj (s1+s2),emap)
        | (~GRconj(ex11, ex12), ex2) => 
          minConj(GRconj(minConj(GRconj(ex11, ex12),emap),ex2),emap)
        | (ex1, ~GRconj(ex21, ex22)) => 
          minConj(GRconj(ex1, minConj(GRconj(ex21, ex22),emap)),emap) 
        | (LL,RR) => GRconj(LL,RR)
        ):GREXP end                
      | ~GRconj(gs) => let
            val glist = genes_listize(gs)
            val mingene = list_min(glist, emap)
            val _ = genes_free(gs)
          in GRconj(genes_make_sing(mingene)) end
      | X => X
      ): GREXP
  end // end of [minConj]

implement
toCNF (bexp, emap): GREXP = let     
  val LR:GREXP = (case+ bexp of 
    | ~GRconj(ex1,ex2) => GRconj (toCNF(ex1,emap),toCNF(ex2,emap))
    | ~GRdisj(ex1,ex2) => GRdisj (toCNF(ex1,emap),toCNF(ex2,emap))   
    | GR => GR):GREXP
  in (case+ LR of  
    | ~GRconj(ex1,ex2) => minConj(GRconj(ex1,ex2),emap) 
    | ~GRdisj(ex1,ex2) => (case+ (ex1,ex2) of
      // Handle disjunctive leaf cases:         
      | (~GRdisj(lx), ~GRgenes(g)) => GRdisj (lx + g) 
      | (~GRgenes(g), ~GRdisj(rx)) => GRdisj (rx + g)
      | (~GRdisj(lx), ~GRdisj(rx)) => GRdisj (lx + rx)
      | (~GRgenes(g1), ~GRgenes(g2)) => GRdisj (g1 + g2)

      // Distribute OR over ANDs:
      | (~GRconj(x1,x2), ~GRconj(g)) => conj1(x1,x2,g,emap) 
      | (~GRconj(g), ~GRconj(x1,x2)) => conj1(x1,x2,g,emap) 
      | (~GRconj(g1), ~GRconj(g2)) => conj2(g1,g2,emap)
      | (~GRconj(lx1,lx2), ~GRconj(rx1,rx2)) => let
        val lx1c = GREXP_copy(lx1)
        val lx2c = GREXP_copy(lx2) 
        val rx1c = GREXP_copy(rx1)
        val rx2c = GREXP_copy(rx2) 
        in GRconj(GRconj(GRconj(toCNF(GRdisj(lx1,rx1),emap), 
          toCNF(GRdisj (lx2, rx1c ),emap)),
          toCNF(GRdisj (lx1c, rx2),emap)), 
          toCNF(GRdisj (lx2c, rx2c),emap)) 
        end

      // Handle e.g.: (.. OR ..) OR (.. AND ...) 
      | (~GRconj(lx1,lx2), RX) => let
        val RXc = GREXP_copy(RX)
        in GRconj(toCNF(GRdisj(lx1,RX),emap),
          toCNF(GRdisj(lx2,RXc),emap)) 
        end
      | (LX ,~GRconj(rx1,rx2)) => let
        val LXc = GREXP_copy(LX)
        in GRconj(toCNF(GRdisj(LX,rx1),emap),
          toCNF(GRdisj(LXc,rx2),emap)) 
        end
      | (~GRconj(gc), RX) => let
        val retGR = toCNF(conjunctivize(RX, gc,emap),emap)
        val _ = genes_free(gc)
        in retGR end
      | (LX, ~GRconj(gc)) => let
        val retGR = toCNF (conjunctivize(LX, gc,emap),emap)
        val _ = genes_free(gc)
        in retGR end
      // All other disjunctive cases
      | (_,_) => GRdisj(toCNF(ex1,emap),toCNF(ex2,emap))
      ):GREXP
    | EX => EX
    ):GREXP

  end
  
%{^
#define __sscanf0(line, gene, exp, std) \
  sscanf(line, "%s\t%lf\t%lf", (char*)gene, exp, std)
%}  
  
implement main0 {n} (argc, argv) = () where  {
  // val () = gc_chunk_count_limit_max_set (~1) // infinite
  val () = assertloc(argc = 3)
  val expInFi = argv[1]
  val rulesInFi = argv[2]

(*
  val pargv = &argv

  
  prval (pf, fpf) = __assert(pargv) where {
    //why extern here? try removing it - maybe for multifile support (but we do have multiple files)
    extern praxi __assert {l:addr} (p: ptr l): (ptrarr n @ l, ptrarr n @ l -<lin,prf> void)
  }
  
  prval (pf1, pf2) = ptrarr_uncons(pf)
  *)
  val inFIEXP = fileref_open_exn (expInFi, file_mode_r)
  val inFIRUL = fileref_open_exn(rulesInFi, file_mode_r)

  var ExpMap = gDMap_make_nil ()
  var STDMap = gDMap_make_nil () //Actually variances
    
  val _ = fileref_get_line_string(inFIEXP) // Assume column name line  
  fun loopDATA(emap: &gDMap, smap: &gDMap):<cloref1> void = let
    val linein = fileref_get_line_string(inFIEXP)
    in 
      if strptr_isnot_null(linein) then let
        #define BSZ 30
        var !p_gene = @[byte][BSZ]()
        var exp: double?
        var std: double?
        var resexp: double?
        var resstd: double?
//
        val _ =
        __sscanf0 (
          linein, !p_gene, exp, std
        ) where {
          extern fun __sscanf0{l:agz}
          (
            line: !strptr l
          , gene: &(@[byte?][BSZ]) >> @[byte][BSZ]
          , exp: &double? >> double, std: &double? >> double
          ) : int = "mac#__sscanf0"
        } (* end of [val] *)
//
        val gene = $UN.cast{String}(p_gene)
        val nstr = string1_length (gene)
        val gene = string_make_substring (gene, 0, nstr)
        // Do we need this in ATS2?:
        //val gene = string_of_strbuf (gene)
//        val _ = println! ("gene = ", gene)
//
        val _ = gDMap_insert(emap, gene, exp) 
        val _ = gDMap_insert(smap, gene, std)
        val _ = if exp < 0.0 orelse std < 0.0 then $raise InvalidData;        
//
        val _ = strptr_free(linein)
      in
        loopDATA(emap, smap)
      end
      else (strptr_free(linein))
    end
  val _ = loopDATA(ExpMap,STDMap)
    
  fun loopCNF(emap: &gDMap, smap: &gDMap, cnt: int):<cloref1> void = let
    val linein = fileref_get_line_string(inFIRUL)
    val cnt = cnt+1
  in
    if stropt_is_some(linein) then
      if stropt_is_GE1(linein) then let
        val rexp = parseGR(GRtokenizer_make(sstream_make(stropt_unsome(linein),16384)))        
        val rexpCNF = minConj(toCNF(rexp, emap),emap)
        val dv = disj_vals(rexpCNF, emap, smap);
      in
        ( (* print_pretty_grexp(rexpCNF); print("\t"); *) 
         print(dv.0); print("\t"); print(sqrt(dv.1)); print("\n");
         GREXP_free(rexpCNF); loopCNF(emap, smap, cnt))
      end
      else (print "nan\tnan\n"; loopCNF(emap, smap, cnt))
    else ()
  end    
  val _ = loopCNF(ExpMap, STDMap, 0)

  (*
  prval pf = ptrarr_cons(pf1, pf2)
  prval () = fpf(pf)  *)
  val () = gDMap_free(ExpMap)

  val () = gDMap_free(STDMap)
  val () = exit(0)
}
