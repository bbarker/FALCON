
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

staload "sstream.sats"


//
// Probably best to include mystring as a separate .dats file?
// Is there a better way?
//

(* ***** GC-Requiring Strings ***** *)
(*         Uncomment to use.       *)

typedef mystring = string
typedef mystring (n:int) = string n


extern
castfn mystr1_of_mystr
(str: mystring):<> [n:nat] mystring n

symintr mystr_length
fn mystr0_length (str: mystring):<> size_t = string0_length(str)
fn mystr1_length {n:nat} (str: mystring n):<> [n1:nat | n == n1] size_t n1 = string1_length(str)
overload mystr_length with mystr0_length
overload mystr_length with mystr1_length


fn mystr_substr {ln,st,m:nat | st+ln <= m} 
  (str: mystring m, st: size_t st, ln:size_t ln):<> [n:nat | n == ln ] mystring n = let 
  val m = mystr1_length(str) 
in
  string1_of_strbuf(string_make_substring(str, st, ln))
end

symintr static_size
extern
prfun static_mystring_size {n:nat} (s: mystring n):<> [n2:nat | n2 == n] void
overload static_size with static_mystring_size

extern
prfun static_size1_size {n:nat} (s: !size_t n):<> [n2:nat | n2 == n] void
overload static_size with static_size1_size

extern
castfn sizeLte_of_size1 {n,m: nat | n <= m} 
  (n: !size_t n, m: !size_t m):<> sizeLte(m) 

fn mystring_append
  (s1: mystring, s2: mystring):<> mystring 
  = string0_append(s1,s2) 
overload + with string0_append

fn mystring1_append {i,j:nat}
  (s1: mystring i, s2: mystring j):<> mystring (i+j) // linear
  = string1_of_strbuf(string1_append(s1, s2))
overload + with mystring1_append 
// overload conflicts with string n overload

(* ***** End of GC-Requiring Strings ***** *)



local // local of gc-friendly sstream

typedef
sstream_struct (n:int, m:int) = @{
data = mystring (n)
, size = size_t (n)
, pos = sizeLte (n)
, posmax = size_t (m)
}            

typedef sstream_struct = [n,m:nat] sstream_struct (n, m)

fun strstream_make (str: mystring, posmax: size_t): sstream_struct = let
  val [n:int] str = mystr1_of_mystr (str)
  val n = mystr_length (str)
  val pos = size1_of_int1(0)
in
  @{data=str, size=n, pos=pos, posmax = size1_of_size(posmax)}
end

assume sstream = ref (sstream_struct)

in (* in of [local] *)
  


implement
sstream_make (str: mystring, posmax: size_t): sstream = 
  ref<sstream_struct> (strstream_make(str, posmax))


implement
sstream_pos (ss) = let
  val (vbox pf | p) = ref_get_view_ptr (ss) in p->pos
end // end of [sstream_pos]
  

implement
sstream_get (ss) = let   
  val (vbox pf | p) = ref_get_view_ptr (ss) 
  val strdat = p->data
  in if p->pos < p->size then strdat[p->pos] else NUL 
end // end of [sstream_get]     


implement
sstream_inc (ss) = let
  val (vbox pf | p) = ref_get_view_ptr (ss)
  in if p->pos < p->size then p->pos := p->pos + 1 
end

implement  
sstream_getinc (ss) = let
  val (vbox pf | p) = ref_get_view_ptr (ss)
  val strdat = p->data
  //Is this optimal with the redundant conditional?
  val c = if p->pos < p->size then strdat[p->pos] else NUL
  val _ = if p->pos < p->size then p->pos := p->pos + 1
in c
end


implement
sstream_read (ss) = let
  val (vbox pf | p) = ref_get_view_ptr (ss)
in p->data
end // end of [sstream_read]  

fun sstream_testsz (ss: sstream): void = let
  val (vbox pf | p) = ref_get_view_ptr (ss)
  prval [n:int] _ = static_size(p->size)
in  () end


implement
sstream_write (ss, appstr) = let
  val (vbox pf | p) = ref_get_view_ptr (ss)
  val [d:int] apsz = mystr1_length(appstr)
  prval [n:int] _ = static_size(p->size)
  prval [p:int] _ = static_size(p->pos)
  prval [pmx:int] _ = static_size(p->posmax)
  prval [ss:int] _ = static_size(p->data)

  // Perform freeing before appending 
  // to avoid high memory use; likely only useful
  // in strbuf version though.
  val _ = (if p->posmax > 0 then let 
    val x = (if p->pos + apsz > p->posmax then let
      val [d2:int] strtmp = mystr_substr(p->data, p->pos , (p->size - p->pos))
      val _ = p->data := strtmp 
      val _ = p->size := (p->size - p->pos) 
      val _ = p->pos := (size1_of_int1(0))   
    in () end) // [end of if p->pos + apsz > p->posmax] 
  in () end) // [end of if p->posmax > 0]  
//  This castfn should be abstracted above:
  val _ = p->data := mystring1_append(p->data,appstr)
  val _ = p->size := p->size + apsz 
//  val _ = p->pos := postmp
in
 () 
end // end of [sstream_write]


implement
sstream_substr (ss: sstream, i: size_t, j: size_t): [n:nat] mystring n =
  let
    val (vbox pf | p) = ref_get_view_ptr (ss)
    val pd = j-i
    val i = size1_of_size(i)
    val pd = size1_of_size(pd)
  in
    if i + pd <= p->size then
      mystr_substr(p->data, i, pd)
      //string1_of_strbuf(string_make_substring(p->data,i,pd))
    else ""
  end       

end (* end of [local] *)  