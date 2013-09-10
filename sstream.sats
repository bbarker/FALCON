(* ****** ****** *)
//
// Some GC (garbage collector) string stream operations. 
// For convenience, the I/O string type is returned where applicable,
// though the user may easily modify the source to use a different
// string type internally if desired (default is just a string).
//
(* ****** ****** *)



#define NUL '\000'



abstype sstream // for string streams


(* ****** ****** *)
//
// Create a string stream starting at position 0
// of the string.
// 
// Posmax is  checked during a write operation; if posmax is
// <= 0 then no action is taken. Otherwise, 
// if the current position in the string stream is greater
// than posmax during a write operation, then the part 
// of the string before pos is freed and pos is reset to 0.
//
(* ****** ****** *)

fun sstream_make (str: string, posmax: size_t): sstream

(* ****** ****** *)
//
// Get the char at the current position,
// without incrementing the position.
//
(* ****** ****** *)

fun sstream_get (ss: sstream): char


(* ****** ****** *)
// 
// Increment the position by 1.
//
(* ****** ****** *)

fun sstream_inc (ss: sstream): void


(* ****** ****** *)
//
// Get the char at the current position, 
// then increment the position by 1.
//
(* ****** ****** *)

fun sstream_getinc (ss: sstream): char


(* ****** ****** *)
//
// Return current position in the sstream
//
(* ****** ****** *)

fun sstream_pos (ss: sstream): size_t


(* ****** ****** *)
//
// Return all characters between i and j inclusive in 
// the sstream without incrementing pos.
//
(* ****** ****** *)

fun sstream_substr (ss: sstream, i: size_t, j: size_t): [n:nat] string n


(* ****** ****** *)
//
// Read all characters from current pos to end of sstream 
// without incrementing pos.
//
(* ****** ****** *)

fun sstream_read (ss: sstream): string


(* ****** ****** *)
//
// Write (append) all characters from the string to the end 
// of the sstream.
//
(* ****** ****** *)

fun sstream_write {n:nat} (ss: sstream, appstr: string n): void 


  