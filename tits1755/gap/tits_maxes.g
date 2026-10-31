# Shared helpers. Run from the directory holding TF42G1-p1755B0.g1/.g2 and TF42G1-max*W1.
LoadPackage("ctbllib");;
LoadPackage("atlasrep");;
OpenOut := function(name) local s; s := OutputTextFile(name, false); SetPrintFormattingStatus(s, false); return s; end;;
CsvCell := function(x) if IsString(x) and (x <> [] or IsStringRep(x)) then return Concatenation("\"", x, "\""); fi; return Concatenation("\"", String(x), "\""); end;;
WriteCsv := function(name, header, rows)
  local s, r;
  s := OpenOut(name);
  WriteAll(s, JoinStringsWithSeparator(List(header, CsvCell), ",")); WriteAll(s, "\n");
  for r in rows do WriteAll(s, JoinStringsWithSeparator(List(r, CsvCell), ",")); WriteAll(s, "\n"); od;
  CloseStream(s);
end;;
t := CharacterTable("2F4(2)'");;
if t = fail then Error("the GAP character table library has no table named 2F4(2)'"); fi;
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;;
G := Group(a, b);;
CLASSWORDS := [
  ["1A", a^2],
  ["2A", a],
  ["2B", Comm(a,b*a*b)^2],
  ["3A", b],
  ["4A", (a*b*a*b*a*b^2)^3],
  ["4B", a*b*a*b*a*b*a*b*a*b^2*a*b*a*b^2*a*b^2],
  ["4C", a*b*a*b*a*b^2*a*b^2],
  ["5A", a*b*a*b^2],
  ["6A", a*b*a*b*a*b*a*b*a*b^2],
  ["8A", (a*b)^5*a*b^2*(a*b*a*b*a*b^2)^2*a*b^2],
  ["8B", (a*b)^5*a*b^2*a*b^2*(a*b*a*b*a*b^2)^2],
  ["8C", a*b*a*b*a*b*a*b^2],
  ["8D", a*b*a*b*(a*b*a*b*a*b*a*b^2)^2],
  ["10A", a*b*a*b*a*b*a*b^2*a*b*a*b^2*a*b^2],
  ["12A", a*b*a*b*a*b^2],
  ["12B", a*b*a*b*a*b*a*b^2*a*b^2],
  ["13A", a*b],
  ["13B", (a*b)^2],
  ["16A", (a*b)^5*a*b^2*a*b*(a*b^2)^3],
  ["16B", (a*b)^5*(a*b^2)^3*a*b*a*b^2],
  ["16C", (a*b)^3*(a*b^2)^3*a*b*a*b*a*b^2],
  ["16D", a*b*a*b*a*b*a*b^2*a*b*a*b*(a*b^2)^3]
];;

if not IsBound(MAXRANGE) then MAXRANGE := [1..8]; fi;;
LOG := OpenOut(Concatenation("tits_maxes_output_", JoinStringsWithSeparator(List(MAXRANGE, String), "_"), ".txt"));;
W := function(arg) local x; for x in arg do if IsString(x) then WriteAll(LOG, x); else WriteAll(LOG, String(x)); fi; od; WriteAll(LOG, "\n"); end;;
W("GAP ", GAPInfo.Version, ", ctbllib ", InstalledPackageVersion("ctbllib"), ", atlasrep ", InstalledPackageVersion("atlasrep"));
cn := ClassNames(t);; ords := OrdersClassRepresentatives(t);; cent := SizesCentralizers(t);;
irr := List(Irr(t), ValuesOfClassFunction);; degs := List(irr, x -> x[1]);; chn := CharacterNames(t);;
ind := Indicator(t, 2);;
pos := List(CLASSWORDS, p -> Position(cn, LowercaseString(p[1])));;
if fail in pos then pos := List(CLASSWORDS, p -> Position(cn, p[1])); fi;
reps := [];; for k in [1..Length(CLASSWORDS)] do reps[pos[k]] := CLASSWORDS[k][2]; od;
ClassOf := function(g)
  local cands, o, c;
  o := Order(g); c := Size(Centralizer(G, g));
  cands := Filtered([1..Length(cn)], k -> ords[k] = o and cent[k] = c);
  if Length(cands) = 1 then return cands[1]; fi;
  return First(cands, k -> IsConjugate(G, g, reps[k]));
end;;
NAMES := ["L3(3):2", "L3(3):2", "2.[2^8].5.4", "L2(25)", "2^2.[2^8].S3", "A6.2^2", "A6.2^2", "5^2:4A4"];;
libpi := List(Maxes(t), nm -> [nm, ValuesOfClassFunction(Induced(CharacterTable(nm), t, [TrivialCharacter(CharacterTable(nm))])[1])]);;
W("Library maxes: ", List(libpi, x -> x[1]));
rows := [];;
M := 0;; cnt := 0;; pi := 0;; mlt := 0;; dc := 0;; k := 0;;
for i in MAXRANGE do
  tt := Runtime();
  prg := ScanStraightLineProgram(Concatenation("TF42G1-max", String(i), "W1"));
  M := Group(ResultOfStraightLineProgram(prg.program, [a, b]));
  d := Size(G) / Size(M);
  # permutation character of G on cosets of M from the class distribution of M
  cnt := List(cn, x -> 0);
  for c in ConjugacyClasses(M) do k := ClassOf(Representative(c)); cnt[k] := cnt[k] + Size(c); od;
  pi := List([1..Length(cn)], k -> cent[k] * cnt[k] / Size(M));
  mlt := List(irr, x -> ScalarProduct(t, ClassFunction(t, x), ClassFunction(t, pi)));
  rk := Sum(mlt, x -> x^2);
  nsp := Sum([1..Length(mlt)], k -> ind[k] * mlt[k]);
  lib := List(Filtered(libpi, x -> x[2] = pi), x -> x[1]);
  # direct count: (M,M)-double cosets = orbitals of the coset action; self-paired test g^-1 in MgM
  dc := DoubleCosetRepsAndSizes(G, M, M);
  sp := Number(dc, x -> x[1]^-1 in DoubleCoset(M, x[1], M));
  cons := List(Filtered([1..Length(mlt)], k -> mlt[k] <> 0), k -> Concatenation(chn[k], "^", String(mlt[k])));
  Add(rows, [i, NAMES[i], Size(M), d, Length(dc), rk, Number(mlt, x -> x <> 0), Maximum(mlt),
             ForAll(mlt, x -> x <= 1), sp, nsp, sp = Length(dc), JoinStringsWithSeparator(cons, " + "), JoinStringsWithSeparator(lib, ";")]);
  W("max", i, " ", NAMES[i], ": |M| = ", Size(M), ", degree ", d, ", rank (double cosets) ", Length(dc), ", <pi,pi> ", rk,
    ", max multiplicity ", Maximum(mlt), ", self-paired orbitals ", sp, " (character count ", nsp, "), library table ", lib,
    ", ", Int((Runtime() - tt) / 1000), " s");
od;
WriteCsv(Concatenation("maxes_rows_", JoinStringsWithSeparator(List(MAXRANGE, String), "_"), ".csv"),
  ["max", "structure", "order", "degree", "rank_double_cosets", "rank_inner_product", "distinct_constituents", "max_multiplicity",
   "commutative", "self_paired_orbitals", "self_paired_by_indicators", "symmetric", "constituents", "library_max_table"], rows);
CloseStream(LOG);
