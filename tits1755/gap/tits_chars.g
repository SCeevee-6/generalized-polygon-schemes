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

LOG := OpenOut("tits_chars_output.txt");;
W := function(arg) local x; for x in arg do if IsString(x) then WriteAll(LOG, x); else WriteAll(LOG, String(x)); fi; od; WriteAll(LOG, "\n"); end;;
W("GAP ", GAPInfo.Version, ", ctbllib ", InstalledPackageVersion("ctbllib"), ", atlasrep ", InstalledPackageVersion("atlasrep"));

# ---------- A. library character table ----------
W("Library table: ", Identifier(t), ", |T| = ", Size(t), ", classes = ", NrConjugacyClasses(t));
cn := ClassNames(t);; ords := OrdersClassRepresentatives(t);; cent := SizesCentralizers(t);;
W("Class names: ", cn);
W("Element orders: ", ords);
W("Centralizer orders: ", cent);
irr := List(Irr(t), ValuesOfClassFunction);;
degs := List(irr, x -> x[1]);;
chn := CharacterNames(t);;
W("Character names: ", chn);
W("Degrees (library order): ", degs);
recalled := [1,26,26,27,27,78,300,325,351,351,351,624,624,650,675,702,702,1300,1300,1728,2048,2048];;
W("Sorted degrees equal the recalled list: ", SortedList(degs) = recalled);
W("Sum of squared degrees = ", Sum(degs, d -> d^2), ", equals |T|: ", Sum(degs, d -> d^2) = Size(t));
ind := Indicator(t, 2);;
W("Frobenius-Schur indicators: ", ind);
W("Characters with all values rational: ", Filtered([1..Length(irr)], i -> ForAll(irr[i], IsRat)));

# ---------- group from ATLAS generators, classes from ATLAS words ----------
W("|G| from generators = ", Size(G), ", equals library |T|: ", Size(G) = Size(t));
W("a has order ", Order(a), ", b has order ", Order(b), ", ab has order ", Order(a*b));
pos := List(CLASSWORDS, p -> Position(cn, LowercaseString(p[1])));;
if fail in pos then pos := List(CLASSWORDS, p -> Position(cn, p[1])); fi;
W("ATLAS class names found in library ClassNames: ", not fail in pos);
reps := [];;
for k in [1..Length(CLASSWORDS)] do reps[pos[k]] := CLASSWORDS[k][2]; od;
okord := ForAll([1..Length(cn)], k -> Order(reps[k]) = ords[k]);;
okcent := ForAll([1..Length(cn)], k -> Size(Centralizer(G, reps[k])) = cent[k]);;
W("Word representatives: element orders agree with library: ", okord);
W("Word representatives: centralizer orders in G agree with library: ", okcent);
noncj := ForAll(Combinations([1..Length(cn)], 2), p -> ords[p[1]] <> ords[p[2]] or cent[p[1]] <> cent[p[2]] or not IsConjugate(G, reps[p[1]], reps[p[2]]));;
W("Word representatives pairwise non-conjugate in G: ", noncj);

# ---------- B. permutation character on 1755 points ----------
if HasMaxes(t) then mx := Maxes(t); else mx := []; fi;;
W("Library maxes: ", mx);
pi_lib := fail;; sname := fail;;
for nm in mx do
  s := CharacterTable(nm);
  if Size(t) / Size(s) = 1755 then sname := nm; pi_lib := ValuesOfClassFunction(Induced(s, t, [TrivialCharacter(s)])[1]); fi;
od;
W("Library max of index 1755: ", sname, " (order ", Size(t)/1755, ")");
pi_grp := List(reps, g -> Number([1..1755], x -> x^g = x));;
W("pi from fixed points of class reps: ", pi_grp);
W("pi from library induction 1_H^T:      ", pi_lib);
W("equal: ", pi_grp = pi_lib);
m := List(irr, x -> ScalarProduct(t, ClassFunction(t, x), ClassFunction(t, pi_grp)));;
W("Multiplicities <pi, chi> (library order): ", m);
W("Constituents: ", List(Filtered([1..Length(m)], i -> m[i] <> 0), i -> [chn[i], degs[i], m[i]]));
W("<pi, pi> = ", ScalarProduct(t, ClassFunction(t, pi_grp), ClassFunction(t, pi_grp)), ", multiplicity-free: ", ForAll(m, x -> x <= 1));
stab := Stabilizer(G, 1);;
W("Point stabilizer order ", Size(stab), ", suborbit lengths ", SortedList(List(Orbits(stab, [1..1755]), Length)));

# ---------- A (cont.). eigenspace characters, exact on every class ----------
n := 1755;;
P := [[1,10,80,640,1024],[1,5,10,0,-16],[1,1,-10,-8,16],[1,-3,2,16,-16],[1,-5,20,-80,64]];;
Q := n * P^-1;;
orbs := ShallowCopy(Orbits(stab, [1..n]));;
Sort(orbs, function(x, y) return Length(x) < Length(y); end);;
idx := [];; for j in [1..5] do for x in orbs[j] do idx[x] := j - 1; od; od;
trans := [];; trans[1] := One(G);; queue := [1];; i := 1;;
while i <= Length(queue) do
  for g in [a, b] do
    y := queue[i]^g;
    if not IsBound(trans[y]) then trans[y] := trans[queue[i]] * g; Add(queue, y); fi;
  od;
  i := i + 1;
od;
invtrans := List(trans, Inverse);;
Rel := function(x, y) return idx[y^invtrans[x]]; end;;
# check that relation 1 is the valency-10 graph and relations are distances: sum of A_1 rows
Ncount := function(g) local N, x; N := [0,0,0,0,0]; for x in [1..n] do N[Rel(x, x^g) + 1] := N[Rel(x, x^g) + 1] + 1; od; return N; end;;
eig := List(reps, g -> Ncount(g) * Q / n);;     # eig[k][e]
eigchars := TransposedMat(eig);;
W("Eigenspace degrees chi_e(1A): ", List(eigchars, x -> x[1]));
W("All values integers: ", ForAll(Flat(eigchars), IsInt));
for e in [1..5] do
  p := Position(irr, eigchars[e]);
  W("chi_", e - 1, " (degree ", eigchars[e][1], ") = library ", chn[p], " (position ", p, "), <chi,chi> = ",
    ScalarProduct(t, ClassFunction(t, eigchars[e]), ClassFunction(t, eigchars[e])));
od;
W("Sum of the five equals pi: ", Sum(eigchars) = pi_grp);

# ---------- CSV output ----------
rows := [];;
Add(rows, Concatenation(["#order", ""], ords));
Add(rows, Concatenation(["#centralizer", ""], cent));
Add(rows, Concatenation(["#classsize", ""], List(cent, c -> Size(t) / c)));
if pi_lib <> fail then Add(rows, Concatenation(["#pi_library", ""], pi_lib)); fi;
for k in [1..Length(irr)] do Add(rows, Concatenation([chn[k], degs[k]], irr[k])); od;
WriteCsv("tits_chartable.csv", Concatenation(["character", "degree"], List(cn, UppercaseString)), rows);
rows := [];;
for k in [1..Length(cn)] do Add(rows, Concatenation([UppercaseString(cn[k]), ords[k], cent[k], pi_grp[k]], List(eigchars, x -> x[k]))); od;
WriteCsv("eigen_characters_gap.csv", ["class", "order", "centralizer", "pi", "chi_1", "chi_351", "chi_650", "chi_675", "chi_78"], rows);

CloseStream(LOG);
