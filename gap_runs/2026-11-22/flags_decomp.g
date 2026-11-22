# flags_decomp.g : decomposition of the permutation character of T on the 8775 flags
LoadPackage("ctbllib");;
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;;
G := Group(a, b);;
deg := 1755;;
CLASSWORDS := [
  ["1A", a^2], ["2A", a], ["2B", Comm(a,b*a*b)^2], ["3A", b],
  ["4A", (a*b*a*b*a*b^2)^3], ["4B", a*b*a*b*a*b*a*b*a*b^2*a*b*a*b^2*a*b^2],
  ["4C", a*b*a*b*a*b^2*a*b^2], ["5A", a*b*a*b^2], ["6A", a*b*a*b*a*b*a*b*a*b^2],
  ["8A", (a*b)^5*a*b^2*(a*b*a*b*a*b^2)^2*a*b^2], ["8B", (a*b)^5*a*b^2*a*b^2*(a*b*a*b*a*b^2)^2],
  ["8C", a*b*a*b*a*b*a*b^2], ["8D", a*b*a*b*(a*b*a*b*a*b*a*b^2)^2],
  ["10A", a*b*a*b*a*b*a*b^2*a*b*a*b^2*a*b^2], ["12A", a*b*a*b*a*b^2],
  ["12B", a*b*a*b*a*b*a*b^2*a*b^2], ["13A", a*b], ["13B", (a*b)^2],
  ["16A", (a*b)^5*a*b^2*a*b*(a*b^2)^3], ["16B", (a*b)^5*(a*b^2)^3*a*b*a*b^2],
  ["16C", (a*b)^3*(a*b^2)^3*a*b*a*b*a*b^2], ["16D", a*b*a*b*a*b*a*b^2*a*b*a*b*(a*b^2)^3]
];;
S := Stabilizer(G, 1);;
nb := First(Orbits(S, [1..deg]), o -> Length(o) = 10);;
NbOf := function(y) local g; g := RepresentativeAction(G, 1, y); return Set(List(nb, x -> x^g)); end;;
pairs := Filtered(Combinations(nb, 2), p -> p[2] in NbOf(p[1]));;
L0 := Set([1, pairs[1][1], pairs[1][2]]);;
Flags := Orbit(G, [1, L0], function(f, g) return [f[1]^g, OnSets(f[2], g)]; end);;
fixedflags := function(g)
  return Number(Flags, f -> [f[1]^g, OnSets(f[2], g)] = f);
end;;
t := CharacterTable("2F4(2)'");;
cn := ClassNames(t);; ords := OrdersClassRepresentatives(t);; cent := SizesCentralizers(t);;
chn := CharacterNames(t);;
pos := List(CLASSWORDS, p -> Position(cn, LowercaseString(p[1])));;
if fail in pos then pos := List(CLASSWORDS, p -> Position(cn, p[1])); fi;
reps := [];; for k in [1..Length(CLASSWORDS)] do reps[pos[k]] := CLASSWORDS[k][2]; od;;
Print("reps match library orders and centralisers: ", ForAll([1..22], k -> Order(reps[k]) = ords[k] and Size(Centralizer(G, reps[k])) = cent[k]), "\n");;
Print("reps pairwise non-conjugate: ", ForAll([1..22], j -> ForAll([1..22], k -> j = k or not IsConjugate(G, reps[j], reps[k]))), "\n");;
pival := List([1..22], k -> fixedflags(reps[k]));;
Print("pi on library classes ", pival, "\n");;
irr := Irr(t);;
csz := List(cent, c -> Size(G)/c);;
mult := List(irr, chi -> Sum(List([1..22], k -> csz[k]*pival[k]*ComplexConjugate(ValuesOfClassFunction(chi)[k]))) / Size(G));;
degs := List(irr, x -> x[1]);;
Print("multiplicities ", mult, "\n");;
Print("degree check (should be 8775): ", Sum(List([1..Length(irr)], i -> mult[i]*degs[i])), "\n");;
Print("sum of squares (should be 18): ", Sum(List(mult, m -> m^2)), "\n");;
Print("constituents ", List(Filtered([1..Length(irr)], i -> mult[i] <> 0), i -> [chn[i], degs[i], mult[i]]), "\n");;
QUIT;
