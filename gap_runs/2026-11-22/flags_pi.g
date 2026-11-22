# flags_pi.g : permutation character of T on the 8775 flags, class by class
LoadPackage("ctbllib");;
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;;
G := Group(a, b);;
deg := 1755;;
S := Stabilizer(G, 1);;
nb := First(Orbits(S, [1..deg]), o -> Length(o) = 10);;
NbOf := function(y) local g; g := RepresentativeAction(G, 1, y); return Set(List(nb, x -> x^g)); end;;
pairs := Filtered(Combinations(nb, 2), p -> p[2] in NbOf(p[1]));;
L0 := Set([1, pairs[1][1], pairs[1][2]]);;
Flags := Orbit(G, [1, L0], function(f, g) return [f[1]^g, OnSets(f[2], g)]; end);;
Print("flags ", Length(Flags), "\n");;
fixedflags := function(g)
  return Number(Flags, f -> [f[1]^g, OnSets(f[2], g)] = f);
end;;
C := ConjugacyClasses(G);;
Print("classes ", Length(C), "\n");;
sizes := List(C, Size);;
vals := List(C, c -> fixedflags(Representative(c)));;
Print("class sizes ", sizes, "\n");;
Print("pi values ", vals, "\n");;
Print("<pi,1> ", Sum(List([1..Length(C)], i -> sizes[i]*vals[i])) / Size(G), "\n");;
Print("<pi,pi> ", Sum(List([1..Length(C)], i -> sizes[i]*vals[i]^2)) / Size(G), "\n");;
t := CharacterTable("2F4(2)'");;
degs := List(Irr(t), x -> x[1]);;
Print("library degrees ", SortedList(degs), "\n");;
Print("flag-algebra degrees in library ", List([1,78,351,624,650,675,2048], d -> d in degs), "\n");;
QUIT;
