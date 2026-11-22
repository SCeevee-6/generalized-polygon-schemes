# flags_check.g : independent check of T = 2F4(2)' on the 8775 flags of GO(2,4)
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;;
G := Group(a, b);;
deg := 1755;;
Print("order ", Size(G), "\n");;
S := Stabilizer(G, 1);;
nb := First(Orbits(S, [1..deg]), o -> Length(o) = 10);;
NbOf := function(y)
  local g;
  g := RepresentativeAction(G, 1, y);
  return Set(List(nb, x -> x^g));
end;;
pairs := Filtered(Combinations(nb, 2), p -> p[2] in NbOf(p[1]));;
Print("lines through point 1: ", Length(pairs), "\n");;
L0 := Set([1, pairs[1][1], pairs[1][2]]);;
Lines := Orbit(G, L0, OnSets);;
Print("lines: ", Length(Lines), "\n");;
act := function(f, g) return [f[1]^g, OnSets(f[2], g)]; end;;
Flags := Orbit(G, [1, L0], act);;
Print("flags: ", Length(Flags), "\n");;
FS := ShallowCopy(Flags);; Sort(FS);;
perm := function(x)
  return PermList(List(FS, f -> PositionSorted(FS, [f[1]^x, OnSets(f[2], x)])));
end;;
Gf := Group(perm(a), perm(b));;
Print("|Gf| ", Size(Gf), "\n");;
S1 := Stabilizer(Gf, 1);; sub := Orbits(S1, [1..Length(FS)]);;
Print("rank ", Length(sub), "\n");;
Print("suborbit lengths ", SortedList(List(sub, Length)), "\n");;
QUIT;
