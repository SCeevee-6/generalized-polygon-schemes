# Rational character data of the nonabelian groups of a given order, for normal_cayley.py.
# For each group: rational class sizes (identity first), the norms <x,x> = |orbit| of the rationalized
# irreducibles (trivial first), and the matrix of rationalized character values on the rational classes.
# Groups with a normal cyclic Sylow subgroup are skipped: every primitive central S-ring over them is trivial
# (Chen-Muzychuk-Ponomarenko, arXiv:1505.01352, Thm 1.3), so they carry no primitive normal distance-regular Cayley graph
# other than K_n.
# Run: for i in [1..NrSmallGroups(n)] do DumpNC(n, i); od;;  (done in batches of 20-40 groups per GAP call)
# Output lines:  A n i  (abelian) | W n i p (normal cyclic Sylow p) | G n i [sizes] [norms] [[values]]
SetPrintFormattingStatus("*stdout*", false);;
DumpNC := function(n, i)
  local G, p, P, tbl, irr, rat, sz, ncl, keys, ratcl, c, key, pos, norms;
  G := SmallGroup(n, i);
  if IsAbelian(G) then Print("A ", n, " ", i, "\n"); return; fi;
  for p in Set(Factors(n)) do
    P := SylowSubgroup(G, p);
    if IsCyclic(P) and IsNormal(G, P) then Print("W ", n, " ", i, " ", p, "\n"); return; fi;
  od;
  tbl := CharacterTable(G);
  irr := Irr(tbl);
  rat := RationalizedMat(irr);
  sz := SizesConjugacyClasses(tbl);
  ncl := Length(sz);
  keys := []; ratcl := [];
  for c in [1..ncl] do
    key := List(rat, x -> x[c]);
    pos := Position(keys, key);
    if pos = fail then Add(keys, key); Add(ratcl, [c]); else Add(ratcl[pos], c); fi;
  od;
  norms := List(rat, x -> ScalarProduct(tbl, x, x));
  Print("G ", n, " ", i, " ", List(ratcl, o -> Sum(sz{o})), " ", norms, " ",
        List(rat, x -> List(ratcl, o -> x[o[1]])), "\n");
end;;
