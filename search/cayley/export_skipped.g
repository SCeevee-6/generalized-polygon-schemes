# Rational character data for the groups of order 750 and 1000 that export_normal.g skipped (normal cyclic Sylow).
# Run in WSL:  gap -q export_skipped.g < /dev/null > tables_skipped.txt
# Rational character data of the nonabelian groups of a given order, for normal_cayley.py.
# For each group: rational class sizes (identity first), the norms <x,x> = |orbit| of the rationalized
# irreducibles (trivial first), and the matrix of rationalized character values on the rational classes.
# Groups with a normal cyclic Sylow subgroup are skipped: every primitive central S-ring over them is trivial
# (Chen-Muzychuk-Ponomarenko, arXiv:1505.01352, Thm 1.3), so they carry no primitive normal distance-regular Cayley graph
# other than K_n.
# Run: for i in [1..NrSmallGroups(n)] do DumpNC(n, i); od;;  (done in batches of 20-40 groups per GAP call)
# Output lines:  A n i  (abelian) | W n i p (normal cyclic Sylow p) | G n i [sizes] [norms] [[values]]
SetPrintFormattingStatus("*stdout*", false);;
DumpAll := function(n, i)
  local G, p, P, tbl, irr, rat, sz, ncl, keys, ratcl, c, key, pos, norms;
  G := SmallGroup(n, i);
  if IsAbelian(G) then Print("A ", n, " ", i, "\n"); return; fi;
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

DumpAll(750,1);;
DumpAll(750,2);;
DumpAll(750,3);;
DumpAll(750,7);;
DumpAll(750,8);;
DumpAll(750,9);;
DumpAll(750,10);;
DumpAll(750,11);;
DumpAll(750,12);;
DumpAll(750,13);;
DumpAll(750,14);;
DumpAll(750,15);;
DumpAll(750,16);;
DumpAll(750,17);;
DumpAll(750,18);;
DumpAll(750,19);;
DumpAll(750,20);;
DumpAll(750,21);;
DumpAll(750,22);;
DumpAll(750,24);;
DumpAll(750,25);;
DumpAll(750,31);;
DumpAll(750,32);;
DumpAll(750,33);;
DumpAll(750,34);;
DumpAll(750,35);;
DumpAll(750,36);;
DumpAll(750,37);;
DumpAll(750,38);;
DumpAll(1000,1);;
DumpAll(1000,3);;
DumpAll(1000,4);;
DumpAll(1000,5);;
DumpAll(1000,6);;
DumpAll(1000,7);;
DumpAll(1000,8);;
DumpAll(1000,10);;
DumpAll(1000,11);;
DumpAll(1000,12);;
DumpAll(1000,13);;
DumpAll(1000,22);;
DumpAll(1000,23);;
QUIT;
