# Primitive groups of the degrees that occur among open SRG parameter sets (Brouwer's table, v <= 1300).
# Prints "v i rank |G|" for every group that is not Sym(v) or Alt(v) and has rank >= 4, and the generators
# ("GRP v i" then one "GEN [images]" line per generator). Set DEGREES before reading this file, e.g.
#   DEGREES := [88, 99, 100];;  Read("export_primitive.g");;
# The list of 430 degrees is the column v of open_srg.csv.
for v in DEGREES do
  for i in [1..NrPrimitiveGroups(v)] do
    G := PrimitiveGroup(v, i);;
    if not (IsNaturalSymmetricGroup(G) or IsNaturalAlternatingGroup(G)) then
      rk := Length(OrbitsDomain(Stabilizer(G, 1), [1..v]));;
      if rk >= 4 then
        Print("RANK ", v, " ", i, " ", rk, " ", Size(G), "\n");;
        Print("GRP ", v, " ", i, "\n");;
        for g in GeneratorsOfGroup(G) do Print("GEN ", ListPerm(g, v), "\n");; od;
      fi;
    fi;
  od;
od;
