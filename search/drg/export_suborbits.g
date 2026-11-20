# Suborbit lengths of all primitive groups of degree V other than Sym(V) and Alt(V).
# Prints "G V i |G| [sorted suborbit lengths]" for each group ("S V i" for Sym and Alt).
# Usage:  V := 729;; Read("export_suborbits.g");;
SetPrintFormattingStatus("*stdout*", false);;
for i in [1..NrPrimitiveGroups(V)] do
  G := PrimitiveGroup(V, i);;
  if not (IsNaturalSymmetricGroup(G) or IsNaturalAlternatingGroup(G)) then
    ol := ShallowCopy(OrbitLengths(Stabilizer(G, 1), [1..V]));; Sort(ol);;
    Print("G ", V, " ", i, " ", Size(G), " ", ol, "\n");;
  else
    Print("S ", V, " ", i, "\n");;
  fi;
od;
