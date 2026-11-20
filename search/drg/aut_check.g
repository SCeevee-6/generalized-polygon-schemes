# Automorphism group of a graph given as adjacency lists (written by known_graphs.py), with GRAPE.
# Usage:  F := "path/to/brouwer3.g";; Read("aut_check.g");;
LoadPackage("grape");; SetPrintFormattingStatus("*stdout*", false);;
Read(F);;
gr := Graph(Group(()), [1..Length(ADJ)], OnPoints, function(x, y) return y in ADJ[x]; end, true);;
A := AutomorphismGroup(gr);;
Print(Length(ADJ), " |Aut| ", Size(A), " transitive ", IsTransitive(A, [1..Length(ADJ)]),
      " primitive ", IsPrimitive(A, [1..Length(ADJ)]), " blocksize ", Length(Blocks(A, [1..Length(ADJ)])[1]),
      " rank ", Length(OrbitsDomain(Stabilizer(A, 1), [1..Length(ADJ)])), "\n");;
