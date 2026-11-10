# Coordinates for a primitive affine group G = PrimitiveGroup(V, I) of degree p^n.
# Prints "HDR |G| |N| elementary-abelian regular n rank" and, for every element x of the socle N,
# "PT point suborbit exponents", where point = 1^x, suborbit is the orbit number of that point under the
# stabiliser of 1, and exponents are the coordinates of x with respect to a Pcgs of N.
# Usage:  V := 625;; I := 109;; Read("export_affine.g");;
G := PrimitiveGroup(V, I);; N := Socle(G);; pc := Pcgs(N);;
S := Stabilizer(G, 1);; orbs := OrbitsDomain(S, [1..V]);; lab := [];;
for j in [1..Length(orbs)] do for x in orbs[j] do lab[x] := j;; od;; od;;
SetPrintFormattingStatus("*stdout*", false);;
Print("HDR ", Size(G), " ", Size(N), " ", IsElementaryAbelian(N), " ", IsRegular(N, [1..V]), " ", Length(pc), " ", Length(orbs), "\n");;
for x in N do Print("PT ", 1^x, " ", lab[1^x], " ", ExponentsOfPcElement(pc, x), "\n");; od;;
