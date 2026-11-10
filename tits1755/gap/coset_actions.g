# Permutation action of 2F4(2)' on the cosets of maximal subgroup number MAXNR.
# Run from the folder holding TF42G1-p1755B0.g1/.g2 and TF42G1-max*W1 (ATLAS v3), e.g.
#   echo 'MAXNR := 5;;' | cat - coset_actions.g | gap -q > ../cosets/max5.txt
# The lines GA/GB are the images of the standard generators a, b as lists on 1..n.
# Used for MAXNR in 1,2,4,5,6,7,8 (max3 is the 1755-point action, handled elsewhere).
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;; G := Group(a, b);;
prg := ScanStraightLineProgram(Concatenation("TF42G1-max", String(MAXNR), "W1"));;
M := Group(ResultOfStraightLineProgram(prg.program, [a, b]));;
hom := FactorCosetAction(G, M);;
n := Index(G, M);;
SetPrintFormattingStatus("*stdout*", false);;
Print("SIZE ", Size(G), " ", Size(M), " ", n, "\n");;
if n < 5000 then Print("RANK ", Length(DoubleCosetRepsAndSizes(G, M, M)), "\n");; fi;
Print("GA ", ListPerm(Image(hom, a), n), "\n");;
Print("GB ", ListPerm(Image(hom, b), n), "\n");;
