# Action of 2F4(2)' on the cosets of TableOfMarks class TOMNR (TomLib), printed as images of the two
# generators of UnderlyingGroup(tom). Prepend 'TOMNR := 416;;' etc. Used for classes 416, 417, 418, 422, 423, 424;
# outputs in tits1755/cosets_tom/tom<TOMNR>.txt.

LoadPackage("tomlib");;
tom := TableOfMarks("2F4(2)'");;
G := UnderlyingGroup(tom);;
gens := GeneratorsOfGroup(G);;
H := RepresentativeTom(tom, TOMNR);;
hom := FactorCosetAction(G, H);;
n := Index(G, H);;
SetPrintFormattingStatus("*stdout*", false);;
Print("TOM ", TOMNR, " NGENS ", Length(gens), "\n");;
Print("SIZE ", Size(G), " ", Size(H), " ", n, "\n");;
Print("GA ", ListPerm(Image(hom, gens[1]), n), "\n");;
Print("GB ", ListPerm(Image(hom, gens[2]), n), "\n");;
if Length(gens) > 2 then Print("GC ", ListPerm(Image(hom, gens[3]), n), "\n");; fi;
