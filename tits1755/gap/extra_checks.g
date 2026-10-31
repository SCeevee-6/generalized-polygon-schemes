# Three checks run with GAP 4.16.1, ctbllib 1.3.11 (2026-10-31). Output: extra_checks_output.txt
LoadPackage("ctbllib");;
SetPrintFormattingStatus("*stdout*", false);;

# 1. Complex conjugation and Frobenius-Schur indicators on the irreducibles of 2F4(2)'
t := CharacterTable("2F4(2)'");;
irr := Irr(t);;
Print("degrees ", List(irr, x->x[1]), "\n");;
Print("conj ", List(irr, x->Position(irr, ComplexConjugate(x))), "\n");;
Print("FS ", Indicator(t,2), "\n");;

# 2. Is the Tits group in the Breuer-Lux database shipped with ctbllib (tst/multfree.dat)?
Read(Filename(DirectoriesPackageLibrary("ctbllib","tst"), "multfree.dat"));;
Print("MULTFREEINFO entries ", List(RecNames(MULTFREEINFO), x -> [x, MULTFREEINFO.(x)[1]]), "\n");;
s := StringFile(Filename(DirectoriesPackageLibrary("ctbllib","tst"), "mferctbl.gap"));;
Print("mferctbl.gap mentions 2F4: ", PositionSublist(s, "2F4") <> fail, "\n");;

# 3. Primitive permutation characters of T.2 = 2F4(2)
t2 := CharacterTable("2F4(2)'.2");;
res := List(Maxes(t2), function(nm) local mm, pp; mm := CharacterTable(nm); pp := TrivialCharacter(mm)^t2;
  return [nm, pp[1], ScalarProduct(t2, pp, pp), Maximum(List(Irr(t2), y -> ScalarProduct(t2, pp, y)))]; end);;
Print("T.2 maxes [name, degree, rank, max mult] ", res, "\n");;
res := List(["L2(25).2_3", "2^2.[2^9]:S3", "2.[2^9]:5:4"], function(nm) local mm, pp; mm := CharacterTable(nm);
  pp := TrivialCharacter(mm)^t2;
  return [nm, Filtered(List(Irr(t2), y -> [y[1], ScalarProduct(t2, pp, y), Indicator(t2, [y], 2)[1]]), z -> z[2] > 0)]; end);;
Print("T.2 constituents [degree, mult, FS] ", res, "\n");;
