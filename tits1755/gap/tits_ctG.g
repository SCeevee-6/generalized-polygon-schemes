# Shared helpers. Run from the directory holding TF42G1-p1755B0.g1/.g2 and TF42G1-max*W1.
LoadPackage("ctbllib");;
LoadPackage("atlasrep");;
OpenOut := function(name) local s; s := OutputTextFile(name, false); SetPrintFormattingStatus(s, false); return s; end;;
CsvCell := function(x) if IsString(x) and (x <> [] or IsStringRep(x)) then return Concatenation("\"", x, "\""); fi; return Concatenation("\"", String(x), "\""); end;;
WriteCsv := function(name, header, rows)
  local s, r;
  s := OpenOut(name);
  WriteAll(s, JoinStringsWithSeparator(List(header, CsvCell), ",")); WriteAll(s, "\n");
  for r in rows do WriteAll(s, JoinStringsWithSeparator(List(r, CsvCell), ",")); WriteAll(s, "\n"); od;
  CloseStream(s);
end;;
t := CharacterTable("2F4(2)'");;
if t = fail then Error("the GAP character table library has no table named 2F4(2)'"); fi;
Read("TF42G1-p1755B0.g1");; Read("TF42G1-p1755B0.g2");;
a := b11;; b := b21;;
G := Group(a, b);;
LOG := OpenOut("tits_ctG_output.txt");;
W := function(arg) local x; for x in arg do if IsString(x) then WriteAll(LOG, x); else WriteAll(LOG, String(x)); fi; od; WriteAll(LOG, "\n"); end;;
# ---------- optional: table of G computed from scratch, compared with the library ----------
tt := Runtime();;
ctG := CharacterTable(G);; Irr(ctG);;
tr := TransformingPermutationsCharacterTables(ctG, t);;
W("Character table of G computed from the generators (Dixon-Schneider) in ", Int((Runtime() - tt) / 1000), " s; equivalent to library table: ", tr <> fail);
CloseStream(LOG);
