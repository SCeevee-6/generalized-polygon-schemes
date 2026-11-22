# GAP runs, 22 Nov 2026

Run in WSL with GAP 4.15.1, ctbllib 1.3.11 and atlasrep 2.1.9. Each script was run from the folder that
holds the ATLAS generator files TF42G1-p1755B0.g1, TF42G1-p1755B0.g2 and TF42G1-max1W1 to TF42G1-max8W1
(not included here; they come from the ATLAS v3 page for 2F4(2)').

Command used for each script:

    gap -q SCRIPT.g < /dev/null

Scripts and outputs:
- tits_chars.g: character table of 2F4(2)' from the library, and the eigenspace characters of the 1755-point action.
- tits_maxes.g: coherent configurations on the cosets of the eight maximal subgroups.
- flags_check.g: flag action of T on the 8775 flags (2925 lines): group order, rank and suborbit lengths.
- flags_pi.g: permutation character of T on the flags, values on the 22 classes, and <pi,1>, <pi,pi>.
- flags_decomp.g: multiplicities of the library irreducibles in the flag permutation character.

The flags_*_output.txt files were written by PowerShell Tee-Object as UTF-16 and are converted to UTF-8 here.

These results come from GAP 4.15.1. The earlier character-table results were computed with GAP 4.16.1 through
the MCP connector, and they have not yet been compared with these.
