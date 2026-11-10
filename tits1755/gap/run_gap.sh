#!/bin/sh
# GAP checks for T = 2F4(2)'. Needs GAP 4 with ctbllib and atlasrep. ATLAS files are fetched into this folder if missing.
# tits_maxes.g runs all 8 maxes by default; set MAXRANGE := [i];; before reading it to run one.
set -e
cd "$(dirname "$0")"
A=https://brauer.maths.qmul.ac.uk/Atlas/exc/TF42
for f in gap/TF42G1-p1755B0.g1 gap/TF42G1-p1755B0.g2 words/TF42G1-max1W1 words/TF42G1-max2W1 words/TF42G1-max3W1 words/TF42G1-max4W1 words/TF42G1-max5W1 words/TF42G1-max6W1 words/TF42G1-max7W1 words/TF42G1-max8W1; do
  [ -f "$(basename $f)" ] || curl -sSf -o "$(basename $f)" "$A/$f"
done
echo 'Read("tits_chars.g");; QUIT;' | gap -q
echo 'Read("tits_ctG.g");; QUIT;' | gap -q
echo 'Read("tits_maxes.g");; QUIT;' | gap -q
