#!/bin/sh
set -e
cd "$(dirname "$0")"
gap -q tits_chars.g < /dev/null
gap -q tits_maxes.g < /dev/null
ls -l tits_chars_output.txt tits_chartable.csv eigen_characters_gap.csv tits_maxes_output.txt maxes_coherent_configurations.csv
