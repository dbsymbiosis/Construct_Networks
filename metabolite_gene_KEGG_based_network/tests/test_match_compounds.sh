#!/usr/bin/env bash

set -eu

../scripts/match_compound_to_known_molWeight.py --known <(awk -F'\t' '$4=="compound" {print $5"\t"$3}' rn00290.nodes.txt) --unknown unknown_metabolites.txt --known_col 1 --unknown_col 2 -o __matched_metabolites.txt

diff matched_metabolites.txt __matched_metabolites.txt

