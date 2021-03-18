#!/usr/bin/env bash

set -eu

R="rn00290"
./../scripts/KEGG_reaction_KGML_to_network_format.py -x <(wget "http://rest.kegg.jp/get/$R/kgml" -O - -o "test_data/$R.wget.log") --edges __$R.edges.txt --nodes __$R.nodes.txt

diff $R.edges.txt __$R.edges.txt
diff $R.nodes.txt __$R.nodes.txt

