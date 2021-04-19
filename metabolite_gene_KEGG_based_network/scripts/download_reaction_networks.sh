#!/bin/bash

SCRIPT=$(readlink -f "$0")      # Absolute path to this script, e.g. /path/to/data/foo.sh
SCRIPTPATH=$(dirname "$SCRIPT") # Absolute path this script is in, thus /path/to/data

## Get all KEGG Reaction Networks
DIR="kgml"
ALL_NODES="KEGG_Pathway_Networks.nodes.txt"
ALL_EDGES="KEGG_Pathway_Networks.edges.txt"

mkdir -p "$DIR"; cd "$DIR"
rm -f "$ALL_NODES" "$ALL_EDGES"

## Get 'br08901.keg' which lists all KEGG Pathways
wget -O "KEGG_Pathway_Maps_br08901.keg" "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir="
awk '$1~"^C"' "KEGG_Pathway_Maps_br08901.keg" | sed -e 's/^C[ ]*//' -e 's/,//' -e 's/(//g' -e 's/)//g' -e 's@/@@g' -e 's@\\@@g' > "KEGG_Pathway_Maps_br08901_lvlC.txt"

## Get kgml reaction network file for each KEGG Pathway (if one exists)
while read line;
do
        ID=$(echo "rn$line" | awk '{print $1}');
        NAME=$(echo "$line" | sed -e 's/ /_/g');
        echo "Processing $ID   $NAME";
        wget "http://rest.kegg.jp/get/$ID/kgml" -O "$NAME.kgml" -o "$NAME.wget.log"
	if [ -s "$NAME.kgml" ]; then
        	"$SCRIPTPATH/KEGG_reaction_KGML_to_network_format.py" -x "$NAME.kgml" --edges "$NAME.edges.txt" --nodes "$NAME.nodes.txt"
		awk -F'\t' -vNAME="$NAME" '{ if(NR==1) {print $0} else {print NAME$0} }' >> "$ALL_NODES"
		awk -F'\t' -vNAME="$NAME" '{ if(NR==1) {print $0} else {print NAME$1"\t"NAME$2} }' >> "$ALL_EDGES"
	else
		echo "   - No reaction network found!"
        	rm -f "$NAME.kgml" "$NAME.wget.log"
	fi
done < "KEGG_Pathway_Maps_br08901_lvlC.txt"


