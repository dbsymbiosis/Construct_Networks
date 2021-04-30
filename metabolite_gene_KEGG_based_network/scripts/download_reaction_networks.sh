#!/bin/bash

SCRIPT=$(readlink -f "$0")      # Absolute path to this script, e.g. /path/to/data/foo.sh
SCRIPTPATH=$(dirname "$SCRIPT") # Absolute path this script is in, thus /path/to/data

## Get all KEGG Reaction Networks
DIR="kgml"
rm -fr "$DIR"; mkdir -p "$DIR"; cd "$DIR"

ALL_NODES="KEGG_Pathway_Networks.nodes.txt"
ALL_EDGES="KEGG_Pathway_Networks.edges.txt"

rm -f "$ALL_NODES" "$ALL_EDGES"
echo -e "node_id\tkegg_id\tname\ttype\tinfo\tlink\tx\ty\twidth\theight\tshape" > "$ALL_NODES"
echo -e "node_1\tnode_2" > "$ALL_EDGES"

## Get 'br08901.keg' which lists all KEGG Pathways
wget -O "KEGG_Pathway_Maps_br08901.keg" "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir="
awk '$1~"^C"' "KEGG_Pathway_Maps_br08901.keg" | sed -e 's/^C[ ]*//' -e 's/,//' -e 's/(//g' -e 's/)//g' -e 's@/@@g' -e 's@\\@@g' > "KEGG_Pathway_Maps_br08901_lvlC.txt"

## Get kgml reaction network file for each KEGG Pathway (if one exists)
download_network() {
	DIR="${1}"; shift
	ALL_NODES="${1}"; shift
	ALL_EDGES="${1}"; shift
	SCRIPTPATH="${1}"; shift
	line="${@}"
        ID=$(echo "$line" | awk '{print $1}')
        NAME=$(echo "$line" | sed -e 's/ /_/g')
	
	# Print all info to log file
	echo "Started processing $ID   $NAME"
	exec 1> "${NAME}.processing.log" 2>&1
        echo "## Processing $ID   $NAME"
	
        wget "http://rest.kegg.jp/get/rn${ID}/kgml" -O "$NAME.kgml" -o "$NAME.wget.log"
	if [ -s "$NAME.kgml" ]; then
        	"$SCRIPTPATH/KEGG_reaction_KGML_to_network_format.py" -x "$NAME.kgml" --edges "$NAME.edges.txt" --nodes "$NAME.nodes.txt"
		awk -F'\t' -vNAME="$NAME" 'NR>1{print NAME"__"$0}' "$NAME.nodes.txt" >> "$ALL_NODES"
		awk -F'\t' -vNAME="$NAME" 'NR>1{print NAME"__"$1"\t"NAME"__"$2}' "$NAME.edges.txt" >> "$ALL_EDGES"
	else
		echo "##    - No reaction network found!"
		echo "##        - Searching for gene network."
        	rm -f "$NAME.kgml" "$NAME.wget.log"
		
		wget "http://rest.kegg.jp/get/hsa${ID}/kgml" -O "$NAME.kgml" -o "$NAME.wget.log"
		if [ -s "$NAME.kgml" ]; then
			"$SCRIPTPATH/KEGG_reaction_KGML_to_network_format.py" -x "$NAME.kgml" --edges "$NAME.edges.txt" --nodes "$NAME.nodes.txt"
			awk -F'\t' -vNAME="$NAME" 'NR>1{print NAME"__"$0}' "$NAME.nodes.txt" >> "$ALL_NODES"
			awk -F'\t' -vNAME="$NAME" 'NR>1{print NAME"__"$1"\t"NAME"__"$2}' "$NAME.edges.txt" >> "$ALL_EDGES"
		else
			echo "##    - No gene network found!"
			echo "##        - Giving up :("
			rm -f "$NAME.kgml" "$NAME.wget.log"
		fi
	fi
	echo "## Done $ID   $NAME"
}

export -f download_network
parallel -j 12 download_network "$DIR" "$ALL_NODES" "$ALL_EDGES" "$SCRIPTPATH" :::: "KEGG_Pathway_Maps_br08901_lvlC.txt"

echo ""; echo "Done processing networks!"
