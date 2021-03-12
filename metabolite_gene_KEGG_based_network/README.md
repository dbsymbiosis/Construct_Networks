# Metabolite-gene network construction using KEGG pathways

MAGI info: https://magi.nersc.gov/help/

## 1. Setup env
Load conda environment using the below command
```
conda env create -f environment.yml
conda activate networks
```

## 2. Filter MAGI output

Get gene ids and annotated compound InChI Key (compound_score >= 1, e_score_r2g > 5, e_score_g2r > 5, reciprocal_score = 2)
```
awk -F'\t' 'NR>1 && $7>=1 && $10==2 && $12>=5 && $14>=5 {print $3"\t"$4"\t"$13}' "test_data/MAGI_output.txt" > "test_data/test.gene_id-original_compound-rhea_id"
```
Returns a file with `gene_id`, `original_compound` InChI Key and `rhea_id` columns.

## 3. RHEA to KEGG reaction IDs

We need to convert the RHEA_IDs from MAGI into KEGG Reaction IDs so we can link them to our network.
A bunch of IDs will fail becuase they are from other databases (such a meta-cyc), will have to try and handle them at a later date.
```
./get_KEGG_reaction_ID_from_RHEA_ID.py -r <(awk -F'\t' '{print $3}' test_data/test.gene_id-original_compound-rhea_id | sort | uniq) -o test_data/test.rhea_id-reaction_id
```
Returns a file with `rhea_id` (only if ID in rhea-db.org) and KEGG `reaction_id` columns

## 4. KEGG Reaction network to Cytoscape style format

Download a KEGG Reaction network in `KGML` format and expract edge and node info. 
Script will fetch `Name`, `Enzyme`, or `Exact mass` for a given reaction/compound/map link

`Nodes` either a compound (metabolite), reaction (enzyme/gene), or link to another KEGG map.
 - `node_id` id of target node
 - `kegg_id` id of node in KEDD database (can be compound [cpd], reactiopn [rn/rc] or map [path] ids)
 - `name` compound name, reaction EC number, or map name
 - `type` compound/reaction/map
 - `info` compound 'Exact mass'; 'NA' if reaction or map
 - `link` link to kegg.jp site for the given compound/reaction/map
 - `x` node x coords for visualization
 - `y` node y coords for visualization
 - `width` node width for visualization
 - `height` node height for visualization
 - `shape` nade shape for visualization
`Edges` undirected connection between two `nodes`
	`node_1` reaction or map link node id
	`node_2` compount node id

```
R="rn00290"
./KEGG_reaction_KGML_to_network_format.py -x <(wget "http://rest.kegg.jp/get/$R/kgml" -O - -o "test_data/$R.wget.log") --edges test_data/$R.edges.txt --nodes test_data/$R.nodes.txt
```



```
./match_compound_to_known_molWeight.py --known <(awk -F'\t' '$4=="compound" {print $5"\t"$3}' test_data/rn00290.nodes.txt) --unknown test_data/unknown.txt --known_col 1 --unknown_col 2 -o test_data/test.matches.txt
```
