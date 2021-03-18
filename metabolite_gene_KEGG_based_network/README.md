# Metabolite-gene network construction using KEGG pathways

MAGI info: https://magi.nersc.gov/help/

## 0. Download and pre-process ID mapping files.
This step may have been completed previously and the mapping files might already be avaliable to you. 
In that case move on to **Step 1**. If the files are not avaliable then below are instructions for how to generate these files.

#### 0.1 InChIKey to KEGG Compound mapping file
Download the `unique_compounds.csv.gz` which is provided by the MAGI creators. 
```
wget https://magi.nersc.gov/files/processed/unique_compounds.csv.gz
```
Open the uncompressed `unique_compounds.csv` file in excel and "save as" a tab delimited .txt file. 
This will convert the comma delimited file into a tab delimited file which is much easier to process.

Next run the below command to extract just the "inchi_key" and "kegg_id" columns from the file. 
It also ignores InChIKeys without compound IDs and reformats the compound IDs by removing the `///` and ` // ` characters that seperatre multiple compound IDs, and by removing extra `"` characters.
```
zcat unique_compounds.txt.gz | awk -F'\t' '$9!=""{print $8"\t"$9}' | sed -e 's@"@@g' -e 's@///@,@g' -e 's@ // @,@g' | gzip -c > InChIKey_2_KEGG_Compound_mapping.txt.gz
```
This file now has two columns: `InChIKey` [tab] `KEGG_Compound_ID`
Multiple `KEGG_Compound_ID`'s associated with the same key are seperated by commas.

#### 0.2 RHEA to KEGG Reaction mapping file

```
wget https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz
python ../scripts/process_RHEA_rdf_file.py -r rhea.rdf.gz -o RHEA_2_KEGG_Reaction_mapping.txt.gz
```
This file now has four columns: `RHEA_ID` [tab] `KEGG_Reaction_IDs` [tab] `Meatcyc_IDs` [tab] `EC_Numbers`
Multiple `KEGG_Reaction_IDs`, `Meatcyc_IDs`, `EC_Numbers` associated with the same `RHEA_ID` are seperated by commas.

#### 0.3 MetaCyc to KEGG Reaction mapping file




## 1. Run MAGI1
First step is to run MAGI1 localy using your metabolite and gene information. 

NOTE: 
 - Need to use MAGI1 **NOT** MAGI2 becuase MAGI2 uses SMILES not m/z (MAGI2 would be easier as it uses Retrorules, where as MAGI1 uses rhea & metacyc and is much harder to connect to KEGG pathways)
 - When running MAGI you need to include a `feature` column in the input metabolite data file otherwise you cant easily link the compound annotations back to your original annotations.

## 2. Setup env
Load conda environment using the below command
```
conda env create -f environment.yml
conda activate networks
```

## 3. Download working data
Download mapping files to link InChiKeys with KEGG compounds, and to link rhea & metacyc reaction IDs with KEGG EC numbers.

```
wget https://magi.nersc.gov/files/processed/unique_compounds.csv.gz

```

## 4. KEGG Reaction network to Cytoscape style format

Download a KEGG Reaction network in `KGML` format and expract edge and node info.
Script will fetch `Name`, `Enzyme`, or `Exact mass` for a given reaction/compound/map link

```
R="rn00290"
./scripts/KEGG_reaction_KGML_to_network_format.py -x <(wget "http://rest.kegg.jp/get/$R/kgml" -O - -o "test_data/$R.wget.log") --edges test_data/$R.edges.txt --nodes test_data/$R.nodes.txt
```
Output files:
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

## 5. Filter `magi_gene_results.csv`



## 6. Filter `magi_compound_results.csv`


