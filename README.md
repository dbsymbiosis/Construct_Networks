# Construct_Networks
Scripts for constructing gene and metabolite networks.

`differentially_accumulated_metabolites`
R script to identify metabolites that are differentially accumulated at a specific timepoint between two conditions.

`dynamic_network_edge_select`
Jupyter notebook that allows the user to select which edge types to remove from the network before replotting the network as an interactive HTML file.
This is useful for removing edges that represent common metabolites (e.g. ATP) that join a large number of nodes into a big group (that is not very informative given how common the metabolite is in most reactions in the cell). 

`gene_co-expression_network`
Construct gene co-expression network

`metabolite_co-expression_network`
Construct metabolite co-expression network

`metabolite_gene_KEGG_based_network`
Scripts for downloading and connecting the RHEA numbers to KEGG Reation IDs and to build a KEGG reaction network that we can annotate with info about the gene or metabolite/compound expression/accumulation profiles. 
This should hopefully be useful for connecting timeseries gene and metabolite data. KEGG gives structure to the data and allows us to see which gene expression patterns might be feeding into a given pool or metabolites. 


