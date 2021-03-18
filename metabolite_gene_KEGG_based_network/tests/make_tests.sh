
head -n 1 __full.magi_compound_results.csv > magi_compound_results.csv
awk 'BEGIN{OFS=FS=","; COUNT=0} NR>1{COUNT=COUNT+1; $2="gene"COUNT; print}' __full.magi_compound_results.csv | shuf -n 100 >> magi_compound_results.csv

head -n 1 __full.magi_gene_results.csv > magi_gene_results.csv
awk 'BEGIN{OFS=FS=","; COUNT=0} NR>1{COUNT=COUNT+1; $2="gene"COUNT; print}' __full.magi_gene_results.csv | shuf -n 100 >> magi_gene_results.csv


## Compound - Filter by compound_score >= 1, reciprocal_score = 2, e_score_r2g > 5, e_score_g2r > 5
# 1 MAGI_score
# 2 gene_id
# 3 original_compound
# 4 neighbor
# 5 note
# 6 compound_score
# 7 level
# 8 homology_score
# 9 reciprocal_score
# 10 reaction_connection
# 11 e_score_r2g
# 12 database_id_r2g
# 13 e_score_g2r
# 14 database_id_g2r
# 15 Unnamed: 0
# 16 original_mz
# 17 feature
# 18 ppm_error
# 19 searched_adduct
# 20 adj
awk -F',' 'BEGIN{print "feature\toriginal_compound\toriginal_mz"} NR>1 && $6>=1 && $9==2 && $11>5 && $13>5 {print $17"\t"$3"\t"$16}' magi_compound_results.csv > magi_compound_results.filtered.csv



## Gene - Filter by compound_score >= 1, reciprocal_score = 2, e_score_r2g > 5, e_score_g2r > 5
# 1 MAGI_score
# 2 gene_id
# 3 original_compound
# 4 neighbor
# 5 note
# 6 compound_score
# 7 level
# 8 homology_score
# 9 reciprocal_score
# 10 reaction_connection
# 11 e_score_r2g
# 12 database_id_r2g
# 13 e_score_g2r
# 14 database_id_g2r
awk -F',' 'BEGIN{print "gene_id\treaction_id"} NR>1 && $6>=1 && $9==2 && $11>5 && $13>5 {print $2"\t"$14}' magi_gene_results.csv > magi_gene_results.filtered.csv



