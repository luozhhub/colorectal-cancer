"""
1. first you need to down GWAS file from: https://www.ebi.ac.uk/gwas/api/search/downloads/full
2. run this script: python ./CRCSNP.py
3. ln -s /home/zhluo/test/GWAS_site_3columns.txt /home/zhluo/Project/CRC/LD_data/ensembl_io/
4. cd /home/zhluo/Project/CRC/LD_data/ensembl_io/
5. cat ./GWAS_site_3columns.txt |sort |uniq > new.list.GWAS.site.txt
6. nohup perl test1.pl &
7. awk -F"\t" '{print $2"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' ld0.8_collection.tab >/home/zhluo/Project/CRC/LD_data/ensembl_io/snpsnap_gene.txt
"""

#!/usr/bin/python
import pandas as pd
import numpy as np

#mix type
df = pd.read_csv("/home/zhluo/test/full", header=0, sep="\t", low_memory=False)

#p value  < 1e-8
new_df = pd.DataFrame()
for index, row in df.iterrows():
    row["p_VALUE"] = float(row["P-VALUE"])
    if(row["DISEASE/TRAIT"].startswith("Colorectal") and row["p_VALUE"]  <= 1E-8):
        new_df = new_df.append(row)
        
#select columns
select_df = new_df.loc[: ,["DISEASE/TRAIT", "REGION", "CHR_ID", "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE", "SNP_GENE_IDS", "SNPS", "P-VALUE", "OR or BETA"]]

ensembl = []
for item in select_df["SNP_GENE_IDS"]:
    array = str(item).split(",")
    ensembl = ensembl + array
    for arr in array:
        print(arr)
select_df.to_csv("CRC_GWAS_site.txt", header=True, sep="\t")


df = pd.read_csv("CRC_GWAS_site.txt", header=0, sep="\t", index_col=0)
new_df = df.loc[:,["CHR_ID",  "CHR_POS", "SNPS"]]
#new_df["end"] = new_df["CHR_POS"]
new_df.to_csv("~/test/GWAS_site_3columns.txt", header=False, sep="\t", index=False)
