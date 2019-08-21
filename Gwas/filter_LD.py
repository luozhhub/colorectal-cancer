"""
this script is used to process Gwas snp
usage:

"""


import os
import pandas as pd
import re
import json

def filter_LD():
    snp_dir = "/home/zhluo/Project/CRC/LD_data/ensembl_io"
    LD_Dict = {}
    population = ["CEU", "JPT"]
    for pop in population:
        #for CEU
        pop_path = os.path.join(snp_dir, pop)
        snp_files = os.listdir(pop_path)    
        for fileO in snp_files:
            print(fileO)
            LD_Dict[fileO] = []
            df = pd.read_csv(os.path.join(pop_path, fileO), sep="\t", header=None)
            for index, row in df.iterrows():
                association = row[1]
                if re.search(fileO, association):
                    #print (association)
                    LD_Dict[fileO] = LD_Dict[fileO] + association.split("-")
            LD_Dict[fileO] = list(set(LD_Dict[fileO]))
            if fileO in LD_Dict[fileO]:
                LD_Dict[fileO].remove(fileO)        
    json.dump(LD_Dict, open('ld_dict.json', 'w'))
        
def read_LD():
    json.load(open('ld_dict.json', 'r'))        



if __name__ == "__main__":
    filter_LD()