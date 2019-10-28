import os, sys
import collections

import pandas as pd

def RNA_share_pathway(type="up"):
    score_dict = {}
    week2_up = "/home/zhihl/Project/CRC/rna_seq_git/Figure/diff_genes_2w_%s.txt" % type
    df_week2_up = pd.read_csv(week2_up, header=0, sep="\t")
    df_week2_up = df_week2_up[['ID','p.adjust']]
    
    
    week4_up = "/home/zhihl/Project/CRC/rna_seq_git/Figure/diff_genes_4w_%s.txt" % type
    df_week4_up = pd.read_csv(week4_up, header=0, sep="\t")
    df_week4_up = df_week4_up[['ID','p.adjust']]
    
    week7_up = "/home/zhihl/Project/CRC/rna_seq_git/Figure/diff_genes_7w_%s.txt" % type
    df_week7_up = pd.read_csv(week7_up, header=0, sep="\t")
    df_week7_up = df_week7_up[['ID','p.adjust']]
    
    week10_up = "/home/zhihl/Project/CRC/rna_seq_git/Figure/diff_genes_10w_%s.txt" % type
    df_week10_up = pd.read_csv(week10_up, header=0, sep="\t")
    df_week10_up = df_week10_up[['ID','p.adjust']]


    df_total =  pd.concat([df_week2_up, df_week4_up, df_week7_up, df_week10_up], ignore_index=True)
    descript =  dict(zip(df_total["ID"], df_total["p.adjust"]))
    words = collections.Counter(list(df_total["ID"]))
    #print(words)
    #print(descript)
    
    key_list = []
    for key in words :
        if words[key] == 4:
            #print ("%s\t%s" % (key , descript[key]))
            key_list.append(key)
            
    return(key_list)
    
    
def enhancer_pathway(type = "up"):
    #enhancer
    enhancer_colits_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/enhancer/colits_%s.txt" % type, header=0, sep="\t")
    enhancer_colits_up = enhancer_colits_up[['ID','p.adjust']]
    
    enhancer_cancer_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/enhancer/cancer_%s.txt" % type, header=0, sep="\t")
    enhancer_cancer_up = enhancer_cancer_up[['ID','p.adjust']]
    
    enhancer_all_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/enhancer/all_%s.txt" % type, header=0, sep="\t")
    enhancer_all_up = enhancer_all_up[['ID','p.adjust']]
    
    df_total =  pd.concat([enhancer_colits_up, enhancer_cancer_up, enhancer_all_up], ignore_index=True)
    descript =  dict(zip(df_total["ID"], df_total["p.adjust"]))
    words = collections.Counter(list(df_total["ID"]))
    
    key_list = []
    for key in words :
        if words[key] == 3:
            #print ("%s\t%s" % (key , descript[key]))
            key_list.append(key)
            
    return(key_list)
   
    
    
    
def promoter_pathway(type = "up"):
    #promoter
    promoter_colits_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/promoter/colits_%s.txt" % type, header=0, sep="\t")
    promoter_colits_up = promoter_colits_up[['ID','p.adjust']]
    
    promoter_cancer_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/promoter/cancer_%s.txt" % type, header=0, sep="\t")
    promoter_cancer_up = promoter_cancer_up[['ID','p.adjust']]
    
    promoter_all_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/promoter/all_%s.txt" % type, header=0, sep="\t")
    promoter_all_up = promoter_all_up[['ID','p.adjust']]
    
    df_total =  pd.concat([promoter_colits_up, promoter_cancer_up, promoter_all_up], ignore_index=True)
    descript =  dict(zip(df_total["ID"], df_total["p.adjust"]))
    words = collections.Counter(list(df_total["ID"]))
    
    key_list = []
    for key in words :
        if words[key] == 3:
            #print ("%s\t%s" % (key , descript[key]))
            key_list.append(key)
            
    return(key_list)
    
def select_all_GO(go_list=None):
    week2_up = "/home/zhihl/Project/CRC/rna_seq_git/Figure/diff_genes_2w_%s.txt" % type
    df_week2_up = pd.read_csv(week2_up, header=0, sep="\t")
    df_week2_up = df_week2_up[['ID','p.adjust']]
    df_week2_up[df_week2_up["ID"]  ]

    



if __name__ == "__main__":
    RNA_list = RNA_share_pathway(type="up")
    #print(len(RNA_list))

    enhancer_list = enhancer_pathway(type="up")
    #print(len(key_list))

    promoter_list = promoter_pathway(type="up")
    #print(len(promoter_list))
    
    share = list(set(RNA_list) & set(enhancer_list) & set(promoter_list))
    print (share)
    
    promoter_all_up = pd.read_csv("/home/zhihl/Project/CRC/rna_seq_git/Figure/promoter/all_%s.txt" % "up", header=0, sep="\t")
    for one in share:
        print(promoter_all_up[["ID", "Description"]][one == promoter_all_up["ID"]])

