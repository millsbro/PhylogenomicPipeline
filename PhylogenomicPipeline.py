### PHYLOGENOMIC PIPELINE PYTHON SCRIPT ###

## 0. PRECODE ##

# Import Modules #

import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

# Create string objects of input and output directory paths #

input_dir="/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
output_dir="/scratch/millsbro/Week6/PhylogenomicPipeline/brass_align/"


# Make a list of files in Input directory #

catalogue=glob.glob(input_dir+"*fasta")

# Create a shorter list for testing purposes #

#catalogue=catalogue[0:100]

#print(catalogue)


## 1. PERFORM MSA WITH MAFFT ##

for item in catalogue:
    aligned_item_path=item.replace(input_dir,output_dir)
    alignment_command='mafft --auto --quiet '+item+' > '+aligned_item_path
    #print(alignment_command)
    os.system(alignment_command)


## 2. PHYLOGENOMIC INFERENCE WITH IQTREE ##

alignment=glob.glob(output_dir+"*fasta")

for aln_temp in alignment:
    tree_command = f"iqtree -s {aln_temp} -m TEST -nt 2"
    os.system(tree_command)
#print(tree_command)


## 3. TREE TOPOLOGY ANALYSIS ##

trees_list=glob.glob(output_dir+"*treefile")
#print(trees_list)

counts=[]

for list_item in trees_list:
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(list_item, "newick")

    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
	    if "Es_" in tip.name:
		    es_tip = tip
		    break
    
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
    
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
             Bs_temp=t 
        elif "Cr_" in t.name:
              Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
        
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
    

    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    counts.append(topo_str)    

#print(counts)

BsCr_counts=(counts.count("12top"))
BsAt_counts=(counts.count("13top"))
CrAt_counts=(counts.count("23top"))

print(f'The number of trees with Bs-Cr sisters is {BsCr_counts}')
print(f'The number of trees with Bs-At sisters is {BsAt_counts}')
print(f'The number of trees with Cr-At sisters is {CrAt_counts}')
