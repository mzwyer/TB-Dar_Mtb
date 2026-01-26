import pandas as pd
import numpy as np

matrix = pd.read_csv("~/Phylogenetic_SNPs/data/outputs/merged_tbdar_frozen102020_1082_filtered_filGATK_variablepos_mac12_matrix",
                     sep = "\t")
info = pd.read_csv("~/data/interim/metadata_tbdar_status1_frozen102020_1082.txt",
                     sep = "\t") # metadata from samples published in the following publication: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010893#sec022
matrix = matrix.merge(info[["G_NUMBER", "Sublineage", "Lineage", 'LastTANAnc', "Sub_short", "Introduction"]], left_index = True, right_on = "G_NUMBER")

# 1 corresponds to 0/0 which is the reference, 0 corresponds to missing and 3 to 1/1 which is the alternative
# thus replace it in that way
matrix = matrix.replace(0, np.nan)
matrix = matrix.replace(1, 0)
matrix = matrix.replace(3, 1)

snp_dic_sublin = {}
snp_dic_lin = {}
imports = []
var_pos_dic = {}
snp_dic = {}
for lineage in matrix.Sub_short.unique(): # all positions mutated in one sublineage
    subset = matrix[matrix.Sub_short == lineage]
    pos = list(subset.loc[:, (subset==1).all()].columns) # get all columns where all members of this sublineage have a mutation
    snp_dic_sublin[lineage] = pos
    
for lineage in matrix.Lineage.unique(): # all positions mutated in one lineage
    subset = matrix[matrix.Lineage == lineage]
    pos = list(subset.loc[:, (subset==1).all()].columns) # get all columns where all members of this sublineage have a mutation
    snp_dic_sublin[lineage] = pos
    
for intro in matrix.Introduction.unique(): # the most successful introductions
    pos_var = []
    subset = matrix[(matrix.Introduction == intro)]
    pos = list(subset.loc[:, (subset==1).all()].columns)
    for posi in matrix.columns:
        if len(matrix[(~pd.isnull(matrix[posi]))& (matrix.Introduction != intro)][posi].unique()) > 1:
            if posi != "G_NUMBER":
                pos_var.append(posi)
    identifier = intro
    imports.append(identifier)
    snp_dic_lin[identifier] = pos
    var_pos_dic[identifier] = pos_var
# check the SNPs that are unique to each of the sublineages
unique_snps = {}
for lineage in matrix.Sub_short.unique():
    subdic = [value for key, value in snp_dic_sublin.items() if key not in [lineage]]
    pos_rest = [item for sublist in subdic for item in sublist]
    pos_unique = np.setdiff1d(snp_dic_sublin[lineage], pos_rest)
    unique_snps[lineage] = pos_unique        
        
# now check, which of the positions are uniquely mutated in that import
snp_dic_lin.update(snp_dic_sublin) # also add the sublineage-specific SNPs
for lineage in imports: 
    subdic = [value for key, value in snp_dic_lin.items() if key not in [lineage]]
    pos_rest = [item for sublist in subdic for item in sublist]
    pos_var = var_pos_dic[lineage]
    pos_unique = np.setdiff1d(snp_dic_lin[lineage], pos_rest)
    print(len(pos_unique))
    # also remove the positions that are variable in another introduction or the remaining strains of a lineage
    pos_unique = np.setdiff1d(pos_unique, pos_var)
    print(len(pos_unique))
    unique_snps[lineage] = pos_unique

pos_dic = {}
i = 0
for key in list(unique_snps.keys()):
    for item in list(unique_snps[key]):
        pos_dic[i] = {"Position_ref" : item.split(":")[1].split("_")[0],
                     "ancestral" : item.split(":")[1].split("_")[1].split("/")[0],
                     "derived" : item.split(":")[1].split("_")[1].split("/")[1],
                     "PhylogeneticSNP": key}
        i += 1

# now check in what gene the positions are located and what the genebased position is
annot = pd.read_csv("AnnotationH37Rv.ptt",
                   sep = " ")
annot[["Start", "Stop"]] = annot[["Start", "Stop"]].apply(pd.to_numeric)
# add whether the SNP is synonymous or nonsynonymous and which codon item it is
with open("MTB_ancestor_reference.fasta") as f:
    mtb_fasta = f.read()
mtb_fasta = mtb_fasta[9:].split()
mtb_fasta = ''.join(mtb_fasta)
codons = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

for i in range(len(pos_dic)):
    pos = int(pos_dic[i]["Position_ref"])
    sub = annot[(annot.Stop >= pos) & (annot.Start <= pos)]
    if len(sub) == 1:
        gene = sub.ID.item()
        strand = sub.Strand.item()
        if strand == "+":
            gene_start = sub.Start.item()
            gene_end = sub.Stop.item()
            pos_gene = int(pos) - int(gene_start)+1
            rest = pos_gene%3
            if rest == 1:
                position_codon = "first"
                codon = mtb_fasta[pos-1:pos+2]
                if codon[0] != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[0] = pos_dic[i]["derived"]
                codon_new = "".join(l)
            elif rest == 0:
                position_codon = "third"
                codon = mtb_fasta[pos-3:pos]
                if codon[2] != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[2] = pos_dic[i]["derived"]
                codon_new = "".join(l)
            elif rest == 2:
                position_codon = "second"
                codon = mtb_fasta[pos-2:pos+1]
                if codon[1] != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[1] = pos_dic[i]["derived"]
                codon_new = "".join(l)
            if codons[codon_new] != codons[codon]:
                mutation = "nonsynonymous"
            elif codons[codon_new] == codons[codon]:
                mutation = "synonymous"
        elif strand == "-": 
            gene_end = sub.Start.item()
            gene_start = sub.Stop.item()
            pos_gene = int(gene_start) - int(pos) + 1
            rest = pos_gene%3
            if rest == 1:
                position_codon = "first"
                codon = mtb_fasta[pos-3:pos][::-1].translate(str.maketrans('ATGC', "TACG"))
                if codon[0].translate(str.maketrans('ATGC', "TACG")) != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[0] = pos_dic[i]["derived"].translate(str.maketrans('ATGC', "TACG"))
                codon_new = "".join(l)
            elif rest == 0:
                position_codon = "third"
                codon = mtb_fasta[pos-1:pos+2][::-1].translate(str.maketrans('ATGC', "TACG"))
                if codon[2].translate(str.maketrans('ATGC', "TACG")) != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[2] = pos_dic[i]["derived"].translate(str.maketrans('ATGC', "TACG"))
                codon_new = "".join(l)
            elif rest == 2:
                position_codon = "second"
                codon = mtb_fasta[pos-2:pos+1][::-1].translate(str.maketrans('ATGC', "TACG"))
                if codon[1].translate(str.maketrans('ATGC', "TACG")) != pos_dic[i]["ancestral"]:
                    print("Something is wrong, make sure you have the correct position!")
                l = list(codon)
                l[1] = pos_dic[i]["derived"].translate(str.maketrans('ATGC', "TACG"))
                codon_new = "".join(l)
            if codons[codon_new] != codons[codon]:
                mutation = "nonsynonymous"
            elif codons[codon_new] == codons[codon]:
                mutation = "synonymous"
        else:
            gene_start = sub.Start.item()
            gene_end = sub.Stop.item()
            pos_gene = int(pos) - int(gene_start)
            mutation = np.nan
            codon_new = np.nan
            codon = np.nan
            position_codon = np.nan
        pos_dic[i]["Start"] = gene_start
        pos_dic[i]["End"] = gene_end
        pos_dic[i]["Strand"] = strand
        pos_dic[i]["Gene"] = gene
        pos_dic[i]["Position_gene"] = pos_gene
        pos_dic[i]["Position_codon"] = position_codon
        pos_dic[i]["Kind_mutation"] = mutation
    else: 
        pos_dic[i]["Start"] = np.nan
        pos_dic[i]["End"] = np.nan
        pos_dic[i]["Strand"] = np.nan
        pos_dic[i]["Gene"] = np.nan
        pos_dic[i]["Position_gene"] = np.nan
        pos_dic[i]["Position_codon"] = np.nan
        pos_dic[i]["Kind_mutation"] = np.nan

phylogenetic_snps = pd.DataFrame(pos_dic).T

phylogenetic_snps.to_csv("~/Phylogenetic_SNPs/data/outputs/phylogenetic_SNPs_tbdar_frozen102020_1082.txt",
                     sep = "\t", index = False)