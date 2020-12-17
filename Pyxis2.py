
# coding: utf-8

from GFF_parser import GFF_parser
from HPArrayOfGeneInfo import HPArrayOfGeneInfo
from HPCombinations import HPCombinations
from HPFindClusters import HPFindClusters
from Combination import Combinations
from Patterns import patterns
from statsmodels.sandbox.stats.multicomp import multipletests

import csv
import argparse
import numpy as np
import itertools
import matplotlib.pyplot as plt
import re
import os
import os.path

# In[28]:
parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')
parser.add_argument('-v', metavar='v', type=str, help="GFF file") #GFF file
parser.add_argument('-s', metavar='s', type=str, help="List of regulated genes") #file with upregulated genes
parser.add_argument('-w', metavar='w', type=str, nargs='?', default='None', const='None', help='Statistical test used to detect clusters')
parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')
parser.add_argument('-k', metavar='k', type=int, nargs='?', help='Search string used to identify Gene name in GFF file (default: 1)')
parser.add_argument('-p', metavar='p', type=float, nargs='?', default=0.01, const=0.01, help='Lower P-value used to make statistical test more rigorous (default: 0.01)')
parser.add_argument('-q', metavar='q', type=float, nargs='?',default=0.05, const=0.05, help='Maximum p-value to make statistical test more rigorous (default:0.05)')
parser.add_argument('-l', metavar='l', type=int, nargs='?', default=5, const=5, help='Step size to detect gene clusters (default: 5)')
args = parser.parse_args()

if args.k == 1:
    type_gene = str("Name=")
    new_gene = str("Protein coding genes")
    
if args.k == 0:
    type_gene = str("")
    new_gene = str("Protein coding, dubious and merged ORFs")

print("*******************************PYXIS2***********************************************************")
print("Developed by Louis Conradie")
print("Centre for Bioinformatics and Computational Biology")
print("Head: Prof Hugh Patterton")
print("Stellenbosch University")
print("2020")
print("************************************************************************************************")
print("The search string for the GFF file is: ", args.i)
print("Pyxis2 is now detecting: ", new_gene)
print("The cut-off p-value is: ", args.q)
print("The lowest p-value for most significant clusters is:", args.p)
print("The step size for the detection of clusters is: ", args.l)
print("************************************************************************************************")

GFF_file = args.v
Text_file = args.s

try:
    if os.stat(GFF_file).st_size > 0:
        print("Code proceeding...")
    else:
        pass
except FileNotFoundError:
    print("GFF3 file not found, check file path")
    exit()  

try:
    if os.stat(Text_file).st_size > 0:
        print("Code proceeding...")
    else:
        pass
except FileNotFoundError:
    print("Text file not found, check file path")
    exit()  

try:
    if os.stat(GFF_file).st_size > 0:
        print("GFF3 file found:", GFF_file)
    else:
        pass
except TypeError:
    print("GFF3 file found not found")
    exit()    
    
try:
    if os.stat(Text_file).st_size > 0:
        print("Text file found:", Text_file)
    else:
        pass
except TypeError:
    print("Text file found not found")
    exit()

print("***********************************************************************************************")

GFF_file = args.v
files = args.s
idstring = args.i
low_sig_value = args.p
high_sig_value = args.q
minimum_window = args.l

#files = file that contains the names of upregulated genes

gene = []
#creates accumulator list called gene
files = open(files,"r")
for line in files:
    #for line in file
    lines = line.split()
    #split each line in the file
    gene.append("".join(lines))
    #append the split lines to a list called gene

gff_parser = GFF_parser()
chromosome_list = gff_parser.chromosome(GFF_file)
gene_names = gff_parser.gene_names(chromosome_list,GFF_file)
clusters = patterns()
patternofclusters = clusters.chrom(gene_names,gene)

bed_data = []
hg_data = []
cluster_info1 = []
data = []

counter = -1

for chromosome in range(len(chromosome_list)):
    print("_______________________________________________________________________________________________")
    print("chromosome:",chromosome_list[chromosome])
    findclusters = HPFindClusters()
    #initialize the cluster called HPFindClusters
    combinations = HPCombinations()
    combination = Combinations()
    #initialize the cluster called HPCombinations
    ball_positions = [i for i, x in enumerate(patternofclusters[chromosome]) if x == 1]
    #append the index of the upregulated genes to a list called ball positions
    print(ball_positions)
    number_balls = len(ball_positions)
    #number of balls represents the number of upregulated genes per chromosome
    print("number balls:",number_balls)
    total_length = len(gene_names[chromosome])
    #total length represents the number of genes per chromosome
    print("total length:",total_length)
    #number of gene distance which upregulated genes must be to be considered a cluster
    number_of_clusters = findclusters.FindClusters(total_length,ball_positions,number_balls,minimum_window)
    print("number of clusters:",number_of_clusters)
    if number_of_clusters >= 1:

        print("===============================================================================================")
        #if there is a gene cluster present on chromosome
        chromosome_length = total_length
        #chromosome length is equal to total number of genes on chromosome
        print("chromosome length:", chromosome_length)
        total_genes = number_balls
        #number of responsive genes on chromosome
        print("total genes:",total_genes)
        genes_in_range = HPArrayOfGeneInfo()
        cluster = genes_in_range.GenesInRange(ball_positions,number_balls,minimum_window)
        #cluster represents list of list that contains the index of clustered responsive genes that are within minimum window:5
        print("cluster:",cluster)
        for window in cluster:
            print("***********************************************************************************************")
            counter += 1
            print("cluster no.", counter+1)
            genes = []
            names = []
            print("gene positions:",window)
            size = window[-1] - window[0] + 1
            #size represents the total number of genes in gene window.
            print("window size:",size)
            genes_in_window = len(window)
            #genes_in_window represents the number of responsive genes in the window
            print("genes in window:",genes_in_window)
            for position in window:
            #adds information of the individual genes
                genes.append(gene_names[chromosome][position])
                bed_data.append(gene_names[chromosome][position])
                gene = GFF_parser()
                gene_region = gene.gene_region(chromosome_list,GFF_file)
                names.append(gene_region[chromosome][position])
                if chromosome_length < 1000:
                    significance = combinations.ClusterSignificance(chromosome_length, total_genes, size, genes_in_window)
                #This applies to larger genomes with larger GFF files that has >1000 ORF on a chromosome. Uses the Combination class with Decimal datatype
                else:
                    significance = combination.ClusterSignificance(chromosome_length, total_genes, size, genes_in_window)
                hg_data.append(significance)
            data.append(significance)
            words = (names[0][0],names[-1][1])
            print("gene names:",genes)
            for i in bed_data:
                cluster_info = [chromosome_list[chromosome],words[0],words[1],data[counter]]
            cluster_info1.append(cluster_info)            
            print("significance:","{:.2E}".format(significance))
            print("***********************************************************************************************")


extract_list = []
#This part of the code collects information of the different clusters and also individual ORFs to be able to write to a text file
#It also checks if the user applied a false discovery rate and appends that p-value accordingly
if args.w == "None":
    sig_data = list(zip(bed_data, hg_data))
    collect_info = list(zip(cluster_info1,data))            
    gt = [list(i) for i in collect_info]
    for i in gt:
        lin = [i[0][0],i[0][1],i[0][2],i[0][3],i[1]]
        extract_list.append(lin)
    
else:
    if len(data) == 0:
        print("No clusters were detected")
        print("Could be possible errors, check input files or parameters")
        exit()
    else:
        pass
    p_adjusted = multipletests(data, method=args.w)
    new_pvalues = list(p_adjusted[1])
    cluster_newpvalues = list(zip(cluster_info1, new_pvalues))
    list_gene_sig = [list(i) for i in cluster_newpvalues]
    for i in list_gene_sig:
        extract_gene_sig = [i[0][0],i[0][1],i[0][2],i[0][3],i[1]]
        extract_list.append(extract_gene_sig)

    collect_info_list = []
    if type_gene == "":
        with open(GFF_file, "r") as f:
            for line in f:
                if str(idstring) in line:
                    for i in bed_data:
                        split_line = line.split()
                        searched_line = split_line[8]                    
                        split_searched_line = [(re.split(r'(:+|;+|=)', searched_line))]
                        for y in split_searched_line:
                            if y[4] == i:
                                info_split_searched_line = [split_line[0],int(split_line[3]),int(split_line[4]),split_line[8]]
                                collect_info_list.append(info_split_searched_line)
    if type_gene == "Name=":
        with open(GFF_file, "r") as f:
            for line in f:
                if str(type_gene) in line:
                    if str(idstring) in line:
                        for i in bed_data:
                            split_line = line.split()
                            searched_line = split_line[8]                       
                            split_searched_line = [(re.split(r'(:+|;+|=)', searched_line))]
                            for y in split_searched_line:
                                if y[4] == i:
                                    info_split_searched_line = [split_line[0],int(split_line[3]),int(split_line[4]),split_line[8]]
                                    collect_info_list.append(info_split_searched_line)

    cluster_info_list = []
    for i in collect_info_list:
        for k in extract_list:
            if i[0] == k[0] and i[1] == k[1]:
                found_info = [i[0],i[1],i[2],i[3],k[3],k[4]]
                cluster_info_list.append(found_info)

            if i[0] == k[0] and i[2] == k[2]:
                u = [i[0],i[1],i[2],i[3],k[3],k[4]]
                cluster_info_list.append(found_info)

            if i[0] == k[0] and i[1] > k[1] and i[2] < k[2]:
                found_info = [i[0],i[1],i[2],i[3],k[3],k[4]]
                cluster_info_list.append(found_info)


    adjusted_value = []
    for i in cluster_info_list:
        p = [i[5]]
        adjusted_value.append(p)
        flattened = [val for sublist in adjusted_value for val in sublist]
        sig_data = list(zip(bed_data, flattened))

if args.w == "None":
    print("_______________________________________________________________________________________________")
            
else:
    print("_______________________________________________________________________________________________")
    print("Correction test:",args.w)
    for i in extract_list:
        print("chromosome:",i[0])
        print("region:",i[1],"-",i[2])
        print("corrected p-value:",i[4])
        print("_______________________________________________________________________________________________")

print("===============================================================================================")

def flatten(list):
  for o in list:
    for r in o:
      yield r


def myround(x, base=10):
    return base * round(x/base)

#This part of the code writes the information of the respective gene clusters and corresponding ORFs to two separate BED files
#This also integrates a colour gradient for the p-value as RGB value as seen by the log scale in the code (Red colour gradient: 255,0,0- most significant p-value)

if len(sig_data) == 0:
    print("No clusters detected")
    exit()

def bed_builder(GFF_file,sig_data):
    fl = 'track name= %s description="Item RGB demonstration" itemRgb="On"' %args.s
    filename = "%s_genes.bed" %args.s
    k = open(filename,"w")
    k.write(fl + '\n')
    k.close()
    k = open(filename,"a+")
    BED = []

    for m in sig_data:
        if m[1] <= low_sig_value:
            colour = str("255,0,0")
            if type_gene == "":
                with open(GFF_file, "r") as f:
                    for line in f:
                        if str(idstring) in line:
                            split_line_1 = line.split()
                            element_in_line = split_line_1[8]
                            char_split = [(re.split(r'(:+|;+|=)', element_in_line))]                        
                            for y in char_split:
                                if y[4] == m[0]:
                                    lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[colour]]
                                    new_lines = flatten(lines)
                                    BED.append(new_lines)

            if type_gene == "Name=":
                with open(GFF_file, "r") as f:
                    for line in f:
                        if str(type_gene) in line:
                            if str(idstring) in line:
                                split_line_1 = line.split()
                                element_in_line = split_line_1[8]
                                char_split = [(re.split(r'(:+|;+|=)', element_in_line))]
                                for y in char_split:
                                    if y[4] == m[0]:
                                        lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[colour]]
                                        new_lines = flatten(lines)
                                        BED.append(new_lines)
                            
                                
        if m[1] > low_sig_value and m[1] < high_sig_value:
            if type_gene == "":
                with open(GFF_file, "r") as f:
                    for line in f:
                        if str(idstring) in line:
                            split_line_1 = line.split()
                            element_in_line = split_line_1[8]
                            char_split = [(re.split(r'(:+|;+|=)', element_in_line))]
                        
                            for y in char_split:
                                if y[4] == m[0]:
                                    x = (np.log(m[1])+2.4740649358032383)/-0.007670633343976754
                                    if x >= 250:
                                        colour = str("250,0,0")
                                        lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[colour]]
                                        new_lines = flatten(lines)
                                        BED.append(new_lines)
                                    if x <250:
                                        x = myround(x,10)
                                        colour = [int(x),0,0]
                                        c1 = ','.join(map(str, colour))
                                        lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[c1]]
                                        new_lines = flatten(lines)
                                        BED.append(new_lines)

            if type_gene == "Name=":
                with open(GFF_file, "r") as f:
                    for line in f:
                        if str(type_gene) in line:
                            if str(idstring) in line:
                                split_line_1 = line.split()
                                element_in_line = split_line_1[8]
                                char_split = [(re.split(r'(:+|;+|=)', element_in_line))]
                        
                                for y in char_split:
                                    if y[4] == m[0]:
                                        x = (np.log(m[1])+2.4740649358032383)/-0.007670633343976754
                                        if x >= 250:
                                            colour = str("250,0,0")
                                            lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[colour]]
                                            new_lines = flatten(lines)
                                            BED.append(new_lines)
                                        if x <250:
                                            x = myround(x,10)
                                            colour = [int(x),0,0]
                                            c1 = ','.join(map(str, colour))
                                            lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[c1]]
                                            new_lines = flatten(lines)
                                            BED.append(new_lines)

        if m[1] > high_sig_value:
            continue
    bedbuilder = csv.writer(k, delimiter='\t')
    bedbuilder.writerows(BED)
    
    return(bedbuilder)

bedbuilder = bed_builder(GFF_file, sig_data)

def clusters(extract_list):
    fl = 'track name= %s_clusters description="Item RGB demonstration" itemRgb="On"' %args.s
    filename = "%s_clusters.bed" %args.s
    k = open(filename,"w")
    k.write(fl + '\n')
    k.close()
    k = open(filename,"a+")
    BD = []
    if args.w == "None":
        for i in extract_list:
            if i[3] <= low_sig_value:
                colour = str("255,0,0")
                BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,colour]
                BD.append(BED_info)

            if i[3] > low_sig_value and i[3] < high_sig_value:
                x = (np.log(i[3])+2.4740649358032383)/-0.007670633343976754
                
                if x >= 250:
                    colour = str("250,0,0")
                    BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,colour]
                    BD.append(BED_info)
                if x < 250:
                    x = myround(x,10)
                    colour = [int(x),0,0]
                    c1 = ','.join(map(str, colour))
                    BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,c1]
                    BD.append(BED_info)
            if i[3] > high_sig_value:
                continue

    else:
        for i in extract_list:
            if i[4] <= low_sig_value:
                colour = str("255,0,0")
                BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,colour]
                BD.append(BED_info)

            if i[4] > low_sig_value and i[4] < high_sig_value:
                x = (np.log(i[4])+2.4740649358032383)/-0.007670633343976754
                
                if x >= 250:
                    colour = str("250,0,0")
                    BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,colour]
                    BD.append(BED_info)
                if x < 250:
                    x = myround(x,10)
                    colour = [int(x),0,0]
                    c1 = ','.join(map(str, colour))
                    BED_info = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,c1]
                    BD.append(BED_info)
            if i[4] > high_sig_value:
                continue
    bedbuilder = csv.writer(k, delimiter='\t')
    bedbuilder.writerows(BD)
    return(bedbuilder)


clusters(extract_list)
