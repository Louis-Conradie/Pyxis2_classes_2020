
# coding: utf-8

from GFF_parser import GFF_parser
from HPArrayOfGeneInfo import HPArrayOfGeneInfo
from HPCombinations import HPCombinations
from HPFindClusters import HPFindClusters
from Combination import Combinations
from Patterns import patterns
from statsmodels.sandbox.stats.multicomp import multipletests
##RNDW
import csv
import argparse
import numpy as np
import itertools
import matplotlib.pyplot as plt
import re
##RNDW
# In[28]:
parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')
parser.add_argument('-v', metavar='v', type=str, help="GFF file") #GFF file
parser.add_argument('-s', metavar='s', type=str, help="List of regulated genes") #file with upregulated genes
parser.add_argument('-w', metavar='w', type=str, nargs='?', default='None', const='None', help='Statistical test used to detect clusters')
parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')
parser.add_argument('-k', metavar='k', type=str, nargs='?', default='Name=', const='Name=', help='Search string used to identify Gene name in GFF file (default: Name=)')
parser.add_argument('-p', metavar='p', type=int, nargs='?', default=0.01, const=0.01, help='P-value used to make statistical test more rigorous (default: 0.05)')
parser.add_argument('-l', metavar='l', type=int, nargs='?', default=5, const=5, help='Minimum window to detect gene clusters (default: 5)')
args = parser.parse_args()
print("The search string for the GFF file is: ", args.i)
print("The database string for the GFF file is: ", args.k)
print("The cut-off p-value is: ", args.p)
print("The minimum window for the detection of gene clusters is: ", args.l)
GFF_file = args.v
files = args.s
idstring = args.i
sig_value = args.p
minimum_window = args.l
##RNDW

#files = file that contains the names of upregulated genes
##LTC

#creates accumulator list called gene
#files = file that contains the names of upregulated genes
##LTC
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

##LTC
##rndw
bed_data = []
hg_data = []
significant_genes = []
unsignificant_genes = []
g = []
bd = []
data = []
w = []
j = []

counter = -1
##rndw
##LTC
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
            counter += 1
            genes = []
            names = []
            print("gene cluster:",window)
            size = window[-1] - window[0] + 1
            #size represents the total number of genes in gene window.
            print("window size:",size)
            genes_in_window = len(window)
            #genes_in_window represents the number of responsive genes in the window
            print("genes in window:",genes_in_window)
            for position in window:
                genes.append(gene_names[chromosome][position])
                bed_data.append(gene_names[chromosome][position])
                gene = GFF_parser()
                gene_region = gene.gene_region(chromosome_list,GFF_file)
                names.append(gene_region[chromosome][position])
                if chromosome_length < 1000:
                    significance = combinations.ClusterSignificance(chromosome_length, total_genes, size, genes_in_window)
                else:
                    significance = combination.ClusterSignificance(chromosome_length, total_genes, size, genes_in_window)
                hg_data.append(significance)
            data.append(significance)


            words = (names[0][0],names[-1][1])
            g.append(words)
            print("gene names:",genes)
            for i in bed_data:
                lk = [chromosome_list[chromosome],words[0],words[1],data[counter]]
            bd.append(lk)
            
            bs = []
            r = []
            print("significance:","{:.2E}".format(significance))
            if args.w == "None":
                sig_data = list(zip(bed_data, hg_data))
                lis = list(zip(bd,data))
             
                k = [list(i) for i in lis]
                for i in k:
                    lin = [i[0][0],i[0][1],i[0][2],i[0][3],i[1]]
                    r.append(lin)
                
            else:
                p_adjusted = multipletests(data, method=args.w)
                s = list(p_adjusted[1])
                lis = list(zip(bd, s))
                k = [list(i) for i in lis]
                for i in k:
                    lin = [i[0][0],i[0][1],i[0][2],i[0][3],i[1]]
                    r.append(lin)
                
                l = []
                with open(GFF_file, "r") as f:
                    for line in f:
                        for i in bed_data:
                            if str('ID=gene:') + i + str(';Name=') in line:
                                lines = line.split()
                                q = [lines[0],int(lines[3]),int(lines[4]),lines[8]]
                                l.append(q)
                
                t = []
                for i in l:
                    for k in r:
                        if i[0] == k[0] and i[1] == k[1]:
                            u = [i[0],i[1],i[2],i[3],k[3],k[4]]
                            t.append(u)
                        if i[0] == k[0] and i[2] == k[2]:
                            u = [i[0],i[1],i[2],i[3],k[3],k[4]]
                            t.append(u)
                        if i[0] == k[0] and i[1] > k[1] and i[2] < k[2]:
                    
                            t.append(u)

                y = []
                for i in t:
                    p = [i[5]]
                    y.append(p)
                    flattened = [val for sublist in y for val in sublist]
                    sig_data = list(zip(bed_data, flattened))
 
print("===============================================================================================")



##RNDW
def flatten(list):
  for o in list:
    for r in o:
      yield r


def myround(x, base=10):
    return base * round(x/base)

def bed_builder(GFF_file,sig_data):
    fl = 'track name= %s description="Item RGB demonstration" itemRgb="On"' %args.s
    filename = "%s_genes.bed" %args.s
    k = open(filename,"w")
    k.write(fl + '\n')
    k.close()
    k = open(filename,"a+")
    BED = []
    for m in sig_data:
        if m[1] < args.p:
            colour = str("255,0,0")
            with open(GFF_file, "r") as f:
                for line in f:
                    if str('ID=gene:') + m[0] + str(';Name=') in line:
                        lines = [np.array(str.split(line))[[0,3,4]],[m[0]],["{:.2E}".format(m[1])],[0],[0],[0],[colour]]
                        new_lines = flatten(lines)
                        BED.append(new_lines)
        if m[1] > args.p and m[1] < args.p+0.04:
            with open(GFF_file, "r") as f:
                for line in f:
                    if str('ID=gene:') + m[0] + str(';Name=') in line:
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
        if m[1] > args.p+0.04:
            continue
    bedbuilder = csv.writer(k, delimiter='\t')
    bedbuilder.writerows(BED)
    return(bedbuilder)

bedbuilder = bed_builder(GFF_file, sig_data)
#rt(h,new,new_list, new_list2)

def clusters(r):
    fl = 'track name= %s_clusters description="Item RGB demonstration" itemRgb="On"' %args.s
    filename = "%s_clusters.bed" %args.s
    k = open(filename,"w")
    k.write(fl + '\n')
    k.close()
    k = open(filename,"a+")
    BD = []
    if args.w == "None":
        for i in r:
            if i[3] <= args.p:
                colour = str("255,0,0")
                lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,colour]
                BD.append(lin)

            if i[3] > args.p and i[3] < args.p+0.04:
                x = (np.log(i[4])+2.4740649358032383)/-0.007670633343976754
                
                if x >= 250:
                    colour = str("250,0,0")
                    lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,colour]
                    BD.append(lin)
                if x < 250:
                    x = myround(x,10)
                    colour = [int(x),0,0]
                    c1 = ','.join(map(str, colour))
                    lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[3]),0,0,0,c1]
                    BD.append(lin)
            if i[3] > args.p+0.04:
                continue
    else:
        for i in r:
            if i[4] <= args.p:
                colour = str("255,0,0")
                lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,colour]
                BD.append(lin)

            if i[4] > args.p and i[4] < args.p+0.04:
                x = (np.log(i[4])+2.4740649358032383)/-0.007670633343976754
                
                if x >= 250:
                    colour = str("250,0,0")
                    lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,colour]
                    BD.append(lin)
                if x < 250:
                    x = myround(x,10)
                    colour = [int(x),0,0]
                    c1 = ','.join(map(str, colour))
                    lin = [i[0],i[1],i[2],"cluster","{:.2E}".format(i[4]),0,0,0,c1]
                    BD.append(lin)
            if i[4] > args.p+0.04:
                continue

    bedbuilder = csv.writer(k, delimiter='\t')
    bedbuilder.writerows(BD)
    return(bedbuilder)


clusters(r)

