import re
from itertools import groupby
import argparse

parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')
parser.add_argument('-v', metavar='v', type=str, help="GFF file") #GFF file
parser.add_argument('-s', metavar='s', type=str, help="List of regulated genes") #file with upregulated genes
parser.add_argument('-w', metavar='w', type=str, nargs='?', default='None', const='None', help='Statistical test used to detect clusters')
parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')
parser.add_argument('-k', metavar='k', type=int, nargs='?', help='Search string used to identify Gene name in GFF file (default: Name=)')
parser.add_argument('-p', metavar='p', type=float, nargs='?', default=0.01, const=0.01, help='Lower P-value used to make statistical test more rigorous (default: 0.01)')
parser.add_argument('-q', metavar='q', type=float, nargs='?',default=0.05, const=0.05, help='Maximum p-value to make statistical test more rigorous (default:0.05)')
parser.add_argument('-l', metavar='l', type=int, nargs='?', default=5, const=5, help='Minimum window to detect gene clusters (default: 5)')
args = parser.parse_args()

if args.k == 1:
    type_gene = str("Name=")
    new_gene = str("Protein coding genes")
    
if args.k == 0:
    type_gene = str("")
    new_gene = str("Protein coding, dubious and merged ORFs")


GFF_file = args.v
files = args.s
idstring = args.i
sig_value = args.p
minimum_window = args.l

class GFF_parser():
    def __init__(self):
        self.chromosome_list = []
        self.sequence_region = []


    def chromosome(self,GFF_file):
        seq = ["##sequence-region"]
        # attribute that are used to identify the lines that contain the name of genes
        lines = []
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for att in seq:
                if att in line:
                    lines = line.split()
                    self.chromosome_list.append(lines[1])
                    # finds lines in the opened files that contains the lines with attribute
                    #split the lines into strings and appends it into list called "lines"
                    #from every list with appended string append the first element(string) to empty list "chromosome_list"
        return(self.chromosome_list)

    def sequence_region(self,GFF_file):
        seq = ["##sequence-region"]
        # attribute that are used to identify the lines that contain the name of genes
        lines = []
        gff_file = (GFF_file,"r")
        for line in gff_file:
            for att in seq:
                if att in line:
                    lines = line.split()
                    self.sequence_region.append(lines[2:4])
                    self.sequence_region = [tuple(int(x) for x in tup) for tup in self.sequence_region]
                    # finds lines in the opened files that contains the lines with attribute
                    #split the lines into strings and appends it into list called "lines"
                    #append second to third element to empty list sequence_region
                    #join every two elements into tuple and change the strings into integers
        return(self.sequence_region)

    def gene_region(self,chromosome_list,GFF_file):
        lines = []
        count = len(chromosome_list)
        gene_region = [[] for i in range(0,count)]
        #creates empty list of lists in range of count
        att = [args.i]
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for name in att:
                if name in line:
                    words = line.split()
                    if words[2] == "gene":
                        if type_gene in words[8]:
                            lines.append(words)
                        else:
                            if type_gene == "":
                                lines.append(words)
        header = [list(g) for k, g in groupby(lines, key=lambda x: x[0])]
        #create list of lists called header and group each lines into lists based on the first element (chromosome number) of the split lines
        genes = [[x[3:5] for x in y] for y in header]
        #create new list of lists called genes and append the third and fourth element from each split lines of nested list header
        for t in range(count):
            while t < count:
                gene_region[t] = [tuple(int(x) for x in tup) for tup in genes[t]]
                #create paired elements from list genes and create tuples and convert elements(strings) to integers
                pass
                t = t+1
        return(gene_region)

    def gene_names(self, chromosome_list, GFF_file):
        lines = []
        att = [args.i]
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for name in att:
                if name in line:
                    words = line.split()
                    if words[2] == "gene":
                        #if third element contains the name gene
                        if type_gene in words[8]:
                            lines.append(words)
                        else:
                            if type_gene == "":
                                lines.append(words)
                                
                        #append the split lines to empty list called lines
        names = []
        gene_names = []

        header = [list(g) for k, g in groupby(lines, key=lambda x: x[0])]
        #create list of lists called header and group each lines into lists based on the first element (chromosome number) of the split lines
        for lines in header:
            for lists in lines:
                for element in lists:
                    if args.i in element:
                        words = [lists[0] ,(re.split(r'(:+|;+|=)', element))]
                        #based on the characters you want to use to split the line containing "gene name"
                        names.append(words)
        for elements in names:
            #print(elements[1]) #- print this line and adjust where the gene name is located in the split line for your GFF file
            word = [elements[0],elements[1][4]]
            #GFF files differ based on the line containing the gene name so the x might differ in elements[1][x](elements[1][x-position of gene name])
            gene_names.append(word)
        gene_names = [list(g) for k, g in groupby(gene_names, key=lambda x: x[0])]

        gene_names = [[x[1] for x in y] for y in gene_names]

        return(gene_names)import re
from itertools import groupby
import argparse

parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')
parser.add_argument('-v', metavar='v', type=str, default="Example.GFF3", help="GFF file") #GFF file
parser.add_argument('-s', metavar='s', type=str, default="Example.txt", help="List of regulated genes") #file with upregulated genes
parser.add_argument('-w', metavar='w', type=str, nargs='?', default='None', const='None', help='Statistical test used to detect clusters')
parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')
parser.add_argument('-k', metavar='k', type=int, nargs='?', default=1, help='Search string used to identify Gene name in GFF file (default: 1)')
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


GFF_file = args.v
files = args.s
idstring = args.i
sig_value = args.p
minimum_window = args.l

class GFF_parser():
    def __init__(self):
        self.chromosome_list = []
        self.sequence_region = []


    def chromosome(self,GFF_file):
        seq = ["##sequence-region"]
        # attribute that are used to identify the lines that contain the name of genes
        lines = []
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for att in seq:
                if att in line:
                    lines = line.split()
                    self.chromosome_list.append(lines[1])
                    # finds lines in the opened files that contains the lines with attribute
                    #split the lines into strings and appends it into list called "lines"
                    #from every list with appended string append the first element(string) to empty list "chromosome_list"
        return(self.chromosome_list)

    def sequence_region(self,GFF_file):
        seq = ["##sequence-region"]
        # attribute that are used to identify the lines that contain the name of genes
        lines = []
        gff_file = (GFF_file,"r")
        for line in gff_file:
            for att in seq:
                if att in line:
                    lines = line.split()
                    self.sequence_region.append(lines[2:4])
                    self.sequence_region = [tuple(int(x) for x in tup) for tup in self.sequence_region]
                    # finds lines in the opened files that contains the lines with attribute
                    #split the lines into strings and appends it into list called "lines"
                    #append second to third element to empty list sequence_region
                    #join every two elements into tuple and change the strings into integers
        return(self.sequence_region)

    def gene_region(self,chromosome_list,GFF_file):
        lines = []
        count = len(chromosome_list)
        gene_region = [[] for i in range(0,count)]
        #creates empty list of lists in range of count
        att = [args.i]
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for name in att:
                if name in line:
                    words = line.split()
                    if words[2] == "gene":
                        if type_gene in words[8]:
                            lines.append(words)
                        else:
                            if type_gene == "":
                                lines.append(words)
        header = [list(g) for k, g in groupby(lines, key=lambda x: x[0])]
        #create list of lists called header and group each lines into lists based on the first element (chromosome number) of the split lines
        genes = [[x[3:5] for x in y] for y in header]
        #create new list of lists called genes and append the third and fourth element from each split lines of nested list header
        for t in range(count):
            while t < count:
                gene_region[t] = [tuple(int(x) for x in tup) for tup in genes[t]]
                #create paired elements from list genes and create tuples and convert elements(strings) to integers
                pass
                t = t+1
        return(gene_region)

    def gene_names(self, chromosome_list, GFF_file):
        lines = []
        att = [args.i]
        gff_file = open(GFF_file,"r")
        for line in gff_file:
            for name in att:
                if name in line:
                    words = line.split()
                    if words[2] == "gene":
                        #if third element contains the name gene
                        if type_gene in words[8]:
                            lines.append(words)
                        else:
                            if type_gene == "":
                                lines.append(words)
                                
                        #append the split lines to empty list called lines
        names = []
        gene_names = []

        header = [list(g) for k, g in groupby(lines, key=lambda x: x[0])]
        #create list of lists called header and group each lines into lists based on the first element (chromosome number) of the split lines
        for lines in header:
            for lists in lines:
                for element in lists:
                    if args.i in element:
                        words = [lists[0] ,(re.split(r'(:+|;+|=)', element))]
                        #based on the characters you want to use to split the line containing "gene name"
                        names.append(words)
        for elements in names:
            #print(elements[1]) #- print this line and adjust where the gene name is located in the split line for your GFF file
            word = [elements[0],elements[1][4]]
            #GFF files differ based on the line containing the gene name so the x might differ in elements[1][x](elements[1][x-position of gene name])
            gene_names.append(word)
        gene_names = [list(g) for k, g in groupby(gene_names, key=lambda x: x[0])]

        gene_names = [[x[1] for x in y] for y in gene_names]

        return(gene_names)
