class patterns():
    
    def __init__(self):
        self.s = 0
    
    def chrom(self, gene_names, gene):
        count = len(gene_names)
        #creates integer that represents number of chromosomes(counts number of lists in gene names)
        patternofclusters = [[] for i in range(0, count)]
        #creates empty nested list "patternofclusters" according to the number of chromosomes

        for number in range(0,count):
            for i in range(len(gene_names[number])):
                i = 0
                patternofclusters[number].append(i)
                #append 0's to nested list patternofclusters that represents balls(genes) on a stick(chromosome)
            for name in gene:
                #gene represents the list of upregulated genes so for element in the list of upregulated genes
                try:
                    position = gene_names[number].index(name)
                    #creates variable that finds index for each upregulated gene in nested list gene_names 
                    patternofclusters[number][position] = 1
                    #add integer 1 to index in nested list patternofclusters 
                except ValueError:
                    #if there is any index errors because gene is not found in gene names let code continue
                      continue
        return(patternofclusters)
        #return the nested list patternofclusters
