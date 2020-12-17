import argparse

parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')
parser.add_argument('-v', metavar='v', type=str, help="GFF file") #GFF file
parser.add_argument('-s', metavar='s', type=str, help="List of regulated genes") #file with upregulated genes
parser.add_argument('-w', metavar='w', type=str, nargs='?', default='None', const='None', help='Statistical test used to detect clusters')
parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')
parser.add_argument('-k', metavar='k', type=int, nargs='?', help='Search string used to identify Gene name in GFF file (default: 1)')
parser.add_argument('-p', metavar='p', type=float, nargs='?', default=0.01, const=0.01, help='Lower P-value used to make statistical test more rigorous (default: 0.01)')
parser.add_argument('-q', metavar='q', type=float, nargs='?',default=0.05, const=0.05, help='Maximum p-value to make statistical test more rigorous (default:0.05)')
parser.add_argument('-l', metavar='l', type=int, nargs='?', default=5, const=5, help='Minimum window to detect gene clusters (default: 5)')
args = parser.parse_args()

GFF_file = args.v
files = args.s
idstring = args.i
sig_value = args.p
minimum_window = args.l

class HPArrayOfGeneInfo():

    def __init__(self):
        self.up_gene = []


    def GenesInRange(self, ball_positions, number_balls, minimum_window):
        start = 0
        begin = start

        while start < number_balls - 1 :
            try:
                if ball_positions[start+1] - ball_positions[start] <= minimum_window:
                    self.up_gene.append(ball_positions[start])
                    while ball_positions[start+1] - ball_positions[start] <= minimum_window:
                        self.up_gene.append(ball_positions[start+1])
                        pass
                        start =start + 1
                else: 
                    begin = start + 1

            except IndexError:
                break

            pass
            start = start + 1


        inds=[0]+[ind for ind,(i,j) in enumerate(zip(self.up_gene,self.up_gene[1:]),1) if j-i>minimum_window]+[len(self.up_gene)+1]
        up_genes = [self.up_gene[i:j] for i,j in zip(inds,inds[1:])]
        #split the elements in list up_gene if they are bigger than 5 to create list of lists up_genes
        #Each list represents a "cluster"
        return(up_genes)
