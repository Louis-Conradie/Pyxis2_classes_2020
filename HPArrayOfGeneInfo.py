#class that are imported into python
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


        inds=[0]+[ind for ind,(i,j) in enumerate(zip(self.up_gene,self.up_gene[1:]),1) if j-i>5]+[len(self.up_gene)+1]
        up_genes = [self.up_gene[i:j] for i,j in zip(inds,inds[1:])]
        #split the elements in list up_gene if they are bigger than 5 to create list of lists up_genes
        #Each list represents a "cluster"
        return(up_genes)
