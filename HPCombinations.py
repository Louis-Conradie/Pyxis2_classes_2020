import numpy as np

class HPCombinations():

    def __init__(self):
        self.s = 0
    def Factorial(self, start, end):
        result = np.longdouble(1)
        for i in range(1,end+1):
            result = result*i
        return(result)

    def Combinations(self, n, r):
        C_n_r = self.Factorial(1,(n-r))
        C_r = self.Factorial(1,(r))
        C_n = self.Factorial(1,(n))
        result = np.longdouble(C_n/(C_n_r*C_r))
        return(result)

    def ClusterSignificance(self, total_length, total_responsive_genes, window_size, genes_in_window):
        window_combinations = np.longdouble(self.Combinations(window_size, genes_in_window))
        chromosome_combinations = np.longdouble(self.Combinations(total_length, total_responsive_genes))
        combinations_of_balls_outside_window = np.longdouble(self.Combinations((total_length - window_size),(total_responsive_genes - genes_in_window)))
        result = np.longdouble(combinations_of_balls_outside_window * (window_combinations/chromosome_combinations))
        return(result)
