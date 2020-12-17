class HPFindClusters():
    #detects the number of clusters

    def __init__(self):
        self.number_of_clusters = 0  
    
    def FindClusters(self, total_length, ball_positions, number_balls, minimum_window):
        self.number_of_clusters = 0
        start = 0
        begin = start
        cluster_number = 0
        count = 0
        while start < number_balls - 1 :
            if ball_positions[start+1] - ball_positions[start] <= minimum_window:  
                if start == begin:    
                    cluster_number = cluster_number +1        
                    self.number_of_clusters = cluster_number   
                else:
                    self.number_of_clusters = cluster_number    
            else: 
                begin = start + 1    
            pass
            start = start + 1  
        return(self.number_of_clusters)
