import numpy as np
import random

def recomb_interval(ancestral_intervals,interval_starts,genome_length,break_rate):
    #Choose a start and end site for the recombination interval conditioned
    #on the fact that the recombinant interval must split the ancestral material

    if len(interval_starts) == 0:
        #No break points, therefore recombinant start is equally weighted in interval
        start_site = random.sample(range(ancestral_intervals[0][0],ancestral_intervals[0][1]+1),1)
        start_site = start_site[0]

        #Sample end point via a truncated geometric distribution - must be less than L
        interval_rand = np.random.rand(1)
        length = np.floor(np.log(1 - interval_rand[0] * (1 - (1-break_rate)**(genome_length-1)))/np.log(1-break_rate))
        end_site = (start_site + length) % genome_length

    else:
        #Start of interval is weighted 1/break_rate more than any nucleotide
        #inside an ancestral interal
        #Calculate the total amount of ancestral material in the node
        total_material = 0
        for i in range(len(ancestral_intervals)):
            total_material = total_material + ancestral_intervals[i][1] - ancestral_intervals[i][0] + 1
        start_rand = np.random.rand(1)

        if start_rand[0] < (len(interval_starts)/break_rate)/(total_material - len(interval_starts) + (len(interval_starts)/break_rate)):
            #Recombination starts at an interval start site

            if (ancestral_intervals[0][0] == 0) and (ancestral_intervals[-1][1] == genome_length - 1):
                #Account for ancestral material that wraps around the end of the genome
                start_index = random.sample(range(len(interval_starts)),1)
                start_site = ancestral_intervals[start_index[0] + 1][0]
                interval_rand = np.random.rand(1)
                #Sample end site conditional on it splitting the ancestral material - must be less than L - (start_i - end_i-1)
                length = np.floor(np.log(1 - interval_rand[0] * (1 - (1-break_rate)**(genome_length - ancestral_intervals[start_index[0]+1][0] + ancestral_intervals[start_index[0]][1])))/np.log(1-break_rate))
                end_site = (start_site + length) % genome_length

            else:
                #All intervals are contained within the genome
                start_index = random.sample(range(len(interval_starts)),1)
                start_site = ancestral_intervals[start_index[0]][0]

                if start_index[0] == 0:
                    interval_rand = np.random.rand(1)
                    #Sample end site conditional on it splitting the ancestral material - must be less than L - (start_i - end_i-1)
                    length = np.floor(np.log(1 - interval_rand[0] * (1 - (1-break_rate)**(ancestral_intervals[-1][1] - ancestral_intervals[0][0])))/np.log(1-break_rate))
                    end_site = (start_site + length) % genome_length

                else:
                    interval_rand = np.random.rand(1)
                    #Sample end site conditional on it splitting the ancestral material - must be less than L - (start_i - end_i-1)
                    length = np.floor(np.log(1 - interval_rand[0] * (1 - (1-break_rate)**(genome_length - ancestral_intervals[start_index[0]][0] + ancestral_intervals[start_index[0]-1][1])))/np.log(1-break_rate))
                    end_site = (start_site + length) % genome_length

        else:
            #Recombination starts inside an ancestral interval

            if (ancestral_intervals[0][0] == 0) and (ancestral_intervals[-1][1] == genome_length - 1):
                #Last interval wraps around end of genome - therefore include first and last nucleotides of genome in sample
                start_site = random.sample(range(int(total_material) - len(interval_starts)),1)[0]
                #Find nucleotide which corresponds with chosen element of ancestral material
                for i in range(len(ancestral_intervals)-1):
                    if start_site > ancestral_intervals[i][1]:
                        start_site = start_site - ancestral_intervals[i][1]
                        start_site = start_site + ancestral_intervals[i+1][0]

                    else:
                        break

            else:
                start_site = random.sample(range(int(total_material) - len(interval_starts)),1)
                start_site = start_site[0] + ancestral_intervals[0][0] + 1
                for i in range(len(ancestral_intervals)-1):

                    if start_site > ancestral_intervals[i][1]:
                        start_site = start_site - ancestral_intervals[i][1]
                        start_site = start_site + ancestral_intervals[i+1][0]

                    else:
                        break

            #Sample end point via a truncated geometric distribution - must be less than L
            interval_rand = np.random.rand(1)
            length = np.floor(np.log(1 - interval_rand[0] * (1 - (1-break_rate)**(genome_length-1)))/np.log(1-break_rate))
            end_site = (start_site + length) % genome_length

    return [int(start_site),int(end_site)]
#END RECOMB_INTERVAL
