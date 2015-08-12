def split_ancestry(child_ancestry,recomb_start,recomb_end,genome_length):
    #Take the ancestral material from the chosen child and split the material
    #based on the recombination interval
    parent1_ancestry = []
    parent2_ancestry = []

    if recomb_start <= recomb_end:
        #Recombination interval does not wrap round the end of the genome

        for i in range(len(child_ancestry)):

            if (recomb_start <= child_ancestry[i][1]) and (recomb_end >= child_ancestry[i][0]):
                #Recombination interval contains some material from the current interval

                if (recomb_start <= child_ancestry[i][0]) and (recomb_end >= child_ancestry[i][1]):
                    #Recombination interval contains the entire interval
                    parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])

                elif (recomb_start <= child_ancestry[i][0]) and (recomb_end < child_ancestry[i][1]):
                    #Recombination interval takes the start of an ancestral interval - give to parent 2
                    parent2_ancestry.append([child_ancestry[i][0],recomb_end])
                    parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])

                elif (recomb_start > child_ancestry[i][0]) and (recomb_end >= child_ancestry[i][1]):
                    #Recombination interval takes the end of an ancestral interval - give to parent 2
                    parent2_ancestry.append([recomb_start,child_ancestry[i][1]])
                    parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])

                else:
                    #Recombination interval takes a central section of the interval
                    parent2_ancestry.append([recomb_start,recomb_end])
                    parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])
                    parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])

            else:
                #Recombination interval does not affect current ancestral interval
                parent1_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])

    else:
        #Recombination interval wraps around the end of the genome

        if (child_ancestry[0][0] == 0) and (child_ancestry[0][1] == genome_length - 1):
            #Ancestral material is entire genome
            parent1_ancestry.append([recomb_end+1,recomb_start-1])
            parent2_ancestry.append([0,recomb_end])
            parent2_ancestry.append([recomb_start,genome_length-1])

        else:
            for i in range(len(child_ancestry)):

                if recomb_start <= child_ancestry[i][1]:
                    #Recombination interval includes ancestral interval

                    if recomb_start <= child_ancestry[i][0]:
                        #Recombination interval includes entire ancestral interval
                        parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])

                    else:
                        #Recombination interval splits ancestral interval
                        parent2_ancestry.append([recomb_start,child_ancestry[i][1]])
                        parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])

                elif recomb_end >= child_ancestry[i][0]:
                    #Recombination interval includes ancestral interval

                    if recomb_end >= child_ancestry[i][1]:
                        #Recombination interval includes entire ancestral interval
                        parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])

                    else:
                        #Recombination interval splits ancestral interval
                        parent2_ancestry.append([child_ancestry[i][0],recomb_end])
                        parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])

                else:
                    #Recombination interval does not affect current ancestral interval
                    parent1_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])


    return [parent1_ancestry,parent2_ancestry]
#END SPLIT_ANCESTRY()

def recomb_intervalStarts(parent1_ancestry, parent2_ancestry,genome_length):
    #Count the number of interval start sites at each node
    parent1_starts = []
    parent2_starts = []

    #Parent 1 start sites
    if (parent1_ancestry[0][0] != 0) or (parent1_ancestry[-1][1] != genome_length - 1):
        #Last interval does not wrap around the end of the genome
        parent1_starts.append(parent1_ancestry[0][0])

    for i in range(1,len(parent1_ancestry)):
        #Count the start site of each interval
        parent1_starts.append(parent1_ancestry[i][0])

    #Parent 2 start sites
    if (parent2_ancestry[0][0] != 0) or (parent2_ancestry[-1][1] != genome_length - 1):
        #Last interval wraps around the end of the genome
        parent2_starts.append(parent2_ancestry[0][0])

    for i in range(1,len(parent2_ancestry)):
        #Count the start site of each interval
        parent2_starts.append(parent2_ancestry[i][0])

    return [parent1_starts,parent2_starts]
#END RECOMB_INTERVALSTARTS
