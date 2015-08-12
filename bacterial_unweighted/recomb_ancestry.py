def split_ancestry(child_ancestry,recomb_start,recomb_end,genome_length):
    parent1_ancestry = []
    parent2_ancestry = []

    if recomb_start <= recomb_end:
        for i in range(len(child_ancestry)):
            if (recomb_start <= child_ancestry[i][1]) and (recomb_end >= child_ancestry[i][0]):
                if (recomb_start <= child_ancestry[i][0]) and (recomb_end >= child_ancestry[i][1]):
                    parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
                elif (recomb_start <= child_ancestry[i][0]) and (recomb_end < child_ancestry[i][1]):
                    parent2_ancestry.append([child_ancestry[i][0],recomb_end])
                    parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])
                elif (recomb_start > child_ancestry[i][0]) and (recomb_end >= child_ancestry[i][1]):
                    parent2_ancestry.append([recomb_start,child_ancestry[i][1]])
                    parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])
                else:
                    parent2_ancestry.append([recomb_start,recomb_end])
                    parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])
                    parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])
            else:
                parent1_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
    else:
        if (child_ancestry[0][0] == 0) and (child_ancestry[0][1] == genome_length - 1):
            parent1_ancestry.append([recomb_end+1,recomb_start-1])
            parent2_ancestry.append([0,recomb_end])
            parent2_ancestry.append([recomb_start,genome_length-1])
        else:
            for i in range(len(child_ancestry)):
                if recomb_start <= child_ancestry[i][1]:
                    if recomb_start <= child_ancestry[i][0]:
                        parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
                    else:
                        parent2_ancestry.append([recomb_start,child_ancestry[i][1]])
                        parent1_ancestry.append([child_ancestry[i][0],recomb_start-1])
                elif recomb_end >= child_ancestry[i][0]:
                    if recomb_end >= child_ancestry[i][1]:
                        parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
                    else:
                        parent2_ancestry.append([child_ancestry[i][0],recomb_end])
                        parent1_ancestry.append([recomb_end+1,child_ancestry[i][1]])
                else:
                    parent1_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])


    return [parent1_ancestry,parent2_ancestry]
#END SPLIT_ANCESTRY()
