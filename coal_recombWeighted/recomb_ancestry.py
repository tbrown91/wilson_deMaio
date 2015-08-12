def split_ancestry(child_ancestry,recomb_point):
    parent1_ancestry = []
    parent2_ancestry = []

    for i in range(len(child_ancestry)):
        if recomb_point > child_ancestry[i][1]:
            parent1_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
        elif recomb_point <= child_ancestry[i][0]:
            parent2_ancestry.append([child_ancestry[i][0],child_ancestry[i][1]])
        else:
            parent1_ancestry.append([child_ancestry[i][0],recomb_point-1])
            parent2_ancestry.append([recomb_point,child_ancestry[i][1]])

    return [parent1_ancestry,parent2_ancestry]
#END SPLIT_ANCESTRY()
