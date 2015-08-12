def create_initialList(data,nucleotide):
    #Create list of nodes which contain the nucleotide of interest in their ancestral material
    list = []
    for i in range(len(data)):
        #Check if node is ancestral

        for j in range(len(data[i].ancestral_material)):

            if (nucleotide >= data[i].ancestral_material[j][0]) and (nucleotide <= data[i].ancestral_material[j][1]):
                list.append(i)
                break

    return list
#END CREATE_INITIALLIST

def remove_unancestral(data,list):
    #Remove any children or parents that do not contain the nucleotide of interest

    for i in list:
        a = []
        #Only include parents if they contain the nucleotide in their ancestry
        for j in range(len(data[i].parent)):
            if data[i].parent[j] in list:
                a.append(data[i].parent[j])
        data[i].parent = a[:]

        a = []
        #Only include children if they contain the nucleotide in their ancestry
        for j in range(len(data[i].children)):
            if data[i].children[j] in list:
                a.append(data[i].children[j])
        data[i].children = a[:]

    return list
#END REMOVE_UNANCESTRAL

def remove_recomb(data,list):
    #Remove any nodes that are a result of recombination
    i = 0
    while i < len(list)-1:

        if len(data[list[i]].children) == 1:
            #Recombinant node
            #Update child node with new parent
            k = list[i]
            data[data[list[i]].children[0]].parent = []
            data[data[list[i]].children[0]].parent.append(data[list[i]].parent[0])

            #Update children of parent of the recombinant node with new child
            a = []
            for j in range(len(data[data[list[i]].parent[0]].children)):
                if data[data[list[i]].parent[0]].children[j] != k:
                    a.append(data[data[list[i]].parent[0]].children[j])
            a.append(data[list[i]].children[0])
            data[data[list[i]].parent[0]].children = []
            data[data[list[i]].parent[0]].children = a[:]

            #Remove recombinant node from the list of ancestral nodes
            a = []
            for j in list:
                if j != k:
                    a.append(j)
            list = a[:]
            #Reset node back to first leaves
            i = -1
        #END IF
        i = i+1
    #END WHILE

    return [data,list]
#END REMOVE_RECOMB

def create_newick(data,list,leaves):
    #Create the newick tree for the nucleotide of interest

    #Create the newick tree at each node above the leaves
    for i in list[leaves:]:
        parent = i

        if len(data[parent].label) == 0:
            child_0 = data[parent].children[0]
            child_1 = data[parent].children[1]
            # Create new label for parent node: (A:age_A,B:age_B)
            data[parent].label = "(" + data[child_0].label + ":" + str(data[parent].age - data[child_0].age)[1:-1] + "," + data[child_1].label + ":" + str(data[parent].age - data[child_1].age)[1:-1] + ")"

    return data[list[-1]].label
#END CREATE_NEWICK
