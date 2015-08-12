def create_initialList(data,nucleotide):
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
    for i in list:
        a = []
        for j in range(len(data[i].parent)):
            if data[i].parent[j] in list:
                a.append(data[i].parent[j])
        data[i].parent = a[:]
        a = []
        for j in range(len(data[i].children)):
            if data[i].children[j] in list:
                a.append(data[i].children[j])
        data[i].children = a[:]
    return list
#END REMOVE_UNANCESTRAL

def remove_recomb(data,list):
    i = 0
    while i < len(list)-1:
        if len(data[list[i]].children) == 1:
            #Recombinant node
            #Update child and parent node and remove from list
            k = list[i]
            data[data[list[i]].children[0]].parent = []
            data[data[list[i]].children[0]].parent.append(data[list[i]].parent[0])

            a = []
            for j in range(len(data[data[list[i]].parent[0]].children)):
                if data[data[list[i]].parent[0]].children[j] != k:
                    a.append(data[data[list[i]].parent[0]].children[j])
            a.append(data[list[i]].children[0])
            data[data[list[i]].parent[0]].children = []
            data[data[list[i]].parent[0]].children = a[:]

            a = []
            for j in list:
                if j != k:
                    a.append(j)
            list = a[:]
            i = -1
        i = i+1

    return [data,list]
#END REMOVE_RECOMB

def create_newick(data,list,leaves):
    for i in list[leaves:]:
        parent = i
        if len(data[parent].label) == 0:
            child_0 = data[parent].children[0]
            child_1 = data[parent].children[1]
            # Create new label for parent node
            data[parent].label = "(" + data[child_0].label + ":" + str(data[parent].age - data[child_0].age)[1:-1] + "," + data[child_1].label + ":" + str(data[parent].age - data[child_1].age)[1:-1] + ")"

    return data[list[-1]].label


#END CREATE_NEWICK
