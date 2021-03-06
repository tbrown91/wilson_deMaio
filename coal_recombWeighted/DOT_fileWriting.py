def DOT_coal(parent_ancestral,parent_number,children_numbers,child_1ancestral,child_2ancestral,nodes,edges):
    #Write the nodes and edges to DOT file following a coalescent event
    #If children and parent ancestries contain the nucleotide of interest, draw bold, otherwise gray
    if  parent_ancestral == 0:
        #Parent is not ancestral, therefore both edges are gray
        nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[color=gray]\n')
        edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[color=gray]\n')
    else:
        #Check if first child is ancestral
        if child_1ancestral == 1:
            nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00]\n')
            edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[style=bold]\n')
        else:
            nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[color=gray]\n')
        #Check if first child is ancestral
        if child_2ancestral == 1:
            nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00]\n')
            edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[style=bold]\n')
        else:
            nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[color=gray]\n')

    return [nodes,edges]
#END DOT_COAL()


def DOT_recomb(parent_1Num,parent_2Num,child_num,child_ancestral,parent_1Ancestral,parent_2Ancestral,nodes,edges):
    #Write recombinant nodes and edges
    if child_ancestral == 0:
        nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[color=gray]\n')
        nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[color=gray]\n')
    else:
        #Check if new parents are ancestral and draw appropriately
        if parent_1Ancestral == 1:
            nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00]\n')
            edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[style=bold]\n')
        else:
            nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[color=gray]\n')
        if parent_2Ancestral == 1:
            nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00]\n')
            edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[style=bold]\n')
        else:
            nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[color=gray]\n')

    return [nodes,edges]
#END DOT_RECOMB
