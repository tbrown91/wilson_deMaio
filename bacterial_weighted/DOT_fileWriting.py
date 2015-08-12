def DOT_coal(parent_ancestral,parent_number,children_numbers,child_1ancestral,child_2ancestral,graph_nodes,graph_edges):
    #Write the nodes and edges to DOT file following a coalescent event
    #If children and parent ancestries contain the nucleotide of interest, draw bold, otherwise gray
    if  parent_ancestral == 0:
        #Parent is not ancestral, therefore both edges are gray
        graph_nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[color=gray]\n')
        graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[color=gray]\n')
    else:
        #Check if first child is ancestral
        if child_1ancestral == 1:
            graph_nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00]\n')
            graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[style=bold]\n')
        else:
            graph_nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[0]) + '[color=gray]\n')
        #Check if second child is ancestral
        if child_2ancestral == 1:
            graph_nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00]\n')
            graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[style=bold]\n')
        else:
            graph_nodes.append(str(parent_number) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            graph_edges.append(str(parent_number) + ' -> ' + str(children_numbers[1]) + '[color=gray]\n')

    return [graph_nodes,graph_edges]
#END DOT_COAL()


def DOT_recomb(parent_1Num,parent_2Num,child_num,child_ancestral,parent_1Ancestral,parent_2Ancestral,graph_nodes,graph_edges):
    #Write recombinant nodes and edges
    if child_ancestral == 0:
        #If child is not ancestral, neither will both parents be
        graph_nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        graph_edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[color=gray]\n')
        graph_nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
        graph_edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[color=gray]\n')
    else:
        #Check which new parent is ancestral and draw appropriately
        if parent_1Ancestral == 1:
            graph_nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00]\n')
            graph_edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[style=bold]\n')
            graph_nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            graph_edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[color=gray]\n')
        else:
            graph_nodes.append(str(parent_1Num) + '[shape=point,width=0.00,height=0.00,color=gray]\n')
            graph_edges.append(str(parent_1Num) + ' -> ' + str(child_num) + '[color=gray]\n')
            graph_nodes.append(str(parent_2Num) + '[shape=point,width=0.00,height=0.00]\n')
            graph_edges.append(str(parent_2Num) + ' -> ' + str(child_num) + '[style=bold]\n')
            
    return [graph_nodes,graph_edges]
#END DOT_RECOMB
