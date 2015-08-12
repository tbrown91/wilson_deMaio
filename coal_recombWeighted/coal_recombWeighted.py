from coal_ancestry import combine_ancestralMaterial
from recomb_ancestry import split_ancestry
from DOT_fileWriting import DOT_coal
from DOT_fileWriting import DOT_recomb
from create_clonalFrame import create_initialList
from create_clonalFrame import remove_unancestral
from create_clonalFrame import remove_recomb
from create_clonalFrame import create_newick
from tree_comparison import compare_trees
from simulate_genome import simulate_nuc
from calculate_ancestry import calculate_totAncestry

import numpy as np
import random
import sys
import math
from copy import deepcopy
import matplotlib.pyplot as plt

#########################################
#           MAIN FUNCTION
#########################################

def main():

    num_leaves = 100
    mutation_rate = 1.0
    genome_length = 1000
    recomb_rate = 10.0
    dot_check = 0
    compare_check = 0
    tree_check = 0

    for i in range(len(sys.argv)):
        if sys.argv[i] == "-N":
            num_leaves = int(sys.argv[i+1])
        elif sys.argv[i] == "-M":
            mutation_rate = float(sys.argv[i+1])
        elif sys.argv[i] == "-G":
            genome_length = int(sys.argv[i+1])
        elif sys.argv[i] == "-R":
            recomb_rate = float(sys.argv[i+1])
        elif sys.argv[i] == "-D":
            dot_fileName = (sys.argv[i+1])
            dot_check = 1
        elif sys.argv[i] == "-n":
            nucleotide = int(sys.argv[i+1])
        elif sys.argv[i] == "-c":
            compare_check = 1
        elif sys.argv[i] == "-o":
            tree_check = 1
        elif (sys.argv[i] == "-h") or (sys.argv[i] == "--help"):
            print "Arguments to coal_recombAncestral.py as follows:"
            print "-N\tNumber of leaves (default 100)"
            print "-M\tMutation rate (default 1.0)"
            print "-G\tGenome length (default 1000)"
            print "-R\tRecombination rate (default 10.0)"
            print "-D\tWrite DOT file <file name> (default no file written)"
            print "-n\tNucleotide of interest for DOT file (default middle nucleotide)"
            print "-c\tCompare marginal trees to find detected recombination points (default FALSE)"
            print "-o\tCalculate marginal trees and simulate genomes (default FALSE)"
            sys.exit()


    print ("Number of leaves: " + str(num_leaves))
    print ("Mutation rate: " + str(mutation_rate))
    print ("Genome length: " + str(genome_length))
    print ("Recombination rate: " + str(recomb_rate))


    if dot_check == 1:
        print ("Writing DOT file...")
        if 'nucleotide' in locals():
            print ("Nucleotide of interest: " + str(nucleotide))
        else:
            nucleotide = int(math.floor((genome_length-1) / 2))
            print ("Nucleotide of interest: " + str(nucleotide) + "(default)")

        #Write DOT file header
        dot_file  = open(dot_fileName,"w")
        dot_file.write('digraph asde91 {fontsize=5;ranksep="0.02";ratio=fill;size="10,10";\n')
        dot_file.write('edge[arrowhead=none];\n')
        dot_file.write('{rank=same;')
        for i in range(num_leaves):
            dot_file.write('%d[shape=point] ' % i)
        dot_file.write('}\n')

        graph_nodes = []
        graph_edges = []
        for i in range(num_leaves):
            graph_nodes.append(str(i) + '[shape=point,width=0.00,height=0.00]\n')

    print graph_nodes[0]

    #File to store recombination points
    recomb_file = open("recombination_points.log","w")

    #Define Node class
    class Node:
        def __init__(self):
            self.age = 0.0
            self.children = []
            self.parent = []
            self.ancestral_material = []
            self.tot_ancestralMaterial = 0
            self.ancestral = 0
            self.label = ''
            self.sequence = ''

    #Define initial leaf nodes
    node_list = []
    for i in range(num_leaves):
        node_list.append(Node())
        node_list[i].ancestral_material = [[0,genome_length-1]]
        node_list[i].ancestral = 1
        node_list[i].label = str(i)
        node_list[i].number = i
        node_list[i].tot_ancestralMaterial = genome_length

    live_nodeList = []
    for i in range(num_leaves):
        live_nodeList.append(i)

    remaining_nodes = num_leaves
    new_node = num_leaves
    current_time = 0.0

    #########################################
    # Simulate coalescent-recombination graph
    #########################################

    while remaining_nodes > 1:
        current_ancestralMaterial = 0
        for i in live_nodeList:
            current_ancestralMaterial = current_ancestralMaterial + node_list[i].tot_ancestralMaterial

        #Simulate time to next event
        current_time = current_time - 2 * np.log(1 - np.random.rand(1)) / (remaining_nodes * (remaining_nodes - 1) + recomb_rate * (float(current_ancestralMaterial / genome_length)))

        #Determine whether event is coalescent or recombination
        reac_rand = np.random.rand(1)
        if reac_rand[0] < (float(remaining_nodes * (remaining_nodes - 1)) / (remaining_nodes * (remaining_nodes - 1)  + recomb_rate * (float(current_ancestralMaterial / genome_length)))):
            #Coalescent event occurs

            #Choose two children nodes to coalesce
            children = random.sample(live_nodeList,2)
            #Set the parent of the two children
            node_list[children[0]].parent.append(new_node)
            node_list[children[1]].parent.append(new_node)

            #Make new parent node
            node_list.append(Node())
            node_list[new_node].children = children
            node_list[new_node].age = current_time
            node_list[new_node].number = new_node

            # Create the combined parental ancestral material
            # Combine the two ancestries into one list
            # Call function to combine two parts of ancestral material
            new_ancestry = combine_ancestralMaterial(node_list[children[0]].ancestral_material,node_list[children[1]].ancestral_material)
            node_list[new_node].ancestral_material = new_ancestry[:]

            node_list[new_node].tot_ancestral_material = calculate_totAncestry(node_list[new_node].ancestral_material)

            if dot_check == 1:
                # Check if new parent is ancestral for nucleotide of interest
                ancestral_check = 0
                for i in range(len(node_list[new_node].ancestral_material)):
                    if (nucleotide >= node_list[new_node].ancestral_material[i][0]) and (nucleotide <= node_list[new_node].ancestral_material[i][1]):
                        ancestral_check = 1
                        break
                if ancestral_check == 1:
                    node_list[new_node].ancestral = 1

                [graph_nodes,graph_edges] = DOT_coal(node_list[new_node].ancestral,new_node,children,node_list[children[0]].ancestral,node_list[children[1]].ancestral,graph_nodes,graph_edges)

            temp_liveNodes = []
            temp_liveNodes.append(new_node)
            for i in live_nodeList:
                if (i != children[0]) and (i != children[1]):
                    temp_liveNodes.append(i)
            live_nodeList = temp_liveNodes[:]

            remaining_nodes = remaining_nodes-1
            new_node = new_node + 1

        else:
            #Recombination event ocurrs

            #Choose a child weighted by amount of remaining ancestral material
            ancestral_array = np.zeros(len(live_nodeList))
            for i in range(len(live_nodeList)):
                ancestral_array[i] = float(node_list[live_nodeList[i]].tot_ancestralMaterial) / current_ancestralMaterial
            cumsum_ancestral = np.cumsum(ancestral_array)
            recomb_rand = np.random.rand(1)
            for i in range(len(live_nodeList)):
                if recomb_rand <= cumsum_ancestral[i]:
                    child = live_nodeList[i]

            #Choose which point of the ancestral material is recombinant
            break_point = random.sample(range(node_list[child].ancestral_material[0][0],node_list[child].ancestral_material[len(node_list[child].ancestral_material)-1][1]),1)
            break_point = break_point[0]

            #Write accepted recombination point to file
            recomb_file.write("%d\n" % break_point)

            #Set the two parents
            node_list[child].parent = [new_node,new_node+1]

            #Make new parent nodes
            node_list.append(Node())
            node_list.append(Node())
            node_list[new_node].children.append(child)
            node_list[new_node+1].children.append(child)
            node_list[new_node].age = current_time
            node_list[new_node+1].age = current_time
            node_list[new_node].number = new_node
            node_list[new_node+1].number = new_node+1

            # Split ancestral material into two parents
            # Call function to split ancestral material
            [node_list[new_node].ancestral_material,node_list[new_node+1].ancestral_material] = split_ancestry(node_list[child].ancestral_material,break_point)

            node_list[new_node].tot_ancestral_material = calculate_totAncestry(node_list[new_node].ancestral_material)
            node_list[new_node+1].tot_ancestral_material = calculate_totAncestry(node_list[new_node+1].ancestral_material)


            if dot_check == 1:
                #Check if two parents are ancestral to nucleotide of interest
                ancestral_check = 0
                for i in range(len(node_list[new_node].ancestral_material)):
                    if (nucleotide >= node_list[new_node].ancestral_material[i][0]) and (nucleotide <= node_list[new_node].ancestral_material[i][1]):
                        ancestral_check = 1
                        break
                if ancestral_check == 1:
                    node_list[new_node].ancestral = 1
                ancestral_check = 0
                for i in range(len(node_list[new_node+1].ancestral_material)):
                    if (nucleotide >= node_list[new_node+1].ancestral_material[i][0]) and (nucleotide <= node_list[new_node+1].ancestral_material[i][1]):
                        ancestral_check = 1
                        break
                if ancestral_check == 1:
                    node_list[new_node+1].ancestral = 1

                #Call function to write recombination nodes and edges to file
                [graph_nodes,graph_edges] = DOT_recomb(new_node,new_node+1,child,node_list[child].ancestral,node_list[new_node].ancestral,node_list[new_node+1].ancestral,graph_nodes,graph_edges)

            if (len(node_list[new_node].ancestral_material) == 0) or (len(node_list[new_node+1].ancestral_material) == 0):
                print ("Child: "  +str(node_list[child].ancestral_material))
                print ("Break point: " + str(break_point))
                print ("Recombination parent 1: " + str(node_list[new_node].ancestral_material))
                print ("Recombination parent 2: " + str(node_list[new_node+1].ancestral_material))

            temp_liveNodes = []
            temp_liveNodes.append(new_node)
            temp_liveNodes.append(new_node+1)
            for i in live_nodeList:
                if (i != child):
                    temp_liveNodes.append(i)
            live_nodeList = temp_liveNodes[:]

            remaining_nodes = remaining_nodes + 1
            new_node = new_node + 2
        #END ELSE
    #END WHILE

    recomb_file.close()

    if dot_check == 1:
        for i in range(len(graph_nodes)):
            dot_file.write(str(graph_nodes[i]))
        for i in range(len(graph_edges)):
            dot_file.write(str(graph_edges[i]))
        dot_file.write('}')
        dot_file.close()

    #Check for any ancestral errors
    for i in range(genome_length):
        check = 0
        for j in range(len(node_list[len(node_list)-1].ancestral_material)):
            if (i >= node_list[len(node_list)-1].ancestral_material[j][0]) and (i <= node_list[len(node_list)-1].ancestral_material[j][1]):
                check = 1
        if check == 0:
            print ("Error, nucleotide: " + str(i))
            print ("Root ancestry: " + str(node_list[len(node_list)-1].ancestral_material))
            sys.exit()


    #########################################
    # Create clonal trees for each nucleotide
    #########################################

    if tree_check == 1:
        clonal_trees = []
        #Create tree for each nucleotide
        tree_file = open('trees.nwk','w')
        #for nuc in range(genome_length):
        print "Calculating marginal trees and simulating sequences..."
        for nuc in range(genome_length):
            #Copy node list
            temp_list = []
            temp_list = deepcopy(node_list)
            #Create a list of all clonal nodes
            clonal_list = []
            clonal_list = create_initialList(temp_list,nuc)

            #Remove un-ancestral parents and children
            clonal_list = remove_unancestral(temp_list,clonal_list)

            #Remove recombinant nodes
            [temp_list,clonal_list] = remove_recomb(temp_list,clonal_list)

            #Remove last node if recombination causes elongation of the tree
            if len(temp_list[clonal_list[-1]].children) == 1:
                clonal_list = clonal_list[0:-1]

            #Set last node as root
            temp_list[clonal_list[-1]].parent = []

            # Calculate Newick tree for nucleotide and write to file
            newick_tree = []
            newick_tree = create_newick(temp_list,clonal_list,num_leaves)
            newick_tree = newick_tree + ";"

            tree_file.write("[%d]%s\n" % (nuc,newick_tree))
            clonal_trees.append(newick_tree)

            #Simulate nucleotide down clonal frame
            node_list = simulate_nuc(nuc,node_list,temp_list,clonal_list,mutation_rate,num_leaves)

        sequence_file = open("genomes.fasta","w")
        for i in range(num_leaves):
            sequence_file.write(">%d\n" % i)
            sequence_file.write("%s\n" % node_list[i].sequence)
        sequence_file.close()
        if compare_check == 1:
            print "Comparing marginal trees..."
            tree_comparisons = np.zeros((genome_length,genome_length))
            for i in range(genome_length):
                for j in range(genome_length):
                    if (j == i+1) or (j == i-1):
                        tree_comparisons[i][j] = compare_trees(clonal_trees[i],clonal_trees[j])
                    elif (i == j):
                        tree_comparisons[i][j] = 0.0
                    else:
                        tree_comparisons[i][j] = float('nan')
            #Plot heat map of distance metric
            img = plt.pcolor(tree_comparisons)
            plt.colorbar(img)
            plt.show(img)

        tree_file.close()



if __name__ == "__main__":
    main()
