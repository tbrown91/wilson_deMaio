import numpy as np
import random

def simulate_nuc(nucleotide,node_list,data,list,mut_rate,N):
    node_list[list[-1]].sequence = node_list[list[-1]].sequence + random.choice('ACGT')
    current_node = list[-1]
    leaves_sim = 0
    while leaves_sim < N:
        if len(data[current_node].children) == 2:
            child_1 = data[current_node].children[0]
            child_2 = data[current_node].children[1]
            if len(node_list[child_1].sequence) == nucleotide:
                #Simulate to first child
                current_time = 0.0
                sim_time = data[current_node].age - data[child_1].age
                new_nuc = node_list[current_node].sequence[nucleotide]
                while current_time < sim_time:
                    current_time = current_time - (2 * np.log(1-np.random.rand(1)) / mut_rate )
                    while new_nuc == node_list[current_node].sequence[nucleotide]:
                        new_nuc = random.choice('ACGT')
                node_list[child_1].sequence = node_list[child_1].sequence + new_nuc
                current_node = child_1
            elif len(node_list[child_2].sequence) == nucleotide:
                #Simulate to second child
                current_time = 0.0
                sim_time = data[current_node].age - data[child_2].age
                new_nuc = node_list[current_node].sequence[nucleotide]
                while current_time < sim_time:
                    current_time = current_time - (2 * np.log(1-np.random.rand(1)) / mut_rate )
                    while new_nuc == node_list[current_node].sequence[nucleotide]:
                        new_nuc = random.choice('ACGT')
                node_list[child_2].sequence = node_list[child_2].sequence + new_nuc
                current_node = child_2
            else:
                current_node = data[current_node].parent[0]
        else:
            leaves_sim = leaves_sim + 1
            current_node = data[current_node].parent[0]
    for i in range(len(node_list)):
        if len(node_list[i].sequence) == nucleotide:
            node_list[i].sequence = node_list[i].sequence + random.choice('ACGT')
    #END WHILE
    return node_list

#END SIMULATE_NUC
