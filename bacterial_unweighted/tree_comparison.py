from ete2 import Tree

def compare_trees(tree_1,tree_2):
    t1 = Tree(tree_1)
    t2 = Tree(tree_2)
    rf = t1.robinson_foulds(t2)[0]
    return rf
#END COMPARE_TREES
