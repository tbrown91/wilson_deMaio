def combine_ancestralMaterial(child_ancestry1,child_ancestry2):
    # Combine the two children ancestries into one parent ancestry
    unordered_ancestry = []
    for i in range(len(child_ancestry1)):
        unordered_ancestry.append(child_ancestry1[i])
    for i in range(len(child_ancestry2)):
        unordered_ancestry.append(child_ancestry2[i])

    #Order the ancestry list and remove duplicate ancestries
    ordered_ancestry = []
    while len(unordered_ancestry) > 0:
        if len(unordered_ancestry) == 1:
            ordered_ancestry.append(unordered_ancestry[0])
            break
        j = 0
        k = len(unordered_ancestry) + 1
        for i in range(1,len(unordered_ancestry)):
            if unordered_ancestry[i][0] < unordered_ancestry[j][0]:
                j = i
                k = len(unordered_ancestry) + 1
            elif (unordered_ancestry[i][0] == unordered_ancestry[j][0]) and (unordered_ancestry[i][1] < unordered_ancestry[j][1]):
                j = i
                k = len(unordered_ancestry) + 1
            elif (unordered_ancestry[i][0] == unordered_ancestry[j][0]) and (unordered_ancestry[i][1] == unordered_ancestry[j][1]):
                k = i
        a = []
        ordered_ancestry.append(unordered_ancestry[j])
        for i in range(len(unordered_ancestry)):
            if (i != j) and (i != k):
                a.append(unordered_ancestry[i])
        unordered_ancestry = a[:]
    #END WHILE

    #Combine overlapping ancestral material
    overlapped_ancestry = []
    while len(ordered_ancestry)  > 0:
        if len(ordered_ancestry) == 1:
            overlapped_ancestry.append(ordered_ancestry[0])
            break
        #k = len(ordered_ancestry)  + 1
        k = 0
        a = []
        start_site = ordered_ancestry[0][0]
        end_site = ordered_ancestry[0][1]
        if (end_site >= (ordered_ancestry[1][0] - 1)) and (end_site < (ordered_ancestry[1][1])):
            end_site = ordered_ancestry[1][1]
            k = 1
            a.append([start_site,end_site])
        elif (end_site >= (ordered_ancestry[1][0] - 1)) and (end_site >= (ordered_ancestry[1][1])):
            k = 1
        if k == 0:
            overlapped_ancestry.append([start_site,end_site])
        for i in range(len(ordered_ancestry)):
            if (i != k):
                a.append(ordered_ancestry[i])
        ordered_ancestry = a[:]
    #END WHILE

    return overlapped_ancestry
#END COMBINE_ANCESTRALMATERIAL
