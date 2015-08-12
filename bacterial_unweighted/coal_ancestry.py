def combine_ancestralMaterial(child_ancestry1,child_ancestry2):

    child_1index = 0
    child_2index = 0
    current_start = 0
    current_end = 0
    new_ancestry = []
    while 1:
        #If one ancestry is finished, append the remaining intervals from the other
        #ancestry onto the new ancestry
        if child_1index == len(child_ancestry1):
            for i in range(child_2index,len(child_ancestry2)):
                new_ancestry.append([child_ancestry2[i][0],child_ancestry2[i][1]])
            child_2index = len(child_ancestry2)
            break
        elif child_2index == len(child_ancestry2):
            for i in range(child_1index,len(child_ancestry1)):
                new_ancestry.append([child_ancestry1[i][0],child_ancestry1[i][1]])
            child_1index = len(child_ancestry1)
            break
        else:
            #If the first two intervals share the same start site, find the
            #longer interval and remove the shorter from the "live list"
            if child_ancestry1[child_1index][0] <= child_ancestry2[child_2index][0]:
                if child_ancestry1[child_1index][1] >= child_ancestry2[child_2index][1]:
                    current_start = child_ancestry1[child_1index][0]
                    current_end = child_ancestry1[child_1index][1]
                    child_1index = child_1index + 1
                    child_2index = child_2index + 1
                else:
                    current_start = child_ancestry1[child_1index][0]
                    current_end = child_ancestry1[child_1index][1]
                    child_1index = child_1index + 1
            else:
                if child_ancestry2[child_2index][1] >= child_ancestry1[child_1index][1]:
                    current_start = child_ancestry2[child_2index][0]
                    current_end = child_ancestry2[child_2index][1]
                    child_1index = child_1index + 1
                    child_2index = child_2index + 1
                else:
                    current_start = child_ancestry2[child_2index][0]
                    current_end = child_ancestry2[child_2index][1]
                    child_2index = child_2index + 1
            while 1:
                if (child_1index < len(child_ancestry1)) and (child_2index < len(child_ancestry2)):
                    if current_end >= child_ancestry2[child_2index][1]:
                        child_2index = child_2index + 1
                    elif current_end >= child_ancestry1[child_1index][1]:
                        child_1index = child_1index + 1
                    elif current_end >= child_ancestry2[child_2index][0] - 1:
                        current_end = child_ancestry2[child_2index][1]
                        child_2index = child_2index + 1
                    elif current_end >= child_ancestry1[child_1index][0] - 1:
                        current_end = child_ancestry1[child_1index][1]
                        child_1index = child_1index + 1
                    else:
                        new_ancestry.append([current_start,current_end])
                        break
                elif (child_1index < len(child_ancestry1)):
                    if current_end >= child_ancestry1[child_1index][1]:
                        child_1index = child_1index + 1
                    elif current_end >= child_ancestry1[child_1index][0] - 1:
                        current_end = child_ancestry1[child_1index][1]
                        child_1index = child_1index + 1
                    else:
                        new_ancestry.append([current_start,current_end])
                        break
                elif (child_2index < len(child_ancestry2)):
                    if current_end >= child_ancestry2[child_2index][1]:
                        child_2index = child_2index + 1
                    elif current_end >= child_ancestry2[child_2index][0] - 1:
                        current_end = child_ancestry2[child_2index][1]
                        child_2index = child_2index + 1
                    else:
                        new_ancestry.append([current_start,current_end])
                        break
                else:
                    new_ancestry.append([current_start,current_end])
                    break
            #END WHILE
        #END ELSE
    #END WHILE
    return new_ancestry
#END COMBINE_ANCESTRALMATERIAL
