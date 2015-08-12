def combine_ancestralMaterial(child_ancestry1,child_ancestry2):
    #Take the two child ancestries and combine them into one parental ancestry
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
            break
        elif child_2index == len(child_ancestry2):
            for i in range(child_1index,len(child_ancestry1)):
                new_ancestry.append([child_ancestry1[i][0],child_ancestry1[i][1]])
            break
        else:
            #If the first two intervals share the same start site, find the
            #longer interval and remove the shorter from the "live list"
            if child_ancestry1[child_1index][0] <= child_ancestry2[child_2index][0]:

                if child_ancestry1[child_1index][1] >= child_ancestry2[child_2index][1]:
                    #Child 1 interval contains child 2 interval
                    current_start = child_ancestry1[child_1index][0]
                    current_end = child_ancestry1[child_1index][1]
                    child_1index = child_1index + 1
                    child_2index = child_2index + 1

                else:
                    #Child 1 interval start and end before child 2 interval start and end
                    current_start = child_ancestry1[child_1index][0]
                    current_end = child_ancestry1[child_1index][1]
                    child_1index = child_1index + 1

            else:

                if child_ancestry2[child_2index][1] >= child_ancestry1[child_1index][1]:
                    #Child 2 interval contains child 1 interval
                    current_start = child_ancestry2[child_2index][0]
                    current_end = child_ancestry2[child_2index][1]
                    child_1index = child_1index + 1
                    child_2index = child_2index + 1

                else:
                    #Child 2 interval start and end before child2 interval start and end
                    current_start = child_ancestry2[child_2index][0]
                    current_end = child_ancestry2[child_2index][1]
                    child_2index = child_2index + 1

            while 1:
                #Find overlapping regions and update new ancestry

                if (child_1index < len(child_ancestry1)) and (child_2index < len(child_ancestry2)):
                    #While there are still intervals to check in both children

                    if current_end >= child_ancestry2[child_2index][1]:
                        #Interval overlaps child 2 interval entirely
                        child_2index = child_2index + 1

                    elif current_end >= child_ancestry1[child_1index][1]:
                        #Interval overlaps child 1 interval entirely
                        child_1index = child_1index + 1

                    elif current_end >= child_ancestry2[child_2index][0] - 1:
                        #Interval overlaps or joins child 2 start site
                        current_end = child_ancestry2[child_2index][1]
                        child_2index = child_2index + 1

                    elif current_end >= child_ancestry1[child_1index][0] - 1:
                        #Interval overlaps or joins child 1 start site
                        current_end = child_ancestry1[child_1index][1]
                        child_1index = child_1index + 1

                    else:
                        #Interval does not overlap any more, write to new ancestry array
                        new_ancestry.append([current_start,current_end])
                        break

                elif (child_1index < len(child_ancestry1)):
                    #While there are still intervals in child 1 and not child 2

                    if current_end >= child_ancestry1[child_1index][1]:
                        #Interval overlaps child 1 interval entirely
                        child_1index = child_1index + 1

                    elif current_end >= child_ancestry1[child_1index][0] - 1:
                        #Interval overlaps or joins child 1 start site
                        current_end = child_ancestry1[child_1index][1]
                        child_1index = child_1index + 1

                    else:
                        #Interval does not overlap any more, write to new ancestry array
                        new_ancestry.append([current_start,current_end])
                        break

                elif (child_2index < len(child_ancestry2)):
                    #While there are still intervals in child 2 and not child 1

                    if current_end >= child_ancestry2[child_2index][1]:
                        #Interval overlaps child 2 interval entirely
                        child_2index = child_2index + 1

                    elif current_end >= child_ancestry2[child_2index][0] - 1:
                        #Interval overlaps or joins child 2 start site
                        current_end = child_ancestry2[child_2index][1]
                        child_2index = child_2index + 1

                    else:
                        #Interval does not overlap any more, write to new ancestry array
                        new_ancestry.append([current_start,current_end])
                        break
                else:
                    #No more intervals to check, write to new ancestry array and exit
                    new_ancestry.append([current_start,current_end])
                    break
            #END WHILE
        #END ELSE
    #END WHILE

    return new_ancestry
#END COMBINE_ANCESTRALMATERIAL

def coal_ancestryStarts (parent_ancestry,genome_length):
    #Count the number of interval start sites
    parent_starts = []

    if (parent_ancestry[0][0] != 0) or (parent_ancestry[len(parent_ancestry)-1][1] != genome_length - 1):
        #Last interval wraps around end of genome, do not include first interval
        parent_starts.append(parent_ancestry[0][0])

    for i in range(1,len(parent_ancestry)):
        #All intervals contained within genome, include all intervals
        parent_starts.append(parent_ancestry[i][0])

    return parent_starts
#END COAL_ANCESTRYSTARTS
