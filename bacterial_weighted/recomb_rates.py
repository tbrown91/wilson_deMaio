def clonal_recomb(ancestry_intervals,interval_starts,recomb_rate,genome_length,interval_rate):
    #For each node remaining in the tree calculate the rate of recombination
    #such that the recombinant interval intersects the ancestral interval

    #Calculate the rate of recombination beginning between ancestral intervals
    #and containing the first element of the ancestral interval
    ancestral_recombRate = 0
    if len(interval_starts) == 1:
        #Only one interval in the ancestral material

        if len(ancestry_intervals) == 2:
            #Interval wraps around the 'end' of the genome
            ancestral_recombRate = ancestral_recombRate + ((recomb_rate * (1 - (1-interval_rate)**(ancestry_intervals[1][0] - ancestry_intervals[0][1]))) / (2 * genome_length * interval_rate))

        else:
            #Interval is contained within the start and end of the genome
            ancestral_recombRate = ancestral_recombRate + ((recomb_rate * (1 - (1-interval_rate)**(genome_length + ancestry_intervals[0][0] - ancestry_intervals[0][1]))) / (2 * genome_length * interval_rate))

    elif len(interval_starts) > 1:
        #More than one interval in the ancestral material

        if (ancestry_intervals[0][0] == 0) and (ancestry_intervals[-1][1] == genome_length - 1):
            #Last interval wraps around the end of the genome
            for i in range(1,len(ancestry_intervals)):
                ancestral_recombRate = ancestral_recombRate + ((recomb_rate * (1 - (1-interval_rate)**(ancestry_intervals[i][0] - ancestry_intervals[i-1][1]))) / (2 * genome_length * interval_rate))

        else:
            #All intervals are contained within the genome - gap between last and
            #first interval wraps around the end of the genome
            ancestral_recombRate = ancestral_recombRate + ((recomb_rate * (1 - (1-interval_rate)**(genome_length + ancestry_intervals[0][0] - ancestry_intervals[-1][1]))) / (2 * genome_length * interval_rate))

            for i in range(1,len(ancestry_intervals)):
                ancestral_recombRate = ancestral_recombRate + ((recomb_rate * (1 - (1-interval_rate)**(ancestry_intervals[i][0] - ancestry_intervals[i-1][1]))) / (2 * genome_length * interval_rate))


    #Calculate the rate of a recombinant break starting within an ancestral
    #interval and splitting the ancestral material
    if len(interval_starts) == 0:
        #Entire genome is in the ancestral material
        ancestral_recombRate = ancestral_recombRate + recomb_rate * (ancestry_intervals[0][1] - ancestry_intervals[0][0]) / (2 * genome_length)

    else:
        if (ancestry_intervals[0][0] == 0) and (ancestry_intervals[-1][1] == genome_length - 1):
            #Last interval wraps around the end of the genome
            ancestral_recombRate = ancestral_recombRate + (recomb_rate * (ancestry_intervals[0][1] - ancestry_intervals[0][0] + 1)) / (2 * genome_length)
            for i in range(1,len(ancestry_intervals)):
                ancestral_recombRate = ancestral_recombRate + (recomb_rate * (ancestry_intervals[i][1] - ancestry_intervals[i][0])) / (2 * genome_length)

        else:
            #All intervals are contained within the genome
            for i in range(len(ancestry_intervals)):
                ancestral_recombRate = ancestral_recombRate + (recomb_rate * (ancestry_intervals[i][1] - ancestry_intervals[i][0])) / (2 * genome_length)

    return ancestral_recombRate
#END CLONAL_RECOMB

def non_clonalRecomb(ancestry_intervals,interval_starts,recomb_rate,genome_length,interval_rate):
    #For each node remaining in the tree calculate the rate of recombination such
    #that any recombination event takes the entire ancestral material
    non_ancestralRate = 0

    #Calculate the rate at which recombination begins and ends in the same
    #region between ancestral interals
    if len(interval_starts) == 1:
        #Only one interval within the genome

        if len(ancestry_intervals) == 2:
            #Interval wraps around the end of the genome
            non_ancestralRate = non_ancestralRate + (recomb_rate * ((1 - interval_rate)**genome_length) * ((1 - interval_rate)**(-(ancestry_intervals[1][0] - ancestry_intervals[0][1])))) / (2 * genome_length * interval_rate)

        else:
            #Interval is contained within the genome
            non_ancestralRate = non_ancestralRate + (recomb_rate * ((1 - interval_rate)**genome_length) * ((1 - interval_rate)**(-(genome_length + ancestry_intervals[0][0] - ancestry_intervals[0][1])))) / (2 * genome_length * interval_rate)

    elif len(interval_starts) > 1:
        #More than one ancestral interval within the genome

        if (ancestry_intervals[0][0] == 0) and (ancestry_intervals[-1][1] == genome_length - 1):
            #Last interval wraps around the end of the genome
            for i in range(1,len(ancestry_intervals)):
                non_ancestralRate = non_ancestralRate + (recomb_rate * ((1 - interval_rate)**genome_length) * ((1 - interval_rate)**(-(ancestry_intervals[i][0] - ancestry_intervals[i-1][1])))) / (2 * genome_length * interval_rate)

        else:
            #All intervals contained within the genome
            non_ancestralRate = non_ancestralRate + (recomb_rate * ((1 - interval_rate)**genome_length) * ((1 - interval_rate)**(-(genome_length + ancestry_intervals[0][0] - ancestry_intervals[-1][1])))) / (2 * genome_length * interval_rate)
            for i in range(1,len(ancestry_intervals)):
                non_ancestralRate = non_ancestralRate + (recomb_rate * ((1 - interval_rate)**genome_length) * ((1 - interval_rate)**(-(ancestry_intervals[i][0] - ancestry_intervals[i-1][1])))) / (2 * genome_length * interval_rate)


    #Calculate the rate of recombination which begins within an ancestral interval
    #and the recombinant interval takes the entire genome
    if len(interval_starts) == 0:
        #Ancestral material contains the entire genome
        non_ancestralRate = non_ancestralRate + recomb_rate * ((1 - interval_rate)**(genome_length-1)) * (ancestry_intervals[0][1] - ancestry_intervals[0][0]) / (2 * genome_length)

    else:

        if (ancestry_intervals[0][0] == 0) and (ancestry_intervals[-1][1] == genome_length - 1):
            #Last interval wraps around the end of the genome
            non_ancestralRate = non_ancestralRate + recomb_rate * ((1 - interval_rate)**(genome_length-1)) * (ancestry_intervals[0][1] - ancestry_intervals[0][0] + 1) / (2 * genome_length)
            for i in range(1,len(ancestry_intervals)):
                non_ancestralRate = non_ancestralRate + recomb_rate * ((1 - interval_rate)**(genome_length-1)) * (ancestry_intervals[i][1] - ancestry_intervals[i][0]) / (2 * genome_length)

        else:
            #All intervals are contained within the genome
            for i in range(len(ancestry_intervals)):
                non_ancestralRate = non_ancestralRate + recomb_rate * ((1 - interval_rate)**(genome_length-1)) * (ancestry_intervals[i][1] - ancestry_intervals[i][0]) / (2 * genome_length)

    return non_ancestralRate
#END NON_CLONALRECOMB
