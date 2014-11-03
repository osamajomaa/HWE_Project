import random
import csv
import sys
from collections import OrderedDict
from collections import deque

'''
This class generates a simulated population under a specified allele distribution
'''
class SimPop:
    def __init__(self, infile, outfile, pop_size, loci_num = None):
        self.infile = infile
        self.pop_size = pop_size
        self.outfile = outfile
        self.allels_dist = OrderedDict()
        self.loci_num = loci_num

    '''
    A weighted random selection algorithm.
    @param weights: Weights (distribution) of allels for each locus
    '''
    def pick_allel(self, weights):
        accProb = random.random()
        for allel in range(0, len(weights)):
            accProb -= weights[allel]
            if (accProb <= 0):
                return allel

    '''
    Read the input file which contains the allels distribution to an ordered dictionary the key is the locus
    name and value is the list of allels distribution in that locus
    '''
    def parse_infile(self):
        with open(self.infile) as f:
            if self.loci_num:
                last_lines = deque(f, self.loci_num)
                for line in last_lines:
                    self.add_allele_dist(line)
            else:
                for line in f:
                    self.add_allele_dist(line)

    def add_allele_dist(self, line):
        dist = line.split(",")
        self.allels_dist[dist[0]] = [float(a.strip()) for a in dist[1:]]

    '''
    Create the simulated population file. First row contains the loci names where each name is written twice
    as each individual will have two allels. Then for each row, the first value is the individual id
    and the rest are the allels he has in each locus.
    '''
    def create_pop_genotype(self):
        title_row = ['SampleCode']
        for locus in self.allels_dist.keys():
            title_row.append(locus)
            title_row.append(locus)
        with open(self.outfile, 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerow(title_row)
            for p in range(0, self.pop_size):
                row = ['person'+str(p)]
                for locus in self.allels_dist:
                    allel1 = self.pick_allel(self.allels_dist[locus])
                    allel2 = self.pick_allel(self.allels_dist[locus])
                    row.append(str(allel1))
                    row.append(str(allel2))
                writer.writerow(row)



'''
Command line arguments:
    argv[1]: input file
    argv[2]: output file
    argv[3]: population size
    argv[4]: # loci to use (so use last # loci)
'''
if __name__ == "__main__":
    if (len(sys.argv) < 4):
        print("Please enter at least three arguments. Input file name, output file name and the population size")
        sys.exit(0)

    if (len(sys.argv) == 4):
        sp = SimPop(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        sp = SimPop(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))

    #Parse the input file first
    sp.parse_infile()
    #Generate the simulated population genotype file
    sp.create_pop_genotype()

    print("Done! Check the file", sys.argv[2])
