import csv
import sys
import scipy
# from sets import set
from collections import OrderedDict
from scipy.stats.stats import chisquare

class HWE_Linkage:
    def __init__(self):
        self.loci = []
        self.alleles = []
        self.matrix = []
        self.pop_size = 0
        self.n = 0
        self.m = 0
        self.alinlocus = OrderedDict()

    def compute(self, gt_file, sig_value):
        self.parse_genotype_file(gt_file)
        self.calc_hwe(sig_value)


    def parse_genotype_file(self, gt_file):
        """
        Parse genotype file. Store as 2d self.matrix self.matrix where self.matrix[i][k] is the
        allele for person i at locus k.
        """
        rows = sum(1 for row in csv.reader( open(gt_file, 'rU'))) -1
        with open(gt_file, 'rU') as f:
            reader = csv.reader(f)
            rownum = 0
            unique_alleles = set()
            for row in reader:
                if rownum == 0:
                    cols = len(row)-1
                    self.matrix = [[0 for k in range(cols)] for j in range(rows)]
                    colnum = 0
                    for col in row:
                        if colnum > 0:
                            self.loci.append(col)
                        colnum += 1
                else:
                    colnum = 0
                    for col in row:
                        if colnum > 0:
                            self.matrix[rownum-1][colnum-1] = col
                            unique_alleles.add(col)
                        colnum += 1
                rownum += 1
        self.pop_size = rownum-1
        self.alleles = list(unique_alleles)
        self.n = len(self.alleles)
        self.m = len(self.loci)

    def calc_hwe(self, sig_value):
        """
        Determine the loci that are self.not in HWE, and store list of self.names. Compute
        3d self.matrix P where P[k][i][j] = percentage of population that has
        alleles i and j at locus k, and L where L[k] = self.number of different alleles
        found at locus k.
        """
        self.calculate_alleles_loci();
        N, P, E, p = self.calculate_allele_frequencies();

        for locus in self.find_not_hwe_loci(N, E, sig_value):
            print(locus)






    def calculate_alleles_loci(self):
        '''
        Calculates an ordered dictionary with loci self.names as keys, self.mapping to the set
        of alleles found at this locus
        '''
        for loc in self.loci:
            self.alinlocus[loc] = {}
        for k in range(0, self.m):
            namek = self.loci[k]
            if not namek in self.alinlocus[namek]:
                self.alinlocus[namek] = {}
            for i in range(0,self.n):
                ival = self.matrix[i][k]
                if not ival in self.alinlocus[namek]:
                    self.alinlocus[namek][ival] = 0
                self.alinlocus[namek][ival] += 1

    def calculate_allele_frequencies(self):
        #N = [[[0 for k in range(self.n)] for j in range(self.n)] for i in range(self.m//2)]
        N = {}
        #P = [[[0 for k in range(self.n)] for j in range(self.n)] for i in range(self.m//2)]
        P = {}
        k = 0
        for namek in self.alinlocus:
            if not namek in N:
                N[namek] = {}
                P[namek] = {}
            for i in range(0, self.n):
                ival = self.matrix[i][k]
                jval = self.matrix[i][k+1]
                if ival in self.alinlocus[namek] and jval in self.alinlocus[namek]:
                    mini = min(ival, jval)
                    maxi = max(ival, jval)
                    if not (mini, maxi) in N[namek]:
                        N[namek][(mini,maxi)] = 0
                    N[namek][(mini,maxi)] += 1
            for (i, j) in N[namek]:
                P[namek][(i, j)] = float(N[namek][(i,j)])/self.pop_size

            k += 2
        E, p = self.calculate_allele_chrom_fractions(P)
        return N, P, E, p


    def calculate_allele_chrom_fractions(self, P):
        '''
        Calculates and returns p, a 2d self.matrix where p[k][i] is the percentage of
        chromosomes that have allele i at locus k
        '''
        #p = [[0 for k in range(self.n)] for j in range(self.m//2)]
        #E = [[[0 for k in range(self.n)] for j in range(self.n)] for i in range(self.m//2)]
        p = {}
        E = {}
        for k in range(0, self.m//2):
            namek = self.loci[k*2]
            if not namek in p:
                p[namek] = {}
                E[namek] = {}
            for ival in self.alinlocus[namek]:
                p[namek][ival] = self.alinlocus[namek][ival]/(self.pop_size*2)
            for ival, jval in P[namek]:
                E[namek][(ival,jval)] = self.calculate_expected(p, ival, jval, k*2, namek)

        return E, p



    def calculate_expected(self, p, ival, jval, k, namek):
        '''
        Calculates and returns a 3d self.matrix E where E[k][i][j] is the percent of
        people expected to have alleles i, j at locus k.
        '''
        pki = float(p[namek][ival]) if ival in p[namek] else 1
        if (ival != jval):
            pkj = float(p[namek][jval]) if jval in p[namek] else 1
            return 2 * pki *pkj * self.pop_size
        else:
            return pki * pki * self.pop_size


    def find_not_hwe_loci(self, N, E, sigValue):
        '''
        Determines and returns list of loci that are self.not in HWE.
        '''
        not_hwe_loci = []
        for k in range(0, self.m//2):
            obs = []
            exp = []
            namek = self.loci[k*2]
            for ival, jval in N[namek]:
                if (E[namek][(ival,jval)] != 0):
                    #print P[k][i][j]
                    exp.append(E[namek][(ival,jval)])
                    obs.append(N[namek][(ival,jval)])

            Lk = len(self.alinlocus[namek])
            ddof = 0.5 * Lk *(Lk-1)
            chisq, p = chisquare(obs, exp, ddof)
            if p < sigValue:
                not_hwe_loci.append(namek)
        return not_hwe_loci






if __name__ == "__main__":
    if len(sys.argv) < 3:
        genotype = "genotypes/genotype.african_american.csv"
        sig_value = 0.5
    else:
        genotype = sys.argv[1]
        if len(sys.argv) > 2:
            sig_value = float(sys.argv[2])
        else:
            sig_value = 0.5

    HWE_Linkage().compute(genotype, sig_value)

        #pop_size, loci, alleles, self.matrix = parse_genotype_file(genotype)
        #self.nhwe_loci, P, L = calc_hwe(pop_size, loci, alleles, self.matrix, sig_value)
        #for locus in self.nhwe_loci:
        #    print(locus)
        # To Karro: this self.mIGHT work, but we ran out of self.memory.. :(
        #ld_loci = calc_lindis(pop_size, loci, alleles, self.matrix, P, L, sig_value)
