import csv
import sys
import scipy
# from sets import set
from collections import OrderedDict
from scipy.stats.stats import chisquare
import scipy.stats.mstats as mst
from scipy import stats
import numpy as np

class HWE_Linkage:
    def __init__(self):
        self.loci = []
        self.alleles = []
        self.matrix = []
        self.pop_size = 0
        self.n = 0
        self.m = 0
        self.alinlocus = {}
        self.detailed = False
        self.locations = {}

    def compute(self, gt_file, sig_value):
        self.parse_genotype_file(gt_file)
        print("Loci NOT in HWE:")
        self.calc_hwe(sig_value)
        print("Loci NOT in LE:")
        self.calc_le(sig_value)


    def parse_genotype_file(self, gt_file):
        """
        Parse genotype file. Store as 2d self.matrix self.matrix where self.matrix[i][k] is the
        allele for person i at locus k.
        """
        rows = sum(1 for row in csv.reader( open(gt_file, 'rU'))) -1
        with open(gt_file, 'rU') as f:
            reader = csv.reader(f)
            rownum = 0
            for row in reader:
                if rownum == 0:
                    cols = len(row)-1
                    self.matrix = [[0 for k in range(cols)] for j in range(rows)]
                    colnum = 0
                    for col in row:
                        if colnum > 0 and colnum%2==0: # and colnum < 3:
                            self.loci.append(col)
                        colnum += 1
                else:
                    colnum = 0
                    for col in row:
                        if colnum > 0: # and colnum < 3:
                            self.matrix[rownum-1][colnum-1] = col
                        colnum += 1
                rownum += 1
        self.pop_size = rownum-1
        self.n = self.pop_size
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
            for i in range(0,self.n):
                ival = self.matrix[i][k*2]
                jval = self.matrix[i][k*2 + 1]
                if not ival in self.alinlocus[namek]:
                    self.alinlocus[namek][ival] = 0
                self.alinlocus[namek][ival] += 1
                if not jval in self.alinlocus[namek]:
                    self.alinlocus[namek][jval] = 0
                self.alinlocus[namek][jval] += 1
                
    def calculate_allele_frequencies(self):
        N = {}
        P = {}
        k = 0
        for k in range(0, self.m):
            namek = self.loci[k]
            if not namek in N:
                N[namek] = OrderedDict()#{}
                P[namek] = OrderedDict()#{}
                self.locations[namek] = {}
            for i in range(0, self.n):
                ival = self.matrix[i][k*2]
                jval = self.matrix[i][k*2+1]
                if not (ival, jval) in N[namek]:
                    N[namek][(ival,jval)] = 0
                    self.locations[namek][(ival, jval)] = []
                N[namek][(ival,jval)] += 1
                self.locations[namek][(ival, jval)].append(i)
            for pair in N[namek]:
                P[namek][pair] = float(N[namek][pair])/self.pop_size
            k += 2
        E, p = self.calculate_allele_chrom_fractions(P)
        return N, P, E, p


    def calculate_allele_chrom_fractions(self, P):
        '''
        Calculates and returns p, a 2d self.matrix where p[k][i] is the percentage of
        chromosomes that have allele i at locus k
        '''
        p = {}
        E = {}
        for k in range(0, self.m):
            namek = self.loci[k]
            if not namek in p:
                p[namek] = {}
                E[namek] = {}
            for ival in self.alinlocus[namek]:
                p[namek][ival] = self.alinlocus[namek][ival]/(self.pop_size*2)
            for pair in P[namek]:
                E[namek][pair] = self.calculate_expected(p, pair[0], pair[1], namek)
        return E, p



    def calculate_expected(self, p, ival, jval, namek):
        '''
        Calculates and returns a 3d self.matrix E where E[k][i][j] is the percent of
        people expected to have alleles i, j at locus k.
        '''
        pki = float(p[namek][ival])# if ival in p[namek] else 1
        if (ival != jval):
            pkj = float(p[namek][jval])# if jval in p[namek] else 1
            return 2 * pki *pkj * self.pop_size
        else:
            return pki * pki * self.pop_size


    def find_not_hwe_loci(self, N, E, sigValue):
        '''
        Determines and returns list of loci that are self.not in HWE.
        '''
        not_hwe_loci = []
        for k in range(0, self.m):
            obs = []
            exp = []
            namek = self.loci[k]
            #sumsquarediff = 0
            for ival, jval in N[namek]:
                if (E[namek][(ival,jval)] != 0):
                    exp.append(E[namek][(ival,jval)])
                    obs.append(N[namek][(ival,jval)])
                    #diff = N[namek][(ival,jval)] - E[namek][(ival,jval)]
                    #sumsquarediff += (diff * diff)/E[namek][(ival,jval)]
            Lk = len(self.alinlocus[namek])
            ddof = 0.5 * Lk *(Lk-1)
            #pval = 1 - stats.chi2.cdf(sumsquarediff, ddof)
            chisq, p = mst.chisquare(np.array(obs), np.array(exp), len(obs)-ddof)
            if self.detailed:
               print(namek,'\t',chisq, '\t',ddof)
            if p < sigValue:
                not_hwe_loci.append(namek)
        return not_hwe_loci

    def calc_le(self, sig_value):
        """
        """
        self.calculate_alleles_loci();
        N, P, E, p = self.calculate_allele_frequencies();

        for locusPair in self.find_not_link_loci(N,P, E, sig_value):
            print(locusPair)
    
    def num_occurrences(self, namek, pairk, names, pairs):
        locations_ijk = self.locations[namek][pairk]
        locations_rts = self.locations[names][pairs]
        return len(list(set(locations_ijk) & set(locations_rts)))

    def num_expected(self, P, namek, pairk, names, pairs):
        return P[namek][pairk] * P[names][pairs] * self.pop_size

    def find_not_link_loci(self, N, P, E, sigValue):
        '''
        Determines and returns a list of pairs of loci that are self.not in linkage equilibrium.
        '''
        not_link_loci = []
        k = 0
        s = 1
        ijk_count = 0
        rts_count = 0
        while k < self.m - 1:
            s = k + 1
            namek = self.loci[k]
            while s < self.m:
                names = self.loci[s]
                ijkrts_obs = []
                ijkrts_exp = []
                sumsquarediff = 0
                for pairk in N[namek]:
                    for pairs in N[names]:
                        N_ijkrts = self.num_occurrences(namek, pairk, names, pairs)
                        E_ijkrts = self.num_expected(P, namek, pairk, names, pairs)
                        if E_ijkrts != 0:
                            ijkrts_obs.append(N_ijkrts)
                            ijkrts_exp.append(E_ijkrts)
                            diff = N_ijkrts - E_ijkrts
                            sumsquarediff += (diff * diff)/E_ijkrts
                LK_ijk = len(self.alinlocus[namek])
                LK_rts = len(self.alinlocus[names])
                ddof_ijk =LK_ijk * (LK_ijk + 1) * .5
                ddof_rts =LK_rts * (LK_rts + 1) * .5
                ddof = (ddof_ijk - 1) * (ddof_rts - 1)
                tmp = 0
                pval = 1 - stats.chi2.cdf(sumsquarediff, ddof)
                
                chisq, p = chisquare(ijkrts_obs, ijkrts_exp, len(ijkrts_obs) - ddof)
                #chisq, p = chisquare(ijkrts_obs, ijkrts_exp, len(ijkrts_obs) - ddof)
                if self.detailed:
                    print((namek, names),'\t',chisq, '\t',ddof)
                if p < sigValue:
                    not_link_loci.append((namek, names))
                s = s + 1
            k = k + 1
        return not_link_loci




if __name__ == "__main__":
    if len(sys.argv) < 2:
        genotypes = ["african_american", "american", "asian_american", "caucasian_american", "hispanic_american"]
        sig_value = 0.05
        for genotype in genotypes:
            hwe = HWE_Linkage()
            print('\n' + genotype.upper() + ":")
            genotype_fname = "genotypes/genotype." + genotype + ".csv"
            hwe.compute(genotype_fname, sig_value)
    else:
        genotype = sys.argv[1]
        if len(sys.argv) > 2:
            sig_value = float(sys.argv[2])
        else:
            sig_value = 0.05
        hwe = HWE_Linkage()
        hwe.compute(genotype, sig_value)



