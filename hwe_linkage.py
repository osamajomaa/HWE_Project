import csv
import sys
import scipy
from sets import Set
from collections import OrderedDict
from scipy.stats.stats import chisquare

def parse_genotype_file(gt_file):
    '''
    Parses genotype file. Returns 2d matrix matrix where matrix[i][k] is the
    allele for person i at locus k.
    '''
    loci = []
    allels = Set()
    matrix = []
    rows = sum(1 for row in csv.reader( open(gt_file, 'rU'))) -1
    with open(gt_file, 'rU') as f:
        reader = csv.reader(f)
        rownum = 0
        for row in reader:
            if rownum == 0:
                cols = len(row)-1
                matrix = [[0 for k in xrange(cols)] for j in xrange(rows)]
                colnum = 0
                for col in row:
                    if colnum > 0:
                        loci.append(col)
                    colnum += 1
            else:
                colnum = 0
                for col in row:
                    if colnum > 0:
                        matrix[rownum-1][colnum-1] = col
                        allels.add(col)
                    colnum += 1
            rownum += 1
    return rownum-1, loci, list(allels), matrix


def allels_loci(matrix, loci):
    '''
    Returns an ordered dictionary with loci names as keys, mapping to the set
    of alleles found at this locus
    '''
    alinlocus = OrderedDict()
    for loc in loci:
        alinlocus[loc] = Set()
    n = len(matrix)
    m = len(matrix[0])
    for i in range(0,n):
        for j in range(0,m):
            alinlocus[loci[j]].add(matrix[i][j])
    return alinlocus


def ppl_allels(matrix, a1, a2, locus):
    '''
    Returns the number of individuals with allels a1, a2 at locus locus
    '''
    ppl = 0
    n = len(matrix)
    for i in range(0,n):
        if ((matrix[i][locus] == a1 and matrix[i][locus+1] == a2) or (matrix[i][locus] == a2 and matrix[i][locus+1] == a1)):
                ppl += 1
    return ppl


def pop_fraction(N, n, m, total_pop):
    '''
    Calculates and returns P, a 3d matrix where P[k][i][j] is the percentage
    of the population with alleles i, j at locus k
    '''
    P = [[[0 for k in xrange(n)] for j in xrange(n)] for i in xrange(m/2)]
    for k in range(0, m/2):
        for i in range (0,n):
            for j in range (0,n):
                if (i <= j):
                    P[k][i][j] = float(N[k][i][j])/total_pop
    return P


def chrom_fraction(P, n, m):
    '''
    Calculates and returns p, a 2d matrix where p[k][i] is the percentage of
    chromosomes that have allele i at locus k
    '''
    p = [[0 for k in xrange(n)] for j in xrange(m/2)]
    for k in range(0, m/2):
        for i in range(0, n):
            Pji = 0
            for j in range(0, i):
                Pji += P[k][j][i]
            Pij = 0
            for j in range(i+1, n):
                Pij += P[k][i][j]
            p[k][i] = P[k][i][i] + 0.5*Pji + 0.5*Pij
    return p


def expected(p, n, m, total_pop, matrix):
    '''
    Calculates and returns a 3d matrix E where E[k][i][j] is the percent of
    people expected to have alleles i, j at locus k.
    '''
    E = [[[0 for k in xrange(n)] for j in xrange(n)] for i in xrange(m/2)]
    for k in range(0, m)[::2]:
        for i in range (0,n):
            for j in range (0,n):
                if (i <= j):
                    if (matrix[i][k] != matrix[j][k+1]):
                        E[k/2][i][j] = 2 * float(p[k/2][i]) * float(p[k/2][j])
                    else:
                        E[k/2][i][j] = float(p[k/2][i]) * float(p[k/2][i])
    return E


def chi_square(E, n, m, P, L, sigValue, loci, pop_size):
    '''
    Determines and returns list of loci that are not in HWE.
    '''
    nhwe_loci = []
    for k in range(0, m/2):
        obs = []
        exp = []
        for i in range(0, n):
            for j in range(0, n):
                if (E[k][i][j] != 0):
                    #print P[k][i][j]
                    exp.append(E[k][i][j]*pop_size)
                    obs.append(P[k][i][j]*pop_size)

        ddof = 0.5 * len(L[loci[k*2]])*(len(L[loci[k*2]])-1)
        chisq, p = chisquare(obs, exp, ddof)
        if (p < sigValue):
            nhwe_loci.append(loci[k*2])

    return nhwe_loci

def ld_ppl_allel(matrix, a1, a2, a3, a4, l1, l2):
    '''
    Calculates the number of people who have alleles a1, a2 at locus l1 and
    alleles a3, a4 at locus l2
    '''
    ppl = 0
    n = len(matrix)
    for i in range(0,n):
        if (((matrix[i][l1] == a1 and matrix[i][l1+1] == a2) or (matrix[i][l1] == a2 and matrix[i][l1+1] == a1))
            and (((matrix[i][l2] == a3 and matrix[i][l2+1] == a4) or (matrix[i][l2] == a4 and matrix[i][l2+1] == a3)))):
                ppl += 1
    return ppl


def ld_actual_value(P, n, m, matrix, pop_size):
    '''
    Determines and returns N, 6d matrix where N[k][s][i][j][r][t] = number of
    people in population with alleles i, j at locus k AND alleles r,t at locus
    s, and E, 6d matrix where E[k][s][i][j][r][t] = number of people
    expected to have  alleles i, j at locus k AND alleles r,t at locus s
    '''
    N = [[[[[[0 for k in xrange(n)]for j in xrange(n)]for i in xrange(n)] for l in xrange(n)]for q in xrange(m/2)] for u in xrange(m/2)]
    E = [[[[[[0 for k in xrange(n)]for j in xrange(n)]for i in xrange(n)] for l in xrange(n)]for q in xrange(m/2)] for u in xrange(m/2)]
    for k in range (0,m/2):
        for s in range (0, m/2):
            for i in range(0, n):
                for j in range(0,n):
                    for r in range(0,n):
                        for t in range(0,n):
                            if (i <=j and r <= t):
                                N[k][s][i][j][r][t] = ld_ppl_allel(matrix, i, j, r, t, k, s)
                                E[k][s][i][j][r][t] = P[k][i][j] * P[s][r][t] * pop_size
    return N, E



def calc_hwe(total_pop, loci, allels, matrix, sigValue):
    '''
    Determine the loci that are not in HWE. Return list of names of loci not
    in HWE, 3d matrix P where P[k][i][j] = percentage of population that has
    alleles i and j at locus k, and L where L[k] = number of different alleles
    found at locus k.
    '''
    L = allels_loci(matrix, loci);

    n = len(allels)
    m = len(loci)
    N = [[[0 for k in xrange(n)] for j in xrange(n)] for i in xrange(m/2)]
    s = 0
    for k in range(0, m)[::2]:
        for i in range (0,n):
            for j in range (0,n):
                if (i <= j):
                    N[s][i][j] = ppl_allels(matrix, allels[i], allels[j], k)
        s += 1

    P = pop_fraction(N, n, m, total_pop)

    p = chrom_fraction(P, n, m)

    E = expected(p, n, m, total_pop, matrix)

    nhwe_loci = chi_square(E, n, m, P, L, sigValue, loci, total_pop)
    total_pop
    return nhwe_loci, P, L


def ld_chi_square(E, n, m, P, L, sigValue, loci, pop_size):
    '''
    Determine the pairs of loci in linkage disequilibrium
    '''
    ld_loci = []
    for k in range(0, m/2):
        for s in range(0,m/2):
            obs = []
            exp = []
            for i in range(0, n):
                for j in range(0, n):
                    for r in range(0,n):
                        for t in range(0,n):
                            if (E[k][s][i][j][r][t] != 0):
                                exp.append(E[k][s][i][j][r][t])
                                obs.append(P[k][s][i][j][r][t])
            ddofk = 0.5 * len(L[loci[k*2]])*(len(L[loci[k*2]])-1)
            ddofs = 0.5 * len(L[loci[s*2]])*(len(L[loci[s*2]])-1)
            ddof = (ddofk-1) * (ddofs-1)
            chisq, p = chisquare(obs, exp, ddof)
            if (p < sigValue):
                ld_loci.append((loci[k*2], loci[s*2]))

    return ld_loci


def calc_lindis(pop_size, loci, allels, matrix, P, L, sigValue):
    '''
    This function returns a list of pairs of loci in linkage disequilibrium.
    (Hypothetically)
    '''
    n = len(allels)
    m = len(loci)
    N, E = ld_actual_value(P, n, m, matrix, pop_size)
    return ld_chi_square(E, n, m, N, L, sigValue, loci, pop_size)



if __name__ == "__main__":
    if len(sys.argv) < 3:
        genotype = "genotypes/genotype.african_american.csv"
        sig_value = 0.5
    else:
        genotype = sys.argv[1]
        if len(sys.argv) > 2:
            sig_value = sys.argv[2]
        else:
            sig_value = 0.5
    pop_size, loci, allels, matrix = parse_genotype_file(genotype)
    nhwe_loci, P, L = calc_hwe(pop_size, loci, allels, matrix, sig_value)
    for locus in nwe_loci:
        print(locus)
    # To Karro: this MIGHT work, but we ran out of memory.. :(
    #ld_loci = calc_lindis(pop_size, loci, allels, matrix, P, L, sig_value)
