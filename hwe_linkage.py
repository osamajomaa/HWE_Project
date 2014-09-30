import csv
import sys
import scipy
from sets import Set
from collections import OrderedDict
from scipy.stats.stats import chisquare

def parse_genotype_file(gt_file):
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
    ppl = 0
    n = len(matrix)
    for i in range(0,n):
        if ((matrix[i][locus] == a1 and matrix[i][locus+1] == a2) or (matrix[i][locus] == a2 and matrix[i][locus+1] == a1)):
                ppl += 1
    return ppl

    
def pop_fraction(N, n, m, total_pop): 
    P = [[[0 for k in xrange(n)] for j in xrange(n)] for i in xrange(m/2)]
    for k in range(0, m/2):
        for i in range (0,n):
            for j in range (0,n):
                if (i <= j):              
                    P[k][i][j] = float(N[k][i][j])/total_pop
    return P
    

def chrom_fraction(P, n, m):
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
        print p
        if (p < sigValue):
            nhwe_loci.append(loci[k*2])     
    
    return nhwe_loci

def ld_ppl_allel(matrix, a1, a2, a3, a4, l1, l2):
    ppl = 0
    n = len(matrix)
    for i in range(0,n):
        if (((matrix[i][l1] == a1 and matrix[i][l1+1] == a2) or (matrix[i][l1] == a2 and matrix[i][l1+1] == a1))
            and (((matrix[i][l2] == a3 and matrix[i][l2+1] == a4) or (matrix[i][l2] == a4 and matrix[i][l2+1] == a3)))):
                ppl += 1
    return ppl


def ld_actual_value(P, n, m, matrix, pop_size):
    
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
            #print p
            if (p < sigValue):
                ld_loci.append((loci[k*2], loci[s*2]))     
    
    return ld_loci

def calc_lindis(pop_size, loci, allels, matrix, P, L, sigValue):
    n = len(allels)
    m = len(loci)
    N, E = ld_actual_value(P, n, m, matrix, pop_size)
    return ld_chi_square(E, n, m, N, L, sigValue, loci, pop_size)
    


if __name__ == "__main__":
    genotype = "genotypes/genotype.african_american.csv"
    sig_value = 0.5
    #gt_file = sys.argv[1]
    #sig_value = sys.argv[2]
    pop_size, loci, allels, matrix = parse_genotype_file(genotype)
    nhwe_loci, P, L = calc_hwe(pop_size, loci, allels, matrix, sig_value)
    ld_loci = calc_lindis(pop_size, loci, allels, matrix, P, L, sig_value)
    for locus in ld_loci:
        print locus
    
    
                        

        