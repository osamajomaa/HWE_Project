import csv
import sys
from sets import Set
from collections import OrderedDict

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
            for j in range(0, i-1):
                Pji += P[k][j][i]
            Pij = 0
            for j in range(i+1, n):
                Pij += P[k][i][j]
            p[k][i] = P[k][i][i] + 0.5*Pji + Pij
    return p
            

def calc_hwe(gt_file, sig_value):
    
    total_pop, loci, allels, matrix = parse_genotype_file(gt_file)
    
    alinlocus = allels_loci(matrix, loci);
    
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
    
    
    


if __name__ == "__main__":
    genotype = "genotypes/genotype.african_american.csv"
    sig_value = 0.5
    #gt_file = sys.argv[1]
    #sig_value = sys.argv[2]
    calc_hwe(genotype, sig_value)
    
                        

        