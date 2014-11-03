import csv
import sys
import scipy
import random
# from sets import set
from collections import OrderedDict
#from scipy.stats.stats import chisquare
#import scipy.stats.mstats as mst
from scipy import stats
import numpy as np
from sim_pop import SimPop
import argparse

class Family:

	# dna would be a 2d matrix (row = person, col = loci halfs, (row, col) = allel name)
	# csv = all possible people
	# parent1 = random row
	# parent2 = random row
	def __init__(self, pop_size, num_sims, verbosity):
		self.infile = "sim_input_sample.csv"
		self.csv = "generated_pop.csv"
		self.num_loci = 13
		self.gen_pop = self.simulate_pop(pop_size)
		self.verbosity = verbosity
		#self.simulate_pop(pop_size)
		self.perform_simulations(num_sims)


	def perform_simulations(self, num_sims):
		i = 0
		random_match_sum = 0
		sibling_match_sum = 0
		cousin_match_sum = 0
		# perform the experiment numTrials times
		while i < num_sims:
			p1, p2 = self.getRandomPair()
			random_match_sum += self.num_match(p1, p2)
			s1, s2 = self.getSiblings(p1,p2)
			sibling_match_sum += self.num_match(s1,s2)
			c1, c2 = self.getCousinsFromSibs(s1,s2)
			cousin_match_sum += self.num_match(c1,c2)
			i += 1
		#total_possible = (num_sims * self.num_loci * 2)
		avg_random = random_match_sum/num_sims #total_possible
		avg_sibs = sibling_match_sum/num_sims #total_possible
		avg_cousins = cousin_match_sum/num_sims #total_possible
		num_alleles = 2*self.num_loci/100
		print("Expected number of alleles shared\n\t1.Between two individuals:\t{:f}\t({:.6}%)\n\t2.Between two cousins:\t\t{:f}\t({:.6}%)\n\t3.Between two siblings:\t\t{:f}\t({:.6}%)\n\t".format(avg_random, avg_random/num_alleles, avg_cousins, avg_cousins/num_alleles, avg_sibs, avg_sibs/num_alleles))



	def getRandomPerson(self):
		rand1 = random.randrange(1, len(self.gen_pop))
		return self.gen_pop[rand1]

	# Get two random rows (not the same because no asexual reproduction)
	def getRandomPair(self):
		rand1 = random.randrange(1, len(self.gen_pop))
		rand2 = random.randrange(1, len(self.gen_pop))
		while(rand1 == rand2):
			rand2 = random.randrange(1, len(self.csv))
		return self.gen_pop[rand1], self.gen_pop[rand2]

	# parent1 and parent2 are arrays.
	# They would be a row of data as in the csv files.
	# the return is an array randomly created from the parents
	def generateChild(self, parent1, parent2):
		child = []
		i = 0
		end = len(parent1)
		while i < end:
			rand = random.choice([0,1])
			if(rand == 0):
				child.append(parent1[i])
			elif(rand == 1):
				child.append(parent2[i])
			else:
				print("There was an error choosing an allele.")
			i = i + 1
		return child

	# the method reurns 2 children
	def getSiblings(self, parent1, parent2):
		return self.generateChild(parent1, parent2), self.generateChild(parent1, parent2)

	# Creates a child and its cusion
	# the two parameters are grand parents
	def getCousins(self, gParent1, gParent2):
		# makes a parent and aunt/uncle
		s = self.getSiblings(gParent1, gParent2)
		# first child's other parent
		r1 = self.getRandomPerson()
		# the first child
		c1 = self.generateChild(s[0], r1)
		# the cusion's other parent
		r2 = self.getRandomPerson()
		# create the cusion
		c2 = self.generateChild(s[1], r2)
		return (c1, c2)

	# Creates a child and its cusion
	# the two parameters are grand parents
	def getCousinsFromSibs(self, s1, s2):
		# first child's other parent
		r1 = self.getRandomPerson()
		# the first child
		c1 = self.generateChild(s1, r1)
		# the cusion's other parent
		r2 = self.getRandomPerson()
		# create the cusion
		c2 = self.generateChild(s2, r2)
		return c1, c2

	def num_match(self, person1, person2):
		i = 0
		numLoci = len(person1) #/ 2
		numMatches = 0
		while i < numLoci:
			if(person1[i] == person2[i] and person1[i + 1] == person2[i + 1]):
				numMatches = numMatches + 1
			i = i + 2
		return numMatches

	def percentMatch(self, person1, person2):
		numLoci = len(person1) / 2
		return (self.num_match(person1, person2) / numLoci)


	# significance will be used to determine how good a percent match is
	def compareManyPeople(self, numTrials, sig):
		i = 0;
		significantMatches = 0
		# perform the experiment numTrials times
		while i < numTrials:
			p1, p2 = self.getRandomPair()
			match = percentMatch(s1, s2)
			if match > sig:
				significantMatches = significantMatches + 1
			i = i + 1
		return (significantMatches / numTrials)

	# significance will be used to determine how good a percent match is
	def compareManySiblings(self, numTrials, sig):
		i = 0;
		significantMatches = 0
		# perform the experiment numTrials times
		while i < numTrials:
			s1 = generateChild(self.parent1, self.parent2)
			s2 = generateChild(self.parent1, self.parent2)
			match = percentMatch(s1, s2)
			if match > sig:
				significantMatches = significantMatches + 1
			i = i + 1
		return (significantMatches / numTrials)

	def compareManyCousins(self, numTrials, sig):
		i = 0;
		significantMatches = 0
		# perform the experiment numTrials times
		while i < numTrials:
			cusions = getCousins(self.parent1, self.parent2)
			match = percentMatch(cusions[0], cusions[1])
			if match > sig:
				significantMatches = significantMatches +1
			i = i + 1
		return (significantMatches / numTrials)

	def simulate_pop(self, pop_size):
		sp = SimPop(self.infile, self.csv, pop_size, self.num_loci)
		#Parse the input file first
		sp.parse_infile()
		#Generate the simulated population genotype file
		sp.create_pop_genotype()
		pop = []
		with open(self.csv) as f:
			first = True
			for line in f:
				if first:
					first = False
					continue
				person = line.split(",")
				pop.append([float(a.strip()) for a in person[1:]])
		return pop


if __name__ == "__main__":
	# population size and number realizations
	parser = argparse.ArgumentParser()
	parser.add_argument("-m","--pop_size", type=int, help="size of simulated population", default=100)
	parser.add_argument("-n","--num_runs", type=int, help="number of realizations", default=100)
	parser.add_argument("-v", "--verbose", action="store_true", default=False)
	args = parser.parse_args()
	fam = Family(args.pop_size, args.num_runs, args.verbose)
