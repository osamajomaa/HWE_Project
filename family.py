import csv
import sys
import scipy
import random
# from sets import set
from collections import OrderedDict
from scipy.stats.stats import chisquare
import scipy.stats.mstats as mst
from scipy import stats
import numpy as np


class Family:

	# dna would be a 2d matrix (row = person, col = loci halfs, (row, col) = allel name)
	# csv = all possible people
	# parent1 = random row
	# parent2 = random row
	def __init__(self, dna):
		self.csv = dna
		self.parent1 = self.getRandomPerson()
		self.parent2 = self.getRandomPerson()
	
	# Gets a random row
	# The elements of the rows are the name's of the allel
	def getRandomPerson(self):
		rand = random.randrange(1, len(self.csv))
		return self.csv[rand]
	
	# parent1 and parent2 are arrays.
	# They would be a row of data as in the csv files.
	# the return is an array randomly created from the parents
	def generateChild(self, parent1, parent2):
		child = []
		#use child.append()
		i = 0
		end = len(parent1)
		while i < end:
			rand = random.choice([0,1])
			if(rand == 0):
				child.append(parent1[i])
			elif(rand == 1):
				child.append(parent2[i])
			else:
				print("There was an error choosing an allel.")
			i = i + 1
		return child
	
	# the method reurns 2 children
	def getSiblings(self, parent1, parent2):
		return [generateChild(parent1, parent2), generateChild(parent1, parent2)]
	
	# Creates a child and its cusion
	# the two parameters are grand parents
	def getCousins(self, gParent1, gParent2):
		# makes a parent and aunt/uncle
		s = getSiblings(gParent1, gParent2)
		# first child's other parent
		r1 = getRandomPerson()
		# the first child
		c1 = generateChild(s[0], r1)
		# the cusion's other parent
		r2 = getRandomPerson()
		# create the cusion
		c2 = generateChild(s[1], r2)
		return (c1, c2)
		
	def percentMatch(self, person1, person2):
		numLoci = len(person1) / 2
		i = 0
		numMatches = 0
		while i < numLoci:
			if(person1[i] == person2[i] and person1[i + 1] == person2[i + 1]):
				numMatches = numMatches + 1
			i = i + 2
		return (numMatches / numLoci)
		
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
				significantMatches = significantMatches +1
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
		
