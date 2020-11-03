#!/usr/bin/env python3
# Author: Colin Naughton
# This code was originally written Spring 2016 for BME160: Research Programming in the Life Sciences with Dr. David Bernick at UCSC.
'''
PSEUDOCODE:
Import sequenceAnalysis
Read FASTA file
	Assign FASTA headers and corresponding sequence objects to a list.
	Iterate through header list.
		Remove substrings found in other tRNA sequences within the FASTA file.
		Remove substrings which contained other substrings from the tRNA.
		Print tRNA header.
		Print tRNA sequence.
		Print unique sequences from the tRNA at its position in the sequence.
'''
import sequenceAnalysis

class tRNAseq:
	'''
	Models a strand of tRNA.

	Keyword arguements:
	strand = Sequence of the tRNA strand.
	otherSet = Set of substrings to be removed the current tRNA's powerset. 
	'''

	def __init__(self, strand):
		'''Initialize the sequence and its substring powerset.'''
		self.seq = strand
		self.seq = ''.join([x for x in self.seq if x!='-' if x!='_']) #Remove unwanted characters from the sequence.
		self.subStrings = set(self.seq[i:k] for i in range(len(self.seq)) for k in range(len(self.seq), i, -1)) # Caluculate substring powerset.
		self.uniqueSet = set()

	def Sequence(self):
		'''Returns sequence.'''
		return self.seq
	
	def powerSet(self):
		'''Returns powerset of sequence's substrings.'''
		return self.subStrings

	def uniqueSet(self, otherSet):
		'''Calculates and returns unique substrings from sequence's substring powerset when compared to another set. '''
		self.uniqueSet = self.subStrings - otherSet #
		return self.uniqueSet


'''Create object of FASTA headers and their complimentary sequences from FASTA file read from directory.'''
thing = sequenceAnalysis.FastAreader()

'''Iterate through FASTA header(s) and sequence(s), making a list of headers and their corresponding tRNA objects''' 
tRNAdict = list()
for head, seq in thing.readFasta() :
	tRNAdict.append((head, tRNAseq(seq))) #Create dictionary value for a tRNAfrom set of substrings.

'''Find and subtract union of other tRNAs from the current tRNA.'''
for tRNAheader in sorted(tRNAdict, key=lambda x: x[0]): #Iterate through dictionary of tRNA and their substrings to be subtracted from.
	subStr = tRNAheader[1].powerSet()
	seq = tRNAheader[1].Sequence()
	unionSet = set() #Create empty set to contain the union of all OTHER substrings than the one to be subtracted from.
	for tRNA in tRNAdict: #Iterate through dictionary of substrings, find all OTHER substrings.
			if tRNA[0] != tRNAheader[0]:   
				unionSet = unionSet.union(tRNA[1].powerSet()) #Add to the set of all other substrings, if substring is not already present.
	union = sorted(list(subStr - unionSet), key=len)
	newList = set(j for i in range(len(union)) for j in union[i+1:len(union)] if union[i] in j)
	union=set(union)
	ham = list(union-newList)

	ham.sort(key=lambda x: seq.index(x))
	print(tRNAheader[0]) #Print the FASTA header for the tRNA.
	print(seq) #Print the sequence of the tRNA.
	'''Print unique sequences of tRNA is position they would be in the original sequence.'''
	for greenEggs in ham:
		print('.'*(seq.index(greenEggs)),end='')
		print(greenEggs)