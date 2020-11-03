# findUnique
Reads a file of FASTA sequences from STDIN, and in each sequence, finds the unique subsequences that don't occur among any of the other FASTA sequences. Each of the sets of subsequences is minimized, such that no member of a sequences' subsequence set is found as a substring of any other member of that set.

INPUT: FASTA file
OUTPUT:
Prints report with the following
Line 1: the sequence name
Line 2: the sequence
Lines 3-xx: the unique elements of the parent sequences ordered by their starting position, with each element on its own line, and spaced so they are directly under that subsequence in the original. Looks similar to an allignment.
This then repeats for each of the FASTA sequences. An example of the output can be seen in Results.txt


This was originally written Spring 2016 as one of the last assignments for BME160: Research Programming for the Life Sciences with Dr. David Bernick.
Prior to assigning it the class, Dr. Bernick made the algorithm to find unique tags within mitochondrial tRNAs to discrimintate between the diffrent types when sequencing pool cia nanopore. As such, the example test file with it's associated results file are mitochondrial tRNA sequences.
