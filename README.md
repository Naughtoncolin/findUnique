# findUnique
This was originally written Spring 2016 as one of the last assignments for BME160: Research Programming for the Life Sciences with Dr. David Bernick.
Prior to assigning it the class, Dr. Bernick made the algorithm to find unique tags within mitochondrial tRNAs to discrimintate between the diffrent types when sequencing pool cia nanopore. As such, the example test file with it's associated results file are mitochondrial tRNA sequences.
Compare sequence pools to find unique subsequences within each.
Reads a file of FASTA sequences from STDIN, and in each sequence, finds the unique subsequences that don't occur among any of the other FASTA sequences. Each of the sets of subsequences is minimized, such that no member of a sequences' subsequence set set is found as a substring of any other member of that set.
