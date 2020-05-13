#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 20:51:30 2020

@author: bs15ansj

This is a small script to help you get used to using biopython. In order to
perform the clustalw multiple sequence alignment, make sure you have clustal
installed. You can do this using the following commands:
    
    sudo apt-get update
    sudo apt-get install clustalw
    
Once installed, make sure that you close the terminal and reopen before 
running the script. 

"""

from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

aas =['V','I','L','E','Q', \
'D','N','H','W','F','Y', \
'R','K','S','T','M','A', \
'G','P','C','-']

def run_blast():
    
    # Read the input file and get the sequence
    
    with open("input.txt", "r") as f:
        in_file = f.read().strip()
        
    # Blast the sequence over the internet (you may have to wait) and 
    # save the results in a qblast object
    
    result_handle = NCBIWWW.qblast("blastp", "nr", 
                                   in_file, hitlist_size=200)
    
    # Open the output file then read the results object and write to 
    # the output file
    
    with open("output.xml", "w") as out_handle:
        out_handle.write(result_handle.read())

def write_sequences():
    
    # Load the blast output file
    
    blast_qresult = SearchIO.read("output.xml", 
                                  "blast-xml")

    # Iterate through ids and sequences and add them to lists. Sequences are
    # only added to the list if they are not already in the list. This is done
    # because sometimes there will be duplicates in the blast output, which 
    # will produce an error when aligning with clustal 
    
    ids = []
    sequences = []
    for hsp in blast_qresult.hsps:
        if hsp.hit.id not in ids:
            ids.append(str(hsp.hit.id))
            sequences.append(str(hsp.hit.seq))
            
    # Open the sequences output file then for each high-scoring pair
    # in the blast results, write the hit ID (proceeded by a ">" for 
    # fasta format), followed by the hit sequence on the next line    
    with open("sequences.fasta", "w") as f:
        for i, s in zip(ids, sequences):
            f.write('> '+i+'\n')
            f.write(s+'\n')

def clustal_alignment():
    
    # Run the command line version of clustal using the sequences.fasta
    # file and output to a clustal format alignment file
    
    cmd = ClustalwCommandline("clustalw", 
                              infile="sequences.fasta", 
                              outfile="alignment.aln")
    cmd()
    
def MSA_matrix():
    
    # Read the clustal alignment
    align = AlignIO.read("alignment.aln",'clustal')
    
    # Read sequences, convert each to list and add to matrix
    matrix = []
    for record in align:
        matrix.append(list(record.seq))
    matrix = np.array(matrix)
    return matrix     
    
def frequency():

    # Generate the MSA matrix and create a new empty frequency matrix
    
    M = MSA_matrix()
    F = []
    
    # Set the number of sequences equal to the number of rows in the matrix
    # (using .shape on a matrix returns a two integer tuple, the first number 
    # is the number of rows, the second is the number of colums)
    
    no_sequences = M.shape[0]
    
    # Enumerate through the transpose of the MSA matrix. 
    # The transpose of a matrix flips the matrix along the diaganol axis e.g.
    # 
    #  A B C               A D G    
    #  D E F  transpose => B E H    
    #  G H I               C F I
    # 
    # This is because when we enumerate, we go through the matrix row by row. 
    # We want to enumerate through the columns so we transpose the matrix, 
    # effectively switching the rows and colums. 
    
    for index, i in enumerate(M.T):
        
        # Create an empty list that will store the frequency of each amino
        # acid for position i
        
        f = []
        
        # Iterate over the amino acids and count their frequency in position i
        # and add this to the list
        
        for aa in aas:
            x = np.sum(np.char.count(i, aa)) / no_sequences
            f.append(x)
        
        # After creating the list of frequencies, f, add this to the frequency
        # matrix, F, for each position i
            
        F.append(f)
    
    # Return the frequency matrix as an array  
        
    return np.array(F)
            
def consensus_seq():
    
    # Generate a frequency matrix
    F = frequency()
    
    # Create empty lists for the consensus sequence, and the frequency of 
    # consensus residues
    con_seq = []
    con_freq = []
    
    # Iterate over rows of the matrix (remember, this has been transposed so 
    # this is equivalent to iterating over columns of the MSA)
    for i in F:
        
        # Find column where the frequency is maximum. This column index
        # corresponds to the list index for a particular amino acid in aas at 
        # the top of this script. 
        m = np.where(i == np.amax(i))
        
        # Use the column index in m to find the amino acid letter in aas, then
        # add this letter to the consensus sequence
        con_seq.append(aas[int(m[0][0])])
        
        # Add the max frequency to the consensus frequency list
        con_freq.append(np.amax(i))

    return con_seq, con_freq

def plot_consensus_frequency_versus_position():
    
    # Generate the consensus frequency list (consensus_seq() returns two
    # variable, so we need to create a variable name for both even though
    # we will only need con_freq)
    con_seq, con_freq = consensus_seq()
    
    plt.plot(range(1, len(con_freq)+1), con_freq)
    plt.xlim(1, len(con_freq))
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.savefig('consensus_frequency_versus_position.png', dpi = 300)


write_sequences()
clustal_alignment()
print(consensus_seq())



